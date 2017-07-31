#include "DICE.h"
#include <ncDim.h>
#include <ncFile.h>
#include <ncVar.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include "DICEClimate.h"
#include "DICEDamage.h"
#include "csv-parser.h"
#include "settingsnode.h"

#define WITH_NLOPT
//#define WITH_PAGMO
//#define WITH_BORG

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wswitch-enum"
#pragma GCC diagnostic ignored "-Wfloat-equal"
#ifdef WITH_NLOPT
#include <nlopt.hpp>
#endif
#ifdef WITH_PAGMO
#include <pagmo/pagmo.hpp>
#endif
#ifdef WITH_BORG
#include "borg.h"
#endif
#pragma GCC diagnostic pop

namespace dice {

template<typename Value, typename Time>
DICE<Value, Time>::DICE(const settings::SettingsNode& settings_p)
    : settings(settings_p), global(settings["parameters"]), control(global.timestep_num), emissions(global, economies) {
}

template<typename Value, typename Time>
void DICE<Value, Time>::initialize() {
    // Initialize climate module
    {
        const settings::SettingsNode& climate_node = settings["climate"];
        const std::string& type = climate_node["type"].as<std::string>();
        if (type == "dice") {
            climate.reset(new climate::DICEClimate<Value, Time>(climate_node["parameters"], global, emissions));
        } else {
            throw std::runtime_error("unknown climate module type '" + type + "'");
        }
        climate->initialize();
    }

    // Initialize damage module
    {
        const settings::SettingsNode& damage_node = settings["damage"];
        const std::string& type = damage_node["type"].as<std::string>();
        if (type == "dice") {
            damage.reset(new damage::DICEDamage<Value, Time>(damage_node["parameters"], global, *climate));
        } else {
            throw std::runtime_error("unknown damage module type '" + type + "'");
        }
        damage->initialize();
    }

    // Initialize regions
    {
        for (const auto&& region_node : settings["regions"].as_sequence()) {
            economies.emplace_back(Economy<Value, Time>(region_node["economy"], global, control, *climate, *damage));
        }
    }

    // Initialize control variables
    if (settings.has("control")) {
        class ControlInputObserver : public Observer<Value, Time> {
          protected:
            const settings::SettingsNode& input_node;

          public:
            ControlInputObserver(const settings::SettingsNode& input_node_p) : input_node(input_node_p){};
            bool observe(const std::string& name, TimeSeries<Value>& v) override {
                if (input_node.has(name)) {
                    const settings::SettingsNode& node = input_node[name];
                    const std::string& format = node["format"].as<std::string>();
                    if (format == "csv") {
                        try {
                            const std::string& filename = node["filename"].as<std::string>();
                            std::ifstream datastream(filename);
                            if (!datastream) {
                                throw std::runtime_error("could not open '" + filename + "'");
                            }
                            csv::Parser parser(datastream);
                            const size_t col = node["column"].as<size_t>();
                            parser.next_row();  // Skip header row
                            for (auto d = v.begin(); d != v.end(); ++d) {
                                for (size_t c = 0; c < col; c++) {
                                    parser.next_col();
                                }
                                *d = parser.read<Value>();
                                parser.next_row();
                            }
                        } catch (const csv::parser_exception& ex) {
                            std::stringstream s;
                            s << ex.what();
                            s << " (line " << ex.row << " col " << ex.col << ")";
                            throw std::runtime_error(s.str());
                        }
                    } else {
                        throw std::runtime_error("unknown format '" + format + "'");
                    }
                }
                return true;
            }
        };
        const settings::SettingsNode& input_node = settings["control"];
        ControlInputObserver observer(input_node);
        control.observe(observer);
    }
}

static const char* get_optimization_results(const int& result_) {
    switch (result_) {
        case 1:
            return "Generic success";
        case 2:
            return "Optimization reached target objective value";  // "Optimization stopped because stopval was reached"
        case 3:
            return "Optimization reached target objective precision";  // "Optimization stopped because ftol_rel or ftol_abs was reached"
        case 4:
            return "Optimization reached target control variable precision";  // "Optimization stopped because xtol_rel or xtol_abs was reached"
        case 5:
            return "Optimization maximum iterations reached";  // "Optimization stopped because maxeval was reached";
        case 6:
            return "Optimization timed out";  // "Optimization stopped because maxtime was reached"
        case -1:
            return "Generic failure code";
        case -2:
            return "Invalid arguments(e.g. lower bounds are bigger than upper bounds, an unknown algorithm was specified, etc.";
        case -3:
            return "Ran out of memory";
        case -4:
            return "Halted because roundoff errors limited progress.(In this case, the optimization still typically returns a useful result.) ";
        case -5:
            return "Halted because of a forced termination: the user called nlopt_force_stop(opt) on the optimization’s nlopt_opt object opt from the user’s "
                   "objective function or constraints.";
        default:
            return "Unknown optimization result";
    }
}

#ifdef WITH_BORG
static DICE<double, size_t>* dice;
static size_t n;
static bool with_constraint;
#endif

template<typename Value, typename Time>
void DICE<Value, Time>::optimize(const settings::SettingsNode& optimization_node, TimeSeries<Value>& initial_values) {
    const std::string& library = optimization_node["library"].as<std::string>();
    const Time optimization_variables = initial_values.size();
    if (false) {  // only for ifdefs below
#ifdef WITH_PAGMO
    } else if (library == "pagmo") {
        std::fill(std::begin(control.s), std::end(control.s), global.optlrsav);
        struct PagmoProblem {
            DICE<Value, Time>* dice;
            Time n;
            bool with_constraint;
            pagmo::vector_double fitness(const pagmo::vector_double& vars) const {
                pagmo::vector_double f(with_constraint ? 2 : 1);
                dice->control.s.assign(std::begin(vars), std::end(vars));
                dice->reset();
                f[0] = -dice->calc_single_utility();
                if (with_constraint) {
                    Value cca = 90;  // TODO
                    for (Time t = 0; t < dice->global.timestep_num; ++t) {
                        cca += dice->global.timestep_length * dice->emissions(t) / 3.666;
                    }

                    f[1] = cca - dice->global.fosslim;
                }
                return f;
            };
            std::pair<pagmo::vector_double, pagmo::vector_double> get_bounds() const {
                pagmo::vector_double lb(n, 0);
                pagmo::vector_double ub(n, 1);
                return {lb, ub};
            };
            pagmo::vector_double::size_type get_nobj() const {
                return 1;
            };
            pagmo::vector_double::size_type get_nec() const {
                return 0;
            };
            pagmo::vector_double::size_type get_nic() const {
                return with_constraint ? 1 : 0;
            };
            pagmo::thread_safety get_thread_safety() const {
                return pagmo::thread_safety::none;
            };
        };
        PagmoProblem pagmo_problem;
        pagmo_problem.dice = this;
        pagmo_problem.n = optimization_variables;
        pagmo_problem.with_constraint = optimization_node["limit_cca"].as<bool>();
        pagmo::problem problem{pagmo_problem};
        pagmo::population population{problem};
        population.push_back(initial_values);
        pagmo::algorithm algorithm;

        const std::string& solver_name = optimization_node["solver"].as<std::string>();
        if (solver_name == "ipopt") {
            pagmo::ipopt solver;
            solver.set_numeric_option("tol", optimization_node["utility_precision"].as<Value>() / 3000);
            if (optimization_node.has("maxiter")) {
                solver.set_integer_option("max_iter", optimization_node["maxiter"].as<size_t>());
            }
            if (optimization_node.has("timeout")) {  // timeout given in sec
                solver.set_integer_option("max_cpu_time", optimization_node["timeout"].as<size_t>());
            }
            solver.set_selection("best");
            algorithm = pagmo::algorithm{solver};
        } else if (solver_name == "nlopt") {
            pagmo::nlopt solver(optimization_node["algorithm"].as<std::string>());
            solver.set_ftol_abs(optimization_node["utility_precision"].as<Value>());
            if (optimization_node.has("maxiter")) {
                solver.set_maxeval(optimization_node["maxiter"].as<size_t>());
            }
            if (optimization_node.has("timeout")) {  // timeout given in sec
                solver.set_maxtime(optimization_node["timeout"].as<size_t>());
            }
            algorithm = pagmo::algorithm{solver};
        } else {
            throw std::runtime_error("unknown solver '" + solver_name + "'");
        }
        for (size_t i = 0; i < optimization_node["iterations"].as<size_t>(1); ++i) {
            population = algorithm.evolve(population);
        }
        pagmo::vector_double vars = population.champion_x();
        control.s.assign(std::begin(vars), std::end(vars));
#endif
#ifdef WITH_BORG
    } else if (library == "borg") {
        dice = this;
        n = optimization_variables;
        with_constraint = optimization_node["limit_cca"].as<bool>();
        std::fill(std::begin(control.s), std::end(control.s), global.optlrsav);
        BORG_Problem opt = BORG_Problem_create(optimization_variables, 1, with_constraint ? 1 : 0, [](double* vars, double* objs, double* consts) {
            dice->control.s.assign(vars, vars + n);
            dice->reset();
            Value utility = dice->calc_single_utility();
            Value cca = 90;  // TODO
            for (Time t = 0; t < dice->global.timestep_num; ++t) {
                cca += dice->global.timestep_length * dice->emissions(t) / 3.666;
            }
            objs[0] = -utility;
            if (with_constraint) {
                consts[0] = std::max(0.0, cca - dice->global.fosslim);
            }
        });

        for (Time t = 0; t < optimization_variables; ++t) {
            BORG_Problem_set_bounds(opt, t, 0, 1);
        }

        BORG_Problem_set_epsilon(opt, 0, optimization_node["utility_precision"].as<Value>());
        // BORG_Random_seed(12345);
        BORG_Archive result = BORG_Algorithm_run(opt, optimization_node["maxiter"].as<size_t>());
        // BORG_Archive_print(result, stdout);
        BORG_Archive_destroy(result);
        BORG_Problem_destroy(opt);
#endif
#ifdef WITH_NLOPT
    } else if (library == "nlopt") {
        // "auglag nlopt::algorithm::LN_AUGLAG
        // "auglag nlopt::algorithm::LD_AUGLAG
        // "auglag_eq nlopt::algorithm::LN_AUGLAG_EQ
        // "auglag_eq nlopt::algorithm::LD_AUGLAG_EQ
        const std::string algorithm_name = optimization_node["algorithm"].as<std::string>();
        nlopt::algorithm algorithm_type;
        if (algorithm_name == "direct") {
            algorithm_type = nlopt::algorithm::GN_DIRECT;
        } else if (algorithm_name == "direct_l") {
            algorithm_type = nlopt::algorithm::GN_DIRECT_L;
        } else if (algorithm_name == "direct_lrand") {
            algorithm_type = nlopt::algorithm::GN_DIRECT_L_RAND;
        } else if (algorithm_name == "direct_noscal") {
            algorithm_type = nlopt::algorithm::GN_DIRECT_NOSCAL;
        } else if (algorithm_name == "direct_lnoscal") {
            algorithm_type = nlopt::algorithm::GN_DIRECT_L_NOSCAL;
        } else if (algorithm_name == "direct_lrand_noscal") {
            algorithm_type = nlopt::algorithm::GN_DIRECT_L_RAND_NOSCAL;
        } else if (algorithm_name == "orig_direct") {
            algorithm_type = nlopt::algorithm::GN_ORIG_DIRECT;
        } else if (algorithm_name == "orig_direct_l") {
            algorithm_type = nlopt::algorithm::GN_ORIG_DIRECT_L;
        } else if (algorithm_name == "stogo") {
            algorithm_type = nlopt::algorithm::GD_STOGO;
        } else if (algorithm_name == "stogo_rand") {
            algorithm_type = nlopt::algorithm::GD_STOGO_RAND;
        } else if (algorithm_name == "lbfgs-nocedal") {
            algorithm_type = nlopt::algorithm::LD_LBFGS_NOCEDAL;
        } else if (algorithm_name == "lbfgs") {
            algorithm_type = nlopt::algorithm::LD_LBFGS;
        } else if (algorithm_name == "praxis") {
            algorithm_type = nlopt::algorithm::LN_PRAXIS;
        } else if (algorithm_name == "var1") {
            algorithm_type = nlopt::algorithm::LD_VAR1;
        } else if (algorithm_name == "var2") {
            algorithm_type = nlopt::algorithm::LD_VAR2;
        } else if (algorithm_name == "tnewton") {
            algorithm_type = nlopt::algorithm::LD_TNEWTON;
        } else if (algorithm_name == "tnewton_restart") {
            algorithm_type = nlopt::algorithm::LD_TNEWTON_RESTART;
        } else if (algorithm_name == "tnewton_precond") {
            algorithm_type = nlopt::algorithm::LD_TNEWTON_PRECOND;
        } else if (algorithm_name == "tnewton_precond_restart") {
            algorithm_type = nlopt::algorithm::LD_TNEWTON_PRECOND_RESTART;
        } else if (algorithm_name == "crs2_lm") {
            algorithm_type = nlopt::algorithm::GN_CRS2_LM;
        } else if (algorithm_name == "gn_mlsl") {
            algorithm_type = nlopt::algorithm::GN_MLSL;
        } else if (algorithm_name == "mlsl") {
            algorithm_type = nlopt::algorithm::GD_MLSL;
        } else if (algorithm_name == "mlsl_lds") {
            algorithm_type = nlopt::algorithm::GN_MLSL_LDS;
        } else if (algorithm_name == "mlsl_lds") {
            algorithm_type = nlopt::algorithm::GD_MLSL_LDS;
        } else if (algorithm_name == "mma") {
            algorithm_type = nlopt::algorithm::LD_MMA;
        } else if (algorithm_name == "cobyla") {
            algorithm_type = nlopt::algorithm::LN_COBYLA;
        } else if (algorithm_name == "newuoa") {
            algorithm_type = nlopt::algorithm::LN_NEWUOA;
        } else if (algorithm_name == "newuoa_bound") {
            algorithm_type = nlopt::algorithm::LN_NEWUOA_BOUND;
        } else if (algorithm_name == "neldermead") {
            algorithm_type = nlopt::algorithm::LN_NELDERMEAD;
        } else if (algorithm_name == "sbplx") {
            algorithm_type = nlopt::algorithm::LN_SBPLX;
        } else if (algorithm_name == "bobyqa") {
            algorithm_type = nlopt::algorithm::LN_BOBYQA;
        } else if (algorithm_name == "isres") {
            algorithm_type = nlopt::algorithm::GN_ISRES;
        } else if (algorithm_name == "slsqp") {
            algorithm_type = nlopt::algorithm::LD_SLSQP;
        } else if (algorithm_name == "ccsaq") {
            algorithm_type = nlopt::algorithm::LD_CCSAQ;
        } else if (algorithm_name == "esch") {
            algorithm_type = nlopt::algorithm::GN_ESCH;
        } else {
            throw std::runtime_error("unknown algorithm '" + algorithm_name + "'");
        }
        /*
        nlopt::opt super_opt(nlopt::LN_AUGLAG, optimization_variables);
        super_opt.add_inequality_constraint(
            [](unsigned n, const double* x, double* grad, void* data) {
                DICE* dice = static_cast<DICE*>(data);
                dice->control.s.assign(x, x + n);
                dice->reset();
                dice->calc_single_utility();
                Value cca = 90;  // TODO
                for (Time t = 0; t < dice->global.timestep_num - 1; ++t) {
                    cca += dice->global.timestep_length * dice->emissions(t) / 3.666;
                }
                return cca - dice->global.fosslim;
            },
            this, 0.1);
        super_opt.set_max_objective(
            [](unsigned n, const double* x, double* grad, void* data) {
                DICE* dice = static_cast<DICE*>(data);
                dice->control.s.assign(x, x + n);
                dice->reset();
                return dice->calc_single_utility();
            },
            this);
        super_opt.set_ftol_abs(optimization_node["utility_precision"].as<Value>());
        super_opt.set_lower_bounds(std::vector<Value>(optimization_variables, 0));
        super_opt.set_upper_bounds(std::vector<Value>(optimization_variables, 1));
        if (optimization_node.has("maxiter")) {
            super_opt.set_maxeval(optimization_node["maxiter"].as<size_t>());
        }
        if (optimization_node.has("timeout")) {  // timeout given in sec
            super_opt.set_maxtime(optimization_node["timeout"].as<size_t>());
        }
        //*/

        //*
        nlopt::opt opt(algorithm_type, optimization_variables);
        opt.set_max_objective(
            [](unsigned n, const double* x, double* grad, void* data) {
                DICE* dice = static_cast<DICE*>(data);
                dice->control.s.assign(x, x + n);
                dice->reset();
                return dice->calc_single_utility();
            },
            this);
        //*/

        /*
        nlopt::opt opt(nlopt::LN_BOBYQA, optimization_variables);
        opt.set_max_objective(
            [](unsigned n, const double* x, double* grad, void* data) {
                DICE* dice = static_cast<DICE*>(data);
                dice->control.s.assign(x, x + n);
                dice->reset();
                Value utility = dice->calc_single_utility();
                Value cca = 90;  // TODO
                for (Time t = 0; t < dice->global.timestep_num; ++t) {
                    cca += dice->global.timestep_length * dice->emissions(t) / 3.666;
                }
                return utility - 1e5 * std::max(0.0, cca - dice->global.fosslim);
            },
            this);
        //*/
        opt.set_ftol_abs(optimization_node["utility_precision"].as<Value>());
        opt.set_lower_bounds(std::vector<Value>(optimization_variables, 0));
        opt.set_upper_bounds(std::vector<Value>(optimization_variables, 1));
        if (optimization_node.has("maxiter")) {
            opt.set_maxeval(optimization_node["maxiter"].as<size_t>());
        }
        if (optimization_node.has("timeout")) {  // timeout given in sec
            opt.set_maxtime(optimization_node["timeout"].as<size_t>());
        }
        // super_opt.set_local_optimizer(opt);
        // nlopt::result result = super_opt.optimize(initial_values, utility);
        Value utility;
        nlopt::result result = opt.optimize(initial_values, utility);
        std::cout << get_optimization_results(result) << std::endl;
#endif
    } else {
        throw std::runtime_error("unknown library '" + library + "'");
    }
}

template<typename Value, typename Time>
void DICE<Value, Time>::run() {
    if (economies.size() == 0) {
        throw std::runtime_error("no economies given");
    }
    if (economies.size() == 1) {
        // Control rate limits
        // MIU.up[t] = limmu * partfract[t];
        // MIU.up[t] $(t.val < 146) = 1;

        Value utility;

        const settings::SettingsNode& optimization_node = settings["optimization"];
        if (optimization_node.has("iterations")) {
            const Time optimization_variables = global.timestep_num - optimization_node["s_fix_steps"].as<Time>(0);
            std::fill(std::begin(control.s), std::end(control.s), global.optlrsav);
            TimeSeries<Value> initial_values(optimization_variables, 0);
            for (const auto& iteration_node : optimization_node["iterations"].as_sequence()) {
                initial_values.assign(std::begin(control.s), std::begin(control.s) + optimization_variables);
                optimize(iteration_node, initial_values);
                reset();
                utility = calc_single_utility();
            }
        } else {
            utility = calc_single_utility();
        }
        std::cout << "Finished with utility = " << utility << std::endl;
    } else {
        throw std::runtime_error("multiple regions not supported yet");
    }
}

template<typename Value, typename Time>
void DICE<Value, Time>::reset() {
    emissions.reset();
    climate->reset();
    damage->reset();
    for (auto&& economy : economies) {
        economy.reset();
    }
}

template<typename Value, typename Time>
Value DICE<Value, Time>::calc_single_utility() {
#ifdef DEBUG
    try {
#endif
        Value utility = 0;
        for (Time t = 0; t < global.timestep_num; ++t) {
            utility += economies[0].utility(t);
        }
        return global.scale1 * utility + global.scale2;
#ifdef DEBUG
    } catch (std::exception& e) {
        std::cerr << "Exception '" << e.what() << "' in optimization" << std::endl;
        throw;
    }
#endif
}

template<typename Value, typename Time>
void DICE<Value, Time>::output() {
    if (settings.has("output")) {
        const settings::SettingsNode& output_node = settings["output"];
        const std::string& type = output_node["type"].as<std::string>();
        if (type == "netcdf") {
            write_netcdf_output(output_node);
        }
    }
}

template<typename Value, typename Time>
void DICE<Value, Time>::write_netcdf_output(const settings::SettingsNode& output_node) {
    if (economies.size() == 1) {
        netCDF::NcFile file(output_node["filename"].as<std::string>(), netCDF::NcFile::replace, netCDF::NcFile::nc4);

        netCDF::NcDim time_dim = file.addDim("time", global.timestep_num);
        netCDF::NcVar time_var = file.addVar("time", netCDF::NcType::nc_UINT, {time_dim});
        for (Time t = 0; t < global.timestep_num; ++t) {
            const Time year = global.start_year + t * global.timestep_length;
            time_var.putVar({t}, (const unsigned int)year);
        }
        class NetCDFOutputObserver : public Observer<Value, Time> {
          protected:
            const netCDF::NcGroup& group;
            const netCDF::NcDim& time_dim;
            const settings::SettingsNode& output_node;

          public:
            NetCDFOutputObserver(const netCDF::NcGroup& group_p, const netCDF::NcDim& time_dim_p, const settings::SettingsNode& output_node_p)
                : group(group_p), time_dim(time_dim_p), output_node(output_node_p){};
            bool observe(const std::string& name, TimeSeries<Value>& v) override {
                netCDF::NcVar var = group.addVar(name, netCDF::NcType::nc_FLOAT, {time_dim});
                var.setCompression(false, true, 7);
                // var.setFill<Value>(true, std::numeric_limits<Value>::quiet_NaN());
                var.putVar(&v[0]);
                return true;
            }
        };
        NetCDFOutputObserver observer(file, time_dim, output_node);
        economies[0].observe(observer);
        climate->observe(observer);
        damage->observe(observer);
        control.observe(observer);
        file.putAtt("utility", netCDF::NcType::nc_FLOAT, calc_single_utility());
    } else {
        throw std::runtime_error("multiple regions not supported yet");
    }
}

template class DICE<double, size_t>;
}
