/*
  Copyright (C) 2017 Sven Willner <sven.willner@gmail.com>

  This file is part of DICE++.

  DICE++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as
  published by the Free Software Foundation, either version 3 of
  the License, or (at your option) any later version.

  DICE++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with DICE++.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "Optimization.h"
#include "DICE.h"
#include "config.h"
#include "settingsnode.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wswitch-enum"
#pragma GCC diagnostic ignored "-Wfloat-equal"
#ifdef DICEPP_WITH_NLOPT
#include <nlopt.hpp>
#endif
#ifdef DICEPP_WITH_PAGMO
#include <pagmo/pagmo.hpp>
#endif
#ifdef DICEPP_WITH_BORG
#include "borg.h"
#endif
#ifdef DICEPP_WITH_MIDACO
extern "C" {
int midaco(long int*,
           long int*,
           long int*,
           long int*,
           long int*,
           long int*,
           double*,
           double*,
           double*,
           double*,
           double*,
           long int*,
           long int*,
           double*,
           double*,
           long int*,
           long int*,
           long int*,
           double*,
           long int*,
           char*);
int midaco_print(int,
                 long int,
                 long int,
                 long int*,
                 long int*,
                 double*,
                 double*,
                 double*,
                 double*,
                 double*,
                 long int,
                 long int,
                 long int,
                 long int,
                 long int,
                 double*,
                 double*,
                 long int,
                 long int,
                 double*,
                 long int,
                 char*);
}
#endif
#pragma GCC diagnostic pop

namespace dice {

#ifdef DICEPP_WITH_NLOPT
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
#endif

#ifdef DICEPP_WITH_BORG
static Optimization<double, size_t>* optimization;
#endif

template<typename Value, typename Time>
void Optimization<Value, Time>::optimize(const settings::SettingsNode& settings, TimeSeries<Value>& initial_values) {
    const std::string& library = settings["library"].as<std::string>();
    if (library == "midaco") {
#ifdef DICEPP_WITH_MIDACO
        long int o, n, ni, m, me, maxeval, maxtime, printeval, save2file, iflag = 0, istop = 0;
        std::vector<double> x(initial_values), xl(variables_num, 0), xu(variables_num, 1), param(12);
        char key[] = "MIDACO_LIMITED_VERSION___[CREATIVE_COMMONS_BY-NC-ND_LICENSE]";

        o = 1;              // Number of objectives
        n = variables_num;  // Number of variables (in total)
        ni = 0;             // Number of integer variables (0 <= ni <= n)
        m = 1;              // Number of constraints (in total)
        me = 0;             // Number of equality constraints (0 <= me <= m)

        long int p = 0;  // parallelization
        long int lrw = 105 * n + m * p + 2 * m + o * o + 4 * o * p + 10 * o + 3 * p + 610;
        std::vector<double> rw(lrw);
        long int liw = 3 * n + p + 110;
        std::vector<long int> iw(liw);
        long int paretomax = 100;
        long int lpf = (o + m + n) * paretomax + 1;
        std::vector<double> pf(lpf);
        x = initial_values;

        printeval = 1000;  // Print-Frequency for current best solution (e.g. 1000)
        save2file = 0;     // Save SCREEN and SOLUTION to TXT-files [ 0=NO/ 1=YES]

        maxeval = 10000;
        maxtime = 60;

        param[0] = 0.0;   // ACCURACY
        param[1] = 0.0;   // SEED
        param[2] = 0.0;   // FSTOP
        param[3] = 0.0;   // ALGOSTOP
        param[4] = 0.0;   // EVALSTOP
        param[5] = 0.0;   // FOCUS
        param[6] = 0.0;   // ANTS
        param[7] = 0.0;   // KERNEL
        param[8] = 0.0;   // ORACLE
        param[9] = 0.0;   // PARETOMAX
        param[10] = 0.0;  // EPSILON
        param[11] = 0.0;  // CHARACTER

        Value f = 0;
        Value g = 0;
        midaco_print(1, printeval, save2file, &iflag, &istop, &f, &g, &x[0], &xl[0], &xu[0], o, n, ni, m, me, &rw[0], &pf[0], maxeval, maxtime, &param[0], p,
                     key);
        while (istop == 0) {
            f = -objective(&x[0], nullptr)[0];   // TODO
            g = -constraint(&x[0], nullptr)[0];  // TODO
            midaco(&p, &o, &n, &ni, &m, &me, &x[0], &f, &g, &xl[0], &xu[0], &iflag, &istop, &param[0], &rw[0], &lrw, &iw[0], &liw, &pf[0], &lpf, key);
            midaco_print(2, printeval, save2file, &iflag, &istop, &f, &g, &x[0], &xl[0], &xu[0], o, n, ni, m, me, &rw[0], &pf[0], maxeval, maxtime, &param[0],
                         p, key);
        }
#else
        throw std::runtime_error("library '" + library + "' not supported by this binary");
#endif
    } else if (library == "pagmo") {
#ifdef DICEPP_WITH_PAGMO
        struct PagmoProblem {
            Optimization<Value, Time>* optimization;
            pagmo::vector_double fitness(const pagmo::vector_double& vars) const {
                pagmo::vector_double f(optimization->objectives_num + optimization->constraints_num);
                f[0] = -optimization->objective(&vars[0], nullptr)[0];  // TODO
                if (optimization->constraints_num > 0) {
                    f[1] = optimization->constraint(&vars[0], nullptr)[0];  // TODO
                }
                return f;
            }
            std::pair<pagmo::vector_double, pagmo::vector_double> get_bounds() const {
                pagmo::vector_double lb(optimization->variables_num, 0);
                pagmo::vector_double ub(optimization->variables_num, 1);
                return {lb, ub};
            }
            bool has_gradient() const {
                return true;
            }
            pagmo::vector_double gradient(const pagmo::vector_double& vars) const {
                pagmo::vector_double grad(optimization->objectives_num + optimization->constraints_num);
                optimization->objective(&vars[0], &grad[0])[0];  // TODO
                if (optimization->constraints_num > 0) {
                    optimization->constraint(&vars[0], nullptr)[0];  // TODO
                }
                return grad;
            }
            pagmo::vector_double::size_type get_nobj() const {
                return optimization->objectives_num;
            }
            pagmo::vector_double::size_type get_nec() const {
                return 0;
            }
            pagmo::vector_double::size_type get_nic() const {
                return optimization->constraints_num;
            }
            pagmo::thread_safety get_thread_safety() const {
                return pagmo::thread_safety::none;
            }
        };
        PagmoProblem pagmo_problem;
        pagmo_problem.optimization = this;
        pagmo::problem problem{pagmo_problem};
        pagmo::population population{problem};
        population.push_back(initial_values);
        pagmo::algorithm algorithm;

        const std::string& solver_name = settings["solver"].as<std::string>();
        if (solver_name == "ipopt") {
            pagmo::ipopt solver;
            solver.set_numeric_option("tol", settings["utility_precision"].as<Value>() / 3000);
            if (settings.has("maxiter")) {
                solver.set_integer_option("max_iter", settings["maxiter"].as<size_t>());
            }
            if (settings.has("timeout")) {  // timeout given in sec
                solver.set_integer_option("max_cpu_time", settings["timeout"].as<size_t>());
            }
            solver.set_selection("best");
            algorithm = pagmo::algorithm{solver};
        } else if (solver_name == "nlopt") {
            pagmo::nlopt solver(settings["algorithm"].as<std::string>());
            solver.set_ftol_abs(settings["utility_precision"].as<Value>());
            if (settings.has("maxiter")) {
                solver.set_maxeval(settings["maxiter"].as<size_t>());
            }
            if (settings.has("timeout")) {  // timeout given in sec
                solver.set_maxtime(settings["timeout"].as<size_t>());
            }
            algorithm = pagmo::algorithm{solver};
        } else {
            throw std::runtime_error("unknown solver '" + solver_name + "'");
        }
        for (size_t i = 0; i < settings["iterations"].as<size_t>(1); ++i) {
            population = algorithm.evolve(population);
        }
        pagmo::vector_double vars = population.champion_x();
        objective(&vars[0], nullptr);
#else
        throw std::runtime_error("library '" + library + "' not supported by this binary");
#endif
    } else if (library == "borg") {
#ifdef DICEPP_WITH_BORG
        optimization = this;
        BORG_Problem opt = BORG_Problem_create(variables_num, objectives_num, constraints_num, [](double* vars, double* objs, double* consts) {
            const std::vector<Value> o = optimization->objective(vars, nullptr);
            const std::vector<Value> c = optimization->constraint(vars, nullptr);
            for (size_t i = 0; i < o.size(); ++i) {
                objs[i] = -o[i];
            }
            for (size_t i = 0; i < c.size(); ++i) {
                consts[i] = std::max(0.0, c[i]);
            }
        });

        for (Time t = 0; t < variables_num; ++t) {
            BORG_Problem_set_bounds(opt, t, 0, 1);
        }

        BORG_Problem_set_epsilon(opt, 0, settings["utility_precision"].as<Value>());
        // BORG_Random_seed(12345);
        BORG_Archive result = BORG_Algorithm_run(opt, settings["maxiter"].as<size_t>());
        // BORG_Archive_print(result, stdout);
        BORG_Archive_destroy(result);
        BORG_Problem_destroy(opt);
#else
        throw std::runtime_error("library '" + library + "' not supported by this binary");
#endif
    } else if (library == "nlopt") {
#ifdef DICEPP_WITH_NLOPT
        // "auglag nlopt::algorithm::LN_AUGLAG
        // "auglag nlopt::algorithm::LD_AUGLAG
        // "auglag_eq nlopt::algorithm::LN_AUGLAG_EQ
        // "auglag_eq nlopt::algorithm::LD_AUGLAG_EQ
        const std::string algorithm_name = settings["algorithm"].as<std::string>();
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
        } else if (algorithm_name == "lbfgs_nocedal") {
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
        } else if (algorithm_name == "gn_mlsl_lds") {
            algorithm_type = nlopt::algorithm::GN_MLSL_LDS;
        } else if (algorithm_name == "gd_mlsl_lds") {
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

        nlopt::opt opt(algorithm_type, variables_num);
        if (constraints_num > 0) {  // TODO
            opt.add_inequality_constraint(
                [](unsigned n, const double* x, double* grad, void* data) {
                    Optimization<Value, Time>* optimization = static_cast<Optimization<Value, Time>*>(data);
                    return optimization->constraint(x, grad)[0];  // TODO
                },
                this, 0.1);
        }
        opt.set_max_objective(
            [](unsigned n, const double* x, double* grad, void* data) {
                Optimization<Value, Time>* optimization = static_cast<Optimization<Value, Time>*>(data);
                return optimization->objective(x, grad)[0];  // TODO
            },
            this);

        if (settings.has("utility_precision")) {
            opt.set_ftol_abs(settings["utility_precision"].as<Value>());
        }
        if (settings.has("rel_var_precision")) {
            opt.set_xtol_rel(settings["rel_var_precision"].as<Value>());
        }
        opt.set_lower_bounds(std::vector<Value>(variables_num, 0));
        opt.set_upper_bounds(std::vector<Value>(variables_num, 1));
        if (settings.has("maxiter")) {
            opt.set_maxeval(settings["maxiter"].as<size_t>());
        }
        if (settings.has("timeout")) {  // timeout given in sec
            opt.set_maxtime(settings["timeout"].as<size_t>());
        }

        Value utility;
        nlopt::result result = opt.optimize(initial_values, utility);
        std::cout << get_optimization_results(result) << std::endl;
#else
        throw std::runtime_error("library '" + library + "' not supported by this binary");
#endif
    } else {
        throw std::runtime_error("unknown library '" + library + "'");
    }
}

template class Optimization<double, size_t>;
}
