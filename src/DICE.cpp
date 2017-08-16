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

#include "DICE.h"
#include <fstream>  // IWYU pragma: keep
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>  // IWYU pragma: keep
#include <stdexcept>
#include <string>
#include <tuple>
#include "BurkeDamage.h"
#include "DICEClimate.h"
#include "DICEDamage.h"
#include "Optimization.h"
#include "autodiff.h"
#include "csv-parser.h"
#include "settingsnode.h"

#ifdef DICEPP_WITH_NETCDF
#include <ncDim.h>
#include <ncFile.h>
#include <ncGroup.h>
#include <ncType.h>
#include <ncVar.h>
#endif

namespace dice {

template<typename Value, typename Time>
DICE<Value, Time>::DICE(const settings::SettingsNode& settings_p)
    : settings(settings_p),
      global(settings_p["parameters"]),
      control(settings_p["parameters"]["timestep_num"].as<Time>()),
      emissions(global, control, economies) {
}

template<typename Value, typename Time>
void DICE<Value, Time>::initialize() {
    std::cout << std::setprecision(13);

    // Initialize climate module
    {
        const settings::SettingsNode& climate_node = settings["climate"];
        const std::string& type = climate_node["type"].as<std::string>();
        if (type == "dice") {
            climate.reset(new climate::DICEClimate<autodiff::Value<Value>, Time, Value, autodiff::Variable<Value>>(climate_node["parameters"], global, control,
                                                                                                                   emissions));
        } else {
            throw std::runtime_error("unknown climate module type '" + type + "'");
        }
        climate->initialize();
    }

    // Initialize damage module
    {
        const settings::SettingsNode& damage_node = settings["damage"];
        const std::string& type = damage_node["type"].as<std::string>();
        if (type == "burke") {
            damage.reset(
                new damage::BurkeDamage<autodiff::Value<Value>, Time, Value, autodiff::Variable<Value>>(damage_node["parameters"], global, control, *climate));
        } else if (type == "dice") {
            damage.reset(
                new damage::DICEDamage<autodiff::Value<Value>, Time, Value, autodiff::Variable<Value>>(damage_node["parameters"], global, control, *climate));
        } else {
            throw std::runtime_error("unknown damage module type '" + type + "'");
        }
        damage->initialize();
    }

    // Initialize regions
    {
        for (const auto&& region_node : settings["regions"].as_sequence()) {
            economies.emplace_back(
                Economy<autodiff::Value<Value>, Time, Value, autodiff::Variable<Value>>(region_node["economy"], global, control, *climate, *damage));
        }
    }

    emissions.initialize();

    // Initialize control variables
    if (settings.has("control")) {
        class ControlInputObserver : public Observer<autodiff::Value<Value>, Time, Value> {
          protected:
            const settings::SettingsNode& input_node;

          public:
            explicit ControlInputObserver(const settings::SettingsNode& input_node_p) : input_node(input_node_p){};
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

template<typename Value, typename Time>
void DICE<Value, Time>::single_optimization(Optimization<Value, Time>& optimization,
                                            const settings::SettingsNode& optimization_node,
                                            TimeSeries<Value>& initial_values,
                                            bool verbose) {
    for (const auto& iteration_node : optimization_node["iterations"].as_sequence()) {
        for (size_t i = 0; i < iteration_node["repeat"].as<size_t>(1); ++i) {
            std::copy(std::begin(control.s.value()), std::begin(control.s.value()) + initial_values.size(), std::begin(initial_values));
            optimization.optimize(iteration_node, initial_values, verbose);
            if (verbose) {
                reset();
                const autodiff::Value<Value> utility = calc_single_utility();
                Value sum = 0;
                for (size_t i = 0; i < initial_values.size(); ++i) {
                    sum += utility.derivative()[i] * utility.derivative()[i];
                }
                std::cout << "Gradient length = " << std::sqrt(sum) << std::endl;
                std::cout << "Finished with utility = " << utility.value() << std::endl;
            }
        }
    }
}

template<typename Value, typename Time>
void DICE<Value, Time>::run() {
    if (economies.empty()) {
        throw std::runtime_error("no economies given");
    }
    if (economies.size() == 1) {
        const settings::SettingsNode& optimization_node = settings["optimization"];
        if (optimization_node.has("iterations")) {
            class DICEOptimization : public Optimization<Value, Time> {
              protected:
                DICE& dice;

              public:
                using Optimization<Value, Time>::variables_num;
                using Optimization<Value, Time>::objectives_num;
                using Optimization<Value, Time>::constraints_num;
                DICEOptimization(size_t variables_num_p, size_t objectives_num_p, size_t constraints_num_p, DICE& dice_p)
                    : Optimization<Value, Time>(variables_num_p, objectives_num_p, constraints_num_p), dice(dice_p){};

                std::vector<Value> objective(const Value* vars, Value* grad) override {
#ifdef DEBUG
                    try {
#endif
                        std::copy(vars, vars + variables_num, std::begin(dice.control.s.value()));
                        dice.reset();
                        const autodiff::Value<Value> utility = dice.calc_single_utility();
                        if (grad) {
                            std::copy(&utility.derivative()[0], &utility.derivative()[variables_num], grad);
                        }
                        return {utility.value()};
#ifdef DEBUG
                    } catch (std::exception& e) {
                        std::cerr << "Exception '" << e.what() << "' in optimization" << std::endl;
                        throw;
                    }
#endif
                }

                std::vector<Value> constraint(const Value* vars, Value* grad) override {
#ifdef DEBUG
                    try {
#endif
                        std::copy(vars, vars + variables_num, std::begin(dice.control.s.value()));
                        dice.reset();
                        const autodiff::Value<Value> c = dice.economies[0].cca(dice.global.timestep_num - 1) - dice.global.fosslim;
                        if (grad) {
                            std::copy(&c.derivative()[0], &c.derivative()[variables_num], grad);
                        }
                        return {c.value()};
#ifdef DEBUG
                    } catch (std::exception& e) {
                        std::cerr << "Exception '" << e.what() << "' in optimization" << std::endl;
                        throw;
                    }
#endif
                }
            };

            const size_t optimization_variables_num = global.timestep_num - optimization_node["s_fix_steps"].as<Time>(0);
            const size_t constraints_num = optimization_node["limit_cca"].as<bool>() ? 1 : 0;
            const bool verbose = optimization_node["verbose"].as<bool>();
            DICEOptimization optimization{optimization_variables_num, 1, constraints_num, *this};
            std::fill(std::begin(control.s.value()), std::end(control.s.value()), global.optlrsav);
            TimeSeries<Value> initial_values(optimization_variables_num, 0);

            const settings::SettingsNode& damage_node = settings["damage"];
            const std::string& type = damage_node["type"].as<std::string>();
            if (type == "burke") {
                const auto& iterations_node = optimization_node["iterations"];
                auto& burke_damage = *static_cast<damage::BurkeDamage<autodiff::Value<Value>, Time, Value, autodiff::Variable<Value>>*>(damage.get());
                auto& economy = economies[0];
                const Value iterstop = iterations_node["iterstop"].as<Value>();
                Value phi_diff = 1;
                single_optimization(optimization, iterations_node["before"], initial_values, verbose);
                TimeSeries<Value> s_fix(control.s.value());
                size_t i = 0;
                while (phi_diff > iterstop) {
                    std::cout << "iterstop = " << iterstop << " < phi_diff = " << phi_diff << std::endl;
                    burke_damage.recalc_f(economy, s_fix);
                    if (i % 5 == 0) {
                        single_optimization(optimization, iterations_node["after"], initial_values, verbose);
                    } else {
                        single_optimization(optimization, iterations_node["during"], initial_values, verbose);
                    }
                    ++i;
                    phi_diff = burke_damage.phi_diff(economy, s_fix);
                    if (phi_diff <= iterstop) {
                        phi_diff = 1;
                        while (phi_diff > iterstop) {
                            std::cout << "iterstop = " << iterstop << " < phi_diff = " << phi_diff << std::endl;
                            burke_damage.recalc_f(economy, s_fix);
                            single_optimization(optimization, iterations_node["after"], initial_values, verbose);
                            phi_diff = burke_damage.phi_diff(economy, s_fix);
                        }
                    }
                }
                burke_damage.calc_final_f();
                std::fill(std::begin(control.s.value()), std::end(control.s.value()), global.optlrsav);
                single_optimization(optimization, iterations_node["after"], initial_values, verbose);
            } else {
                single_optimization(optimization, optimization_node, initial_values, verbose);
            }
        } else {
            const size_t optimization_variables_num = global.timestep_num - 10;
            const autodiff::Value<Value> utility = calc_single_utility();
            Value sum = 0;
            for (size_t i = 0; i < optimization_variables_num; ++i) {
                sum += utility.derivative()[i] * utility.derivative()[i];
            }
            std::cout << "Gradient length = " << std::sqrt(sum) << std::endl;
            std::cout << "Finished with utility = " << utility.value() << std::endl;
        }
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
autodiff::Value<Value> DICE<Value, Time>::calc_single_utility() {
    autodiff::Value<Value> utility{control.variables_num, 0};
    for (Time t = 0; t < global.timestep_num; ++t) {
        utility += economies[0].utility(t);
    }
    return global.scale1 * utility + global.scale2;
}

template<typename Value, typename Time>
void DICE<Value, Time>::output() {
    if (settings.has("output")) {
        const settings::SettingsNode& output_node = settings["output"];
        const std::string& type = output_node["type"].as<std::string>();
        if (type == "netcdf") {
#ifdef DICEPP_WITH_NETCDF
            write_netcdf_output(output_node);
#else
            throw std::runtime_error("output type '" + type + "' not supported by this binary");
#endif
        } else if (type == "csv") {
            write_csv_output(output_node);
        } else {
            throw std::runtime_error("unknown output type '" + type + "'");
        }
    }
}

#ifdef DICEPP_WITH_NETCDF
template<typename Value, typename Time>
void DICE<Value, Time>::write_netcdf_output(const settings::SettingsNode& output_node) {
    if (economies.size() == 1) {
        netCDF::NcFile file(output_node["filename"].as<std::string>(), netCDF::NcFile::replace, netCDF::NcFile::nc4);

        netCDF::NcDim time_dim = file.addDim("time", global.timestep_num);
        netCDF::NcVar time_var = file.addVar("time", netCDF::NcType::nc_UINT, {time_dim});
        for (Time t = 0; t < global.timestep_num; ++t) {
            const Time year = global.start_year + t * global.timestep_length;
            time_var.putVar({t}, static_cast<const unsigned int>(year));
        }
        class NetCDFOutputObserver : public Observer<autodiff::Value<Value>, Time, Value> {
          protected:
            const netCDF::NcGroup& group;
            const netCDF::NcDim& time_dim;
            const settings::SettingsNode& output_node;

          public:
            NetCDFOutputObserver(const netCDF::NcGroup& group_p, const netCDF::NcDim& time_dim_p, const settings::SettingsNode& output_node_p)
                : group(group_p), time_dim(time_dim_p), output_node(output_node_p){};
            std::tuple<bool, bool, Time> want(const std::string& name) override {
                (void)name;
                return {true, true, 0};
            }
            bool observe(const std::string& name, TimeSeries<Value>& v) override {
                netCDF::NcVar var = group.addVar(name, netCDF::NcType::nc_FLOAT, {time_dim});
                var.setCompression(false, true, 7);
                var.putVar(&v[0]);
                return true;
            }
        };
        NetCDFOutputObserver observer(file, time_dim, output_node);
        economies[0].observe(observer);
        climate->observe(observer);
        damage->observe(observer);
        control.observe(observer);
        emissions.observe(observer);
        file.putAtt("utility", netCDF::NcType::nc_FLOAT, calc_single_utility().value());
    } else {
        throw std::runtime_error("multiple regions not supported yet");
    }
}
#endif

template<typename Value, typename Time>
void DICE<Value, Time>::write_csv_output(const settings::SettingsNode& output_node) {
    if (economies.size() == 1) {
        const std::string& filename = output_node["filename"].as<std::string>();
        std::ofstream file(filename);
        if (!file) {
            throw std::runtime_error("could not write to '" + filename + "'");
        }

        class CSVOutputObserver : public Observer<autodiff::Value<Value>, Time, Value> {
          protected:
            std::ofstream& file;

          public:
            Time t;
            std::string var;

            explicit CSVOutputObserver(std::ofstream& file_p) : file(file_p){};
            std::tuple<bool, bool, Time> want(const std::string& name) override {
                return std::tuple<bool, bool, Time>(name == var, false, t);
            }
            bool observe(const std::string& name, const autodiff::Value<Value>& v) override {
                (void)name;
                file << v.value();
                return false;
            }
            bool observe(const std::string& name, const Value& v) override {
                (void)name;
                file << v;
                return false;
            }
            bool observe(const std::string& name, TimeSeries<Value>& v) override {
                if (name == var) {
                    file << v[t];
                    return false;
                }
                return true;
            }
        };
        const auto& variables = output_node["columns"].as_sequence();
        CSVOutputObserver observer(file);
        for (auto&& var = std::begin(variables); var != std::end(variables); ++var) {
            if (var != std::begin(variables)) {
                file << ",";
            }
            file << "\"" << (*var).as<std::string>() << "\"";
        }
        file << "\n";
        const auto utility = calc_single_utility();
        const auto& dev = utility.derivative();
        for (Time t = 0; t < global.timestep_num; ++t) {
            observer.t = t;
            for (auto&& var = std::begin(variables); var != std::end(variables); ++var) {
                if (var != std::begin(variables)) {
                    file << ",";
                }
                const std::string& name = (*var).as<std::string>();
                if (name == "t") {
                    file << t;
                } else if (name == "year") {
                    file << (global.start_year + t * global.timestep_length);
                } else if (name == "utility") {
                    file << utility.value();
                } else if (name == "gradient") {
                    file << dev[t];
                } else {
                    observer.var = name;
                    if (economies[0].observe(observer) && climate->observe(observer) && damage->observe(observer) && control.observe(observer)
                        && emissions.observe(observer)) {
                        throw std::runtime_error("variable '" + name + "' not found");
                    }
                }
            }
            file << "\n";
        }
    } else {
        throw std::runtime_error("multiple regions not supported yet");
    }
}

template class DICE<double, size_t>;
}  // namespace dice
