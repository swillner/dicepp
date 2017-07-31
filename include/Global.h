#ifndef GLOBAL_H
#define GLOBAL_H

#include "settingsnode.h"

namespace dice {
template<typename Value, typename Time>
class Global {
  protected:
    const settings::SettingsNode& settings;

  public:
    const Value dK{settings["dK"].as<Value>()};              // Depreciation rate on capital (per year)
    const Value elasmu{settings["elasmu"].as<Value>()};      // Elasticity of marginal utility of consumption
    const Value expcost2{settings["expcost2"].as<Value>()};  // Exponent of control cost function
    const Value fosslim{settings["fosslim"].as<Value>()};    // Maximum cumulative extraction fossil fuels (GtC)
    const Value gamma{settings["gamma"].as<Value>()};        // Capital elasticity in production function
    const Value prstp{settings["prstp"].as<Value>()};        // Initial rate of social time preference per year
    const Value scale1{settings["scale1"].as<Value>()};      // Multiplicative scaling coefficient
    const Value scale2{settings["scale2"].as<Value>()};      // Additive scaling coefficient

    const Time timestep_length{settings["timestep_length"].as<Time>()};
    const Time start_year{settings["start_year"].as<Time>()};
    const Time timestep_num{settings["timestep_num"].as<Time>()};

    const Value optlrsav = (dK + 0.004) / (dK + 0.004 * elasmu + prstp) * gamma;  // Optimal long-run savings rate used for transversality

    Global(const settings::SettingsNode& settings_p) : settings(settings_p){};
};
}

#endif
