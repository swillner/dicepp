#ifndef GLOBAL_H
#define GLOBAL_H

#include "settingsnode.h"

namespace dice {
template<typename Constant, typename Time>
class Global {
  protected:
    const settings::SettingsNode& settings;

  public:
    const Constant dK{settings["dK"].as<Constant>()};              // Depreciation rate on capital (per year)
    const Constant elasmu{settings["elasmu"].as<Constant>()};      // Elasticity of marginal utility of consumption
    const Constant expcost2{settings["expcost2"].as<Constant>()};  // Exponent of control cost function
    const Constant fosslim{settings["fosslim"].as<Constant>()};    // Maximum cumulative extraction fossil fuels (GtC)
    const Constant gamma{settings["gamma"].as<Constant>()};        // Capital elasticity in production function
    const Constant prstp{settings["prstp"].as<Constant>()};        // Initial rate of social time preference per year
    const Constant scale1{settings["scale1"].as<Constant>()};      // Multiplicative scaling coefficient
    const Constant scale2{settings["scale2"].as<Constant>()};      // Additive scaling coefficient

    const Time timestep_length{settings["timestep_length"].as<Time>()};
    const Time start_year{settings["start_year"].as<Time>()};
    const Time timestep_num{settings["timestep_num"].as<Time>()};

    const Constant optlrsav = (dK + 0.004) / (dK + 0.004 * elasmu + prstp) * gamma;  // Optimal long-run savings rate used for transversality

    Global(const settings::SettingsNode& settings_p) : settings(settings_p){};
};
}

#endif
