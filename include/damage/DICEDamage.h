#ifndef DICEDAMAGE_H
#define DICEDAMAGE_H

#include "Climate.h"
#include "Damage.h"

namespace dice {
namespace damage {

template<typename Value, typename Time, typename Constant = Value, typename Variable = TimeSeries<Value>>
class DICEDamage : public Damage<Value, Time, Constant, Variable> {
  protected:
    using Damage<Value, Time, Constant, Variable>::global;
    using Damage<Value, Time, Constant, Variable>::climate;
    const settings::SettingsNode& settings;

    const Constant a1{settings["a1"].as<Constant>()};  // Damage intercept
    const Constant a2{settings["a2"].as<Constant>()};  // Damage quadratic term
    const Constant a3{settings["a3"].as<Constant>()};  // Damage exponent

  public:
    DICEDamage(const settings::SettingsNode& settings_p, const Global<Constant, Time>& global_p, climate::Climate<Value, Time, Constant, Variable>& climate_p)
        : Damage<Value, Time, Constant, Variable>(global_p, climate_p), settings(settings_p) {
    }
    Value damfrac(Time t) override {
        return a1 * climate.T_atm(t) + a2 * std::pow(climate.T_atm(t), a3);
    }
};
}
}

#endif
