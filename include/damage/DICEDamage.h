#ifndef DICEDAMAGE_H
#define DICEDAMAGE_H

#include "Climate.h"
#include "Damage.h"

namespace dice {
namespace damage {

template<typename Value, typename Time>
class DICEDamage : public Damage<Value, Time> {
  protected:
    using Damage<Value, Time>::global;
    using Damage<Value, Time>::climate;
    const settings::SettingsNode& settings;

    const Value a1{settings["a1"].as<Value>()};  // Damage intercept
    const Value a2{settings["a2"].as<Value>()};  // Damage quadratic term
    const Value a3{settings["a3"].as<Value>()};  // Damage exponent

  public:
    DICEDamage(const settings::SettingsNode& settings_p, const Global<Value, Time>& global_p, climate::Climate<Value, Time>& climate_p)
        : Damage<Value, Time>(global_p, climate_p), settings(settings_p) {
    }
    Value damfrac(Time t) override {
        return a1 * climate.T_atm(t) + a2 * std::pow(climate.T_atm(t), a3);
    }
};
}
}

#endif
