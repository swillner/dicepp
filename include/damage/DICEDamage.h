/*
  Copyright (C) 2017-2020 Sven Willner <sven.willner@pik-potsdam.de>

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

    const Constant a1{settings["a1"].template as<Constant>()};  // Damage intercept
    const Constant a2{settings["a2"].template as<Constant>()};  // Damage quadratic term
    const Constant a3{settings["a3"].template as<Constant>()};  // Damage exponent

  public:
    DICEDamage(const settings::SettingsNode& settings_p,
               const Global<Constant, Time>& global_p,
               const Control<Value, Time, Constant, Variable>& control_p,
               climate::Climate<Value, Time, Constant, Variable>& climate_p)
        : Damage<Value, Time, Constant, Variable>(global_p, control_p, climate_p), settings(settings_p) {}
    Value damfrac(Time t) override { return a1 * climate.T_atm(t) + a2 * std::pow(climate.T_atm(t), a3); }
};
}  // namespace damage
}  // namespace dice

#endif
