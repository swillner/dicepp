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

#ifndef EMISSIONS_H
#define EMISSIONS_H

#include <vector>
#include "Control.h"
#include "Economy.h"
#include "Global.h"
#include "types.h"

namespace dice {

template<typename Value, typename Time, typename Constant = Value, typename Variable = TimeSeries<Value>>
class Emissions {
  protected:
    const Global<Constant, Time>& global;
    const Control<Value, Time, Constant, Variable>& control;
    std::vector<Economy<Value, Time, Constant, Variable>>& economies;
    StepwiseBackwardLookingTimeSeries<Value, Time> E_series{global.timestep_num, {control.variables_num, 0}};

  public:
    Emissions(const Global<Constant, Time>& global_p,
              const Control<Value, Time, Constant, Variable>& control_p,
              std::vector<Economy<Value, Time, Constant, Variable>>& economies_p)
        : global(global_p), control(control_p), economies(economies_p){};

    // Total CO2 emissions (GtCO2 per year)
    Value operator()(Time t) {
        return E_series.get(t, [this](Time t, Value E_last) {
            (void)E_last;

            Value E{control.variables_num, 0};
            for (auto&& economy : economies) {
                E += economy.E(t);
            }
            return E;

        });
    }
    void initialize() {
        Value E{control.variables_num, 0};
        for (auto&& economy : economies) {
            E += economy.E(0);
        }
        E_series.set_first_value(E);
    }
    void reset() {
        E_series.reset();
    }
    bool observe(Observer<Value, Time, Constant>& observer) {
        return observer.observe("E_total", *this, global.timestep_num);
    }
};
}  // namespace dice

#endif
