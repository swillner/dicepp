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

#ifndef DAMAGE_H
#define DAMAGE_H

#include "types.h"

namespace dice {

template<typename Value, typename Time>
class Global;
template<typename Value, typename Time, typename Constant, typename Variable>
class Control;

namespace climate {
template<typename Value, typename Time, typename Constant, typename Variable>
class Climate;
}

namespace damage {

template<typename Value, typename Time, typename Constant = Value, typename Variable = TimeSeries<Value>>
class Damage {
  protected:
    const Global<Constant, Time>& global;
    const Control<Value, Time, Constant, Variable>& control;
    climate::Climate<Value, Time, Constant, Variable>& climate;

  public:
    Damage(const Global<Constant, Time>& global_p,
           const Control<Value, Time, Constant, Variable>& control_p,
           climate::Climate<Value, Time, Constant, Variable>& climate_p)
        : global(global_p), climate(climate_p), control(control_p){};
    virtual ~Damage(){};
    virtual bool observe(Observer<Value, Time, Constant>& observer) {
        OBSERVE_VAR(damfrac);
        return true;
    }
    virtual void initialize(){};
    virtual void reset(){};
    virtual Value damfrac(Time t) = 0;  // Damages as fraction of gross output
};
}
}

#endif
