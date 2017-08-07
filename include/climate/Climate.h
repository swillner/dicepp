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

#ifndef CLIMATE_H
#define CLIMATE_H

#include "Control.h"
#include "Emissions.h"
#include "Global.h"
#include "types.h"

namespace dice {
namespace climate {

template<typename Value, typename Time, typename Constant = Value, typename Variable = TimeSeries<Value>>
class Climate {
  protected:
    const Global<Constant, Time>& global;
    const Control<Value, Time, Constant, Variable>& control;
    Emissions<Value, Time, Constant, Variable>& E;

  public:
    Climate(const Global<Constant, Time>& global_p, const Control<Value, Time, Constant, Variable>& control_p, Emissions<Value, Time, Constant, Variable>& E_p)
        : global(global_p), control(control_p), E(E_p){};
    virtual ~Climate(){};
    virtual bool observe(Observer<Value, Time, Constant>& observer) {
        OBSERVE_VAR(T_atm);
        return true;
    }
    virtual void initialize(){};
    virtual void reset(){};
    virtual Value T_atm(Time t) = 0;
};
}
}

#endif
