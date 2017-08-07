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

#ifndef CONTROL_H
#define CONTROL_H

#include "types.h"

namespace dice {
template<typename Value, typename Time, typename Constant = Value, typename Variable = TimeSeries<Value>>
class Control {
  public:
    const size_t variables_num;
    Variable mu{variables_num, variables_num, variables_num, 0};  // Emission control rate GHGs
    Variable s{0, variables_num, variables_num, 0};               // Gross savings rate as fraction of gross world product

    Control(Time length) : variables_num(length){};  // TODO mu

    bool observe(Observer<Value, Time, Constant>& observer) {
        OBSERVE_VARIABLE(mu);
        OBSERVE_VARIABLE(s);
        return true;
    }
};
}

#endif
