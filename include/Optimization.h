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

#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include <cstddef>
#include <vector>

#include "types.h"

namespace settings {
class SettingsNode;
}

namespace dice {
template<typename Value, typename Time>
class Optimization {
  public:
    const unsigned int variables_num;
    const unsigned int objectives_num;
    const unsigned int constraints_num;

    Optimization(unsigned int variables_num_p, unsigned int objectives_num_p, unsigned int constraints_num_p)
        : variables_num(variables_num_p), objectives_num(objectives_num_p), constraints_num(constraints_num_p) {}
    virtual ~Optimization() = default;
    void optimize(const settings::SettingsNode& settings, TimeSeries<Value>& initial_values, bool verbose);
    virtual std::vector<Value> objective(const Value* vars, Value* grad) = 0;   // to be maximized
    virtual std::vector<Value> constraint(const Value* vars, Value* grad) = 0;  // to be <= 0
};
}  // namespace dice

#endif
