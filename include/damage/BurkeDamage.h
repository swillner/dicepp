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

#ifndef BURKEDAMAGE_H
#define BURKEDAMAGE_H

#include "Climate.h"
#include "Damage.h"
#include "Economy.h"

namespace dice {
namespace damage {

template<typename Value, typename Time, typename Constant = Value, typename Variable = TimeSeries<Value>>
class BurkeDamage : public Damage<Value, Time, Constant, Variable> {
  protected:
    using Damage<Value, Time, Constant, Variable>::global;
    using Damage<Value, Time, Constant, Variable>::control;
    using Damage<Value, Time, Constant, Variable>::climate;
    const settings::SettingsNode& settings;

    const Constant impact1{settings["impact1"].template as<Constant>()};  // Multiplicative coefficient for linear temperature impact
    const Constant impact2{settings["impact2"].template as<Constant>()};  // Multiplicative coefficient for linear temperature impact
    const Constant T_avg{settings["T_avg"].template as<Constant>()};      // Global mean temperature average from 1980 to 2010
    const Constant T_2010{settings["T_2010"].template as<Constant>()};    // Global mean temperature in 2010 (average of the last five years)
    const Constant T_2009{settings["T_2009"].template as<Constant>()};    // Global mean temperature in 2009 (average of the last five years)

    const Constant iterstep{settings["iterstep"].template as<Constant>()};

    const Constant href = impact1 * T_avg + impact2 * T_avg * T_avg;  // Historical response function (reference)
    const Constant ddcc0 = impact1 * T_2009 + impact2 * T_2009 * T_2009 - href;
    const Constant produc0 = 1 + ddcc0;

    bool has_final_f = false;
    TimeSeries<Constant> f;
    LinearInterpolator<Value> final_f;

    Constant T_abs(Time t) {
        return climate.T_atm(t).value() - climate.T_atm(0).value() + T_2010;
    }
    Constant h(Time t) {
        return impact1 * T_abs(t) + impact2 * T_abs(t) * T_abs(t);
    }
    Constant phi(Time t) {
        return h(t) - href;
    }
    Constant K_gross(Time t, Economy<Value, Time, Constant, Variable>& economy, const TimeSeries<Constant>& s_fix) {
        if (t == 0) {
            return economy.K(0).value();
        } else {
            return s_fix[t - 1] * economy.Y_gross(t - 1).value() + (1 - global.dK) * economy.K(t - 1).value();
        }
    }
    Constant Y_gross(Time t, Economy<Value, Time, Constant, Variable>& economy, const TimeSeries<Constant>& s_fix) {
        return economy.A(t) * std::pow(economy.L(t) / 1000, 1 - global.gamma) * std::pow(K_gross(t, economy, s_fix), global.gamma);
    }
    Constant eta(Time t, Economy<Value, Time, Constant, Variable>& economy, const TimeSeries<Constant>& s_fix) {
        if (t < global.timestep_num - 1) {
            return Y_gross(t + 1, economy, s_fix) * economy.L(t) / economy.Y_gross(t).value() / economy.L(t + 1) - 1;
        } else {
            return eta(t - 1, economy, s_fix);
        }
    }
    Constant f_fit(Time t, Economy<Value, Time, Constant, Variable>& economy, const TimeSeries<Constant>& s_fix) {
        if (t == 0) {
            return produc0;
        } else {
            return economy.L(t) * economy.Y_net(t - 1).value() * (1 + eta(t - 1, economy, s_fix) + phi(t - 1)) / economy.L(t - 1) / economy.Y_gross(t).value();
        }
    }

  public:
    BurkeDamage(const settings::SettingsNode& settings_p,
                const Global<Constant, Time>& global_p,
                const Control<Value, Time, Constant, Variable>& control_p,
                climate::Climate<Value, Time, Constant, Variable>& climate_p)
        : Damage<Value, Time, Constant, Variable>(global_p, control_p, climate_p), settings(settings_p), f(global_p.timestep_num, 1){};

    Constant phi_diff(Economy<Value, Time, Constant, Variable>& economy, const TimeSeries<Constant>& s_fix) {
        Constant sum = 0;
        for (Time t = 0; t < global.timestep_num - 1; ++t) {
            sum += std::abs(phi(t)
                            - (economy.L(t) * (economy.Y_net(t + 1) / economy.L(t + 1) - economy.Y_net(t) * (1 + eta(t, economy, s_fix)) / economy.L(t))
                               / (economy.Y_net(t)))
                                  .value());
        }
        return sum / global.timestep_num;
    }
    void recalc_f(Economy<Value, Time, Constant, Variable>& economy, const TimeSeries<Constant>& s_fix) {
        f[0] = produc0;
        for (Time t = 1; t < global.timestep_num; ++t) {
            f[t] = std::min(f[t] + (f_fit(t, economy, s_fix) - f[t]) / iterstep, 1.0);
        }
    }
    void calc_final_f() {
        final_f.data.reserve(global.timestep_num);
        Constant T_last = 0;
        for (Time t = 0; t < global.timestep_num; ++t) {
            const Constant T_t = climate.T_atm(t).value();
            if (T_t <= T_last) {
                break;
            } else {
                final_f.data.emplace_back(std::make_pair<Value, Value>({control.variables_num, T_t}, {control.variables_num, f[t]}));
                T_last = T_t;
            }
        }
        has_final_f = true;
    }
    Value damfrac(Time t) override {
        if (has_final_f) {
            return 1 - final_f(climate.T_atm(t));
        } else {
            return {control.variables_num, 1 - f[t]};
        }
    }
};
}
}

#endif
