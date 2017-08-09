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

#ifndef TYPES_H
#define TYPES_H

#include <functional>
#include <string>
#include <vector>
#include "settingsnode.h"

namespace dice {

template<typename Value>
using TimeSeries = std::vector<Value>;

#define OBSERVE_VAR(v)                                                                 \
    if (!observer.observe(#v, [this](Time t) { return v(t); }, global.timestep_num)) { \
        return false;                                                                  \
    }
#define OBSERVE_VARIABLE(v)                 \
    if (!observer.observe(#v, v.value())) { \
        return false;                       \
    }

template<typename Value, typename Time, typename Constant = Value>
class Observer {
  public:
    template<typename Var>
    bool observe(const std::string& name, Var v, Time length) {
        bool want_it;
        bool want_all;
        Time want_t;
        std::tie(want_it, want_all, want_t) = want(name);
        if (want_it) {
            if (want_all) {
                TimeSeries<Constant> series;
                series.reserve(length);
                for (Time t = 0; t < length; ++t) {
                    series.emplace_back(static_cast<Constant>(v(t)));
                }
                return observe(name, series);
            } else {
                return observe(name, v(want_t));
            }
        } else {
            return true;
        }
    }
    virtual std::tuple<bool, bool, Time> want(const std::string& name) {
        return std::tuple<bool, bool, Time>(false, false, 0);
    }
    virtual bool observe(const std::string& name, TimeSeries<Constant>& v) = 0;
    virtual bool observe(const std::string& name, const Value& v) {
        return true;
    }
    virtual bool observe(const std::string& name, const Constant& v) {
        return true;
    }
};

template<typename... Types>
inline void debug_(Types... args);
template<typename Type1, typename... Types>
inline void debug_(Type1 arg1, Types... rest) {
#ifdef DEBUG
    std::cout << arg1;
    debug_(rest...);
#endif
}
template<>
inline void debug_() {
#ifdef DEBUG
    std::cout << std::flush;
#endif
}

template<typename... Types>
inline void debug(Types... args);
template<typename Type1, typename... Types>
inline void debug(Type1 arg1, Types... rest) {
#ifdef DEBUG
    std::cout << arg1;
    debug(rest...);
#endif
}
template<>
inline void debug() {
#ifdef DEBUG
    std::cout << std::endl;
#endif
}

template<typename Value>
class Bounded {
  private:
    Value val;
    Value lower;
    Value upper;
    inline void check() {
        if (val < lower) {
            val = lower;
        } else if (val > upper) {
            val = upper;
        }
    }

  public:
    Bounded(const Value& val_p, const Value& lower_p, const Value& upper_p) : val(val_p), lower(lower_p), upper(upper_p) {
        check();
    }
    inline Bounded& operator=(const Value& val_p) {
        val = val_p;
        check();
        return *this;
    }
    inline operator Value() {
        return val;
    }
    inline operator Value() const {
        return val;
    }
};

template<typename Value>
class LowerBounded {
  private:
    Value val;
    Value lower;
    inline void check() {
        if (val < lower) {
            val = lower;
        }
    }

  public:
    LowerBounded(const Value& val_p, const Value& lower_p) : val(val_p), lower(lower_p) {
        check();
    }
    inline LowerBounded& operator=(const Value& val_p) {
        val = val_p;
        check();
        return *this;
    }
    inline operator Value() {
        return val;
    }
    inline operator Value() const {
        return val;
    }
};

template<typename Value, typename Time>
class BackwardLookingTimeSeries {
  protected:
    TimeSeries<Value> series;
    Time largest_valid_t = 0;
#ifdef DEBUG
    Time calculating_t = 0;
#endif

  public:
    const Value initial_value;
    BackwardLookingTimeSeries(Time size, const Value& initial_value_p) : series(size, initial_value_p), initial_value(initial_value_p){};
    void set_first_value(const Value& first_value) {
        series[0] = first_value;
    }
    const Value& get_first_value() const {
        return series[0];
    }

    template<typename Function>
    inline const Value& get(Time t, Function func) {
        if (t > largest_valid_t) {
#ifdef DEBUG
            if (calculating_t > 0 && calculating_t <= t) {
                throw std::runtime_error("equation loop");
            }
            calculating_t = t;
#endif
            largest_valid_t = func(largest_valid_t, t, series);
#ifdef DEBUG
            calculating_t = 0;
#endif
        }
        return series[t];
    }

    inline void invalidate_after(Time t) {
        largest_valid_t = t;
    }
    inline void invalidate() {
        invalidate_after(0);
    }
    inline void reset() {
        invalidate();
        const Value first_value = series[0];
        std::fill(series.begin(), series.end(), initial_value);
        series[0] = first_value;
    }
};

template<typename Value, typename Time>
class StepwiseBackwardLookingTimeSeries {
  protected:
    TimeSeries<Value> series;
    Time largest_valid_t = 0;
#ifdef DEBUG
    Time calculating_t = 0;
#endif

  public:
    const Value initial_value;
    StepwiseBackwardLookingTimeSeries(Time size, const Value& initial_value_p) : series(size, initial_value_p), initial_value(initial_value_p){};
    void set_first_value(const Value& first_value) {
        series[0] = first_value;
    }
    const Value& get_first_value() const {
        return series[0];
    }

    template<typename Function>
    inline const Value& get(Time t, Function func) {
#ifdef DEBUG
        if (calculating_t > 0 && calculating_t <= t && t > largest_valid_t) {
            throw std::runtime_error("equation loop");
        }
        calculating_t = t;
#endif
        for (; t > largest_valid_t; ++largest_valid_t) {
            series[largest_valid_t + 1] = func(largest_valid_t + 1, series[largest_valid_t]);
        }
#ifdef DEBUG
        calculating_t = 0;
#endif
        return series[t];
    }

    inline void invalidate_after(Time t) {
        largest_valid_t = t;
    }
    inline void invalidate() {
        invalidate_after(0);
    }
    inline void reset() {
        invalidate();
        const Value first_value = series[0];
        std::fill(series.begin(), series.end(), initial_value);
        series[0] = first_value;
    }
};
}

#endif
