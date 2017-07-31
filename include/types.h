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
#define OBSERVE_SERIES(v)           \
    if (!observer.observe(#v, v)) { \
        return false;               \
    }

template<typename Value, typename Time>
class Observer {
  public:
    template<typename Var>
    bool observe(const std::string& name, Var v, Time length) {
        TimeSeries<Value> series(length);
        for (Time t = 0; t < length; ++t) {
            series[t] = v(t);
        }
        return observe(name, series);
    }
    virtual bool observe(const std::string& name, TimeSeries<Value>& v) = 0;
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

template<typename Type>
class Parameter {
  protected:
    Type value;

  public:
    Parameter(const std::string& name_p, const settings::SettingsNode& settings) : value(settings[name_p].as<Type>()){};
    Parameter(const std::string& name_p, const Type& value_p) : value(value_p){};
    operator Type() const {
        return value;
    }
};

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
    Bounded(Value val_p, Value lower_p, Value upper_p) : val(val_p), lower(lower_p), upper(upper_p) {
        check();
    }
    inline Bounded& operator=(Value val_p) {
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
    LowerBounded(Value val_p, Value lower_p) : val(val_p), lower(lower_p) {
        check();
    }
    inline LowerBounded& operator=(Value val_p) {
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
    BackwardLookingTimeSeries(Time size, Value initial_value_p) : series(size, initial_value_p), initial_value(initial_value_p) {
        series[0] = initial_value;
    }

    template<typename Function>
    inline Value get(Time t, Function func) {
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
        std::fill(series.begin(), series.end(), initial_value);
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
    StepwiseBackwardLookingTimeSeries(Time size, Value initial_value_p) : series(size, initial_value_p), initial_value(initial_value_p) {
        series[0] = initial_value;
    }

    template<typename Function>
    inline Value get(Time t, Function func) {
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
        std::fill(series.begin(), series.end(), initial_value);
    }
};
}

#endif
