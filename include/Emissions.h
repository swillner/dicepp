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
}

#endif
