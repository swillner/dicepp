#ifndef EMISSIONS_H
#define EMISSIONS_H

#include <vector>
#include "Economy.h"
#include "Global.h"
#include "types.h"

namespace dice {

template<typename Value, typename Time>
class Emissions {
  protected:
    const Global<Value, Time>& global;
    std::vector<Economy<Value, Time>>& economies;
    StepwiseBackwardLookingTimeSeries<Value, Time> E_series{global.timestep_num, 0};

  public:
    Emissions(const Global<Value, Time>& global_p, std::vector<Economy<Value, Time>>& economies_p)
        : global(global_p), economies(economies_p){};

    // Total CO2 emissions (GtCO2 per year)
    Value operator()(Time t) {
        return E_series.get(t, [this](Time t, Value E_last) {

            Value E = 0;
            for (auto&& economy : economies) {
                E += economy.E(t);
            }
            return E;

        });
    }
    void reset() {
        E_series.reset();
    }
};
}

#endif
