#ifndef DAMAGE_H
#define DAMAGE_H

#include "types.h"

namespace dice {

template<typename Value, typename Time>
class Global;

namespace climate {

template<typename Value, typename Time, typename Constant, typename Variable>
class Climate;
}

namespace damage {

template<typename Value, typename Time, typename Constant = Value, typename Variable = TimeSeries<Value>>
class Damage {
  protected:
    const Global<Constant, Time>& global;
    climate::Climate<Value, Time, Constant, Variable>& climate;

  public:
    Damage(const Global<Constant, Time>& global_p, climate::Climate<Value, Time, Constant, Variable>& climate_p) : global(global_p), climate(climate_p){};
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
