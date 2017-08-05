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
