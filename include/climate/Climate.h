#ifndef CLIMATE_H
#define CLIMATE_H

#include "Emissions.h"
#include "Global.h"
#include "types.h"

namespace dice {
namespace climate {

template<typename Value, typename Time>
class Climate {
  protected:
    const Global<Value, Time>& global;
    Emissions<Value, Time>& E;

  public:
    Climate(const Global<Value, Time>& global_p, Emissions<Value, Time>& E_p) : global(global_p), E(E_p){};
    virtual ~Climate(){};
    virtual bool observe(Observer<Value, Time>& observer) {
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
