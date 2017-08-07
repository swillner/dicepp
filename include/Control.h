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
