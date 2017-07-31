#ifndef CONTROL_H
#define CONTROL_H

#include "types.h"

namespace dice {
template<typename Value, typename Time>
class Control {
  public:
    TimeSeries<Value> mu;  // Emission control rate GHGs
    TimeSeries<Value> s;   // Gross savings rate as fraction of gross world product

    Control(Time length) : mu(length, 0), s(length, 0){};

    bool observe(Observer<Value, Time>& observer) {
        OBSERVE_SERIES(mu);
        OBSERVE_SERIES(s);
        return true;
    }
};
}

#endif
