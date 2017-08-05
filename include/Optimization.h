#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include <vector>
#include "types.h"

namespace settings {
class SettingsNode;
}

namespace dice {
template<typename Value, typename Time>
class Optimization {
  public:
    const size_t variables_num;
    const size_t objectives_num;
    const size_t constraints_num;

    Optimization(size_t variables_num_p, size_t objectives_num_p, size_t constraints_num_p)
        : variables_num(variables_num_p), objectives_num(objectives_num_p), constraints_num(constraints_num_p){};
    virtual ~Optimization(){};
    void optimize(const settings::SettingsNode& settings, TimeSeries<Value>& initial_values);
    virtual std::vector<Value> objective(const Value* vars, Value* grad) = 0;  // to be maximized
    virtual std::vector<Value> constraint(const Value* vars, Value* grad) {    // to be <= 0
        return {0};                                                            // TODO
    }
};
}

#endif
