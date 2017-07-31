#ifndef DICE_H
#define DICE_H

#include <memory>
#include <vector>
#include "Climate.h"
#include "Control.h"
#include "Damage.h"
#include "Economy.h"
#include "Global.h"

namespace settings {
class SettingsNode;
}

namespace netCDF {
class NcGroup;
class NcDim;
}

namespace dice {

template<typename Value, typename Time>
class Emissions;

template<typename Value, typename Time>
class DICE {
  protected:
    const settings::SettingsNode& settings;
    const Global<Value, Time> global;

    std::vector<Economy<Value, Time>> economies;
    std::unique_ptr<climate::Climate<Value, Time>> climate;
    std::unique_ptr<damage::Damage<Value, Time>> damage;

    void write_netcdf_output(const settings::SettingsNode& output_node);
    void optimize(const settings::SettingsNode& optimization_node, TimeSeries<Value>& initial_values);

  public:
    DICE(const settings::SettingsNode& settings_p);
    Control<Value, Time> control;
    inline Value calc_single_utility();
    Emissions<Value, Time> emissions;
    void reset();
    void initialize();
    void output();
    void run();
};
}

#endif
