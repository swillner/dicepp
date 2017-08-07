#ifndef DICE_H
#define DICE_H

#include <autodiff.h>
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

template<typename Value, typename Time, typename Constant, typename Variable>
class Emissions;

template<typename Value, typename Time>
class DICE {
  protected:
    const settings::SettingsNode& settings;
    const Global<Value, Time> global;

  public:
    Control<autodiff::Value<Value>, Time, Value, autodiff::Variable<Value>> control;

  protected:
    std::vector<Economy<autodiff::Value<Value>, Time, Value, autodiff::Variable<Value>>> economies;
    std::unique_ptr<climate::Climate<autodiff::Value<Value>, Time, Value, autodiff::Variable<Value>>> climate;
    std::unique_ptr<damage::Damage<autodiff::Value<Value>, Time, Value, autodiff::Variable<Value>>> damage;

    void write_netcdf_output(const settings::SettingsNode& output_node);
    void write_csv_output(const settings::SettingsNode& output_node);

  public:
    DICE(const settings::SettingsNode& settings_p);
    inline autodiff::Value<Value> calc_single_utility();
    Emissions<autodiff::Value<Value>, Time, Value, autodiff::Variable<Value>> emissions;
    void reset();
    void initialize();
    void output();
    void run();
};
}

#endif
