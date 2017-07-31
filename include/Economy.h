#ifndef ECONOMY_H
#define ECONOMY_H

#include <math.h>
#include <iostream>
#include "Climate.h"
#include "Control.h"
#include "Damage.h"
#include "Global.h"
#include "settingsnode.h"
#include "types.h"

namespace dice {

template<typename Value, typename Time>
class Economy {
  protected:
    const settings::SettingsNode& settings;
    const Control<Value, Time>& control;
    const Global<Value, Time>& global;
    climate::Climate<Value, Time>& climate;
    damage::Damage<Value, Time>& damage;

    const Value C_lower{settings["C_lower"].as<Value>()};
    const Value C_pc_lower{settings["C_pc_lower"].as<Value>()};
    const Value E0{settings["E0"].as<Value>()};                          // Industrial emissions 2010 (GtCO2 per year)
    const Value E_land0{settings["E_land0"].as<Value>()};                // Carbon emissions from land 2010 (GtCO2 per year)
    const Value cprice0{settings["cprice0"].as<Value>()};                // Initial base carbon price (2005$ per tCO2)
    const Value dE_land{settings["dE_land"].as<Value>()};                // Decline rate of land emissions (per 5 years)
    const Value dA{settings["dA"].as<Value>()};                          // Decline rate of TFP per 5 years
    const Value dsig{settings["dsig"].as<Value>()};                      // Decline rate of decarbonization (per period)
    const Value gA0{settings["gA0"].as<Value>()};                        // Initial growth rate for TFP per year
    const Value gback{settings["gback"].as<Value>()};                    // Initial cost decline backstop cost 5 years
    const Value gcprice{settings["gcprice"].as<Value>()};                // Growth rate of base carbon price per year
    const Value gsigma1{settings["gsigma1"].as<Value>()};                // Initial growth of sigma (per year)
    const Value lim_mu{settings["lim_mu"].as<Value>()};                  // Upper limit on control rate after 2150
    const Value mu0{settings["mu0"].as<Value>()};                        // Initial emissions control rate for base case 2010
    const Value partfract2010{settings["partfract2010"].as<Value>()};    // Fraction of emissions under control in 2010
    const Value partfractfull{settings["partfractfull"].as<Value>()};    // Fraction of emissions under control at full time
    const Value pback{settings["pback"].as<Value>()};                    // Cost of backstop 2005$ per tCO2 2010
    const Value periodfullpart{settings["periodfullpart"].as<Value>()};  // Period at which have full participation
    const Value pop_adj{settings["pop_adj"].as<Value>()};                // Growth rate to calibrate to 2050 pop projection
    const Value pop_asym{settings["pop_asym"].as<Value>()};              // Asymptotic population (millions)
    const Value Q0{settings["Q0"].as<Value>()};                          // Initial gross output (trill 2005 USD)
    const Value tnopol{settings["tnopol"].as<Value>()};                  // Period before which no emissions controls base

    StepwiseBackwardLookingTimeSeries<Value, Time> L_series{global.timestep_num, settings["L0"].as<Value>()};
    StepwiseBackwardLookingTimeSeries<Value, Time> A_series{global.timestep_num, settings["A0"].as<Value>()};
    StepwiseBackwardLookingTimeSeries<LowerBounded<Value>, Time> K_series{global.timestep_num, {settings["K0"].as<Value>(), settings["K_lower"].as<Value>()}};
    StepwiseBackwardLookingTimeSeries<Value, Time> cca_series{global.timestep_num, settings["cca0"].as<Value>()};

  public:
    Economy(const settings::SettingsNode& settings_p,
            const Global<Value, Time>& global_p,
            const Control<Value, Time>& control_p,
            climate::Climate<Value, Time>& climate_p,
            damage::Damage<Value, Time>& damage_p)
        : settings(settings_p), global(global_p), control(control_p), climate(climate_p), damage(damage_p) {
    }

    void reset() {
        L_series.reset();
        A_series.reset();
        K_series.reset();
        cca_series.reset();
    }

    bool observe(Observer<Value, Time>& observer) {
        OBSERVE_VAR(A);
        OBSERVE_VAR(C);
        OBSERVE_VAR(C_pc);
        OBSERVE_VAR(E);
        OBSERVE_VAR(E_ind);
        OBSERVE_VAR(E_tree);
        OBSERVE_VAR(I);
        OBSERVE_VAR(K);
        OBSERVE_VAR(L);
        OBSERVE_VAR(Y);
        OBSERVE_VAR(Y_gross);
        OBSERVE_VAR(Y_net);
        OBSERVE_VAR(abatecost);
        OBSERVE_VAR(cca);
        OBSERVE_VAR(cost1);
        OBSERVE_VAR(cprice);
        OBSERVE_VAR(cpricebase);
        OBSERVE_VAR(damages);
        OBSERVE_VAR(mcabate);
        OBSERVE_VAR(partfract);
        OBSERVE_VAR(pbacktime);
        OBSERVE_VAR(periodu);
        OBSERVE_VAR(ri);
        OBSERVE_VAR(rr);
        OBSERVE_VAR(sigma);
        OBSERVE_VAR(utility);
        return true;
    }

    // Population (millions)
    Value L(Time t) {
        return L_series.get(t, [this](Time t, Value L_last) {

            return std::pow(L_series.initial_value, std::pow(1 - pop_adj, global.timestep_length * 0.2 * t)) * pop_asym
                   * std::pow(pop_asym, -std::pow(1 - pop_adj, global.timestep_length * 0.2 * t));

        });
    }

    // Level of total factor productivity
    Value A(Time t) {
        return A_series.get(t, [this](Time t, Value A_last) {

            const Value gA_t_m1 = gA0 * std::exp(-dA * global.timestep_length * (t - 1));
            return A_last / (1 - gA_t_m1);

        });
    }

    // Capital stock (trillions 2005 US dollars)
    Value K(Time t) {
        return K_series.get(t, [this](Time t, Value K_last) {

            return std::pow(1 - global.dK, global.timestep_length) * K_last + global.timestep_length * I(t - 1);

        });
    }

    // Cumulative industrial carbon emissions (GTC)
    Value cca(Time t) {
        return cca_series.get(t, [this](Time t, Value cca_last) {

            return cca_last + global.timestep_length * E_ind(t - 1) / 3.666;

        });
    }

    // CO2-equivalent-emissions output ratio
    Value sigma(Time t) {
        const Value sig0 = E0 / (Q0 * (1 - mu0));  // Carbon intensity 2010 (kgCO2 per output 2005 USD 2010)
        return sig0 * std::exp(5 / global.timestep_length * gsigma1 * (1 - std::pow(1 + dsig, t)) / (1 - std::pow(1 + dsig, 5 / global.timestep_length)));
    }

    // Average utility social discount rate
    Value rr(Time t) {
        return 1 / std::pow(1 + global.prstp, t);
    }

    // Emissions from deforestation
    Value E_tree(Time t) {
        return E_land0 * std::pow(1 - dE_land, 0.2 * global.timestep_length * t);
    }

    // Adjusted cost for backstop
    Value cost1(Time t) {
        const Value pbacktime_t = pback * std::pow(1 - gback, 0.2 * global.timestep_length * t);
        return pbacktime_t * sigma(t) / global.expcost2 / 1000;
    }

    // Fraction of emissions in control regime
    Value partfract(Time t) {
        if (t > periodfullpart) {
            return partfractfull;
        } else {
            return partfract2010 + (partfractfull - partfract2010) * t / periodfullpart;
        }
    }

    // Base Case Carbon Price
    Value cpricebase(Time t) {
        return cprice0 * std::pow(1 + gcprice, global.timestep_length * t);
    }

    // Backstop price
    Value pbacktime(Time t) {
        return pback * std::pow(1 - gback, 0.2 * global.timestep_length * t);
    }

    // Industrial emissions (GtCO2 per year)
    Value E_ind(Time t) {
        return sigma(t) * Y_gross(t) * (1 - control.mu[t]);
    }

    // Investment (trillions 2005 USD per year)
    Value I(Time t) {
        return control.s[t] * Y(t);
    }

    // Consumption (trillions 2005 US dollars per year)
    Value C(Time t) {
        return std::max(C_lower, Y(t) - I(t));
    }

    // Per capita consumption (thousands 2005 USD per year)
    Value C_pc(Time t) {
        return std::max(C_pc_lower, 1000 * C(t) / L(t));
    }

    // Real interest rate (per annum)
    Value ri(Time t) {
        if (t < global.timestep_num - 1) {
            return (1 + global.prstp) * std::pow(C_pc(t + 1) / C_pc(t), global.elasmu) - 1;
        } else {
            return 0.0;
        }
    }

    // Gross world product net of abatement and damages (trillions 2005 USD per year)
    Value Y(Time t) {
        return Y_net(t) - abatecost(t);
    }

    // Gross world product GROSS of abatement and damages (trillions 2005 USD per year)
    Value Y_gross(Time t) {
        return A(t) * std::pow(L(t) / 1000, 1 - global.gamma) * std::pow(K(t), global.gamma);
    }

    // Output net of damages equation (trillions 2005 USD per year)
    Value Y_net(Time t) {
        return Y_gross(t) * (1 - damage.damfrac(t));
    }

    // Damages (trillions 2005 USD per year)
    Value damages(Time t) {
        return Y_gross(t) * damage.damfrac(t);
    }

    // Cost of emissions reductions  (trillions 2005 USD per year)
    Value abatecost(Time t) {
        return Y_gross(t) * cost1(t) * std::pow(control.mu[t], global.expcost2) * std::pow(partfract(t), 1 - global.expcost2);
    }

    // Marginal cost of abatement (2005$ per ton CO2)
    Value mcabate(Time t) {
        return pbacktime(t) * std::pow(control.mu[t], global.expcost2 - 1);
    }

    // Carbon price (2005$ per ton of CO2)
    Value cprice(Time t) {
        return pbacktime(t) * std::pow(control.mu[t] / partfract(t), global.expcost2 - 1);
    }

    // One period utility function
    Value periodu(Time t) {
        return (std::pow(C(t) * 1000 / L(t), 1 - global.elasmu) - 1) / (1 - global.elasmu) - 1;
    }

    // Total CO2 emissions (GtCO2 per year)
    Value E(Time t) {
        return E_ind(t) + E_tree(t);
    }

    Value utility(Time t) {
        return periodu(t) * L(t) * rr(t);
    }
};
}

#endif
