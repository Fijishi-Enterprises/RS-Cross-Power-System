# GridCal
# Copyright (C) 2015 - 2023 Santiago Peñate Vera
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 3 of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

"""
This file implements a DC-OPF for time series
That means that solves the OPF problem for a complete time series at once
"""
import numpy as np
from typing import List, Union
import ortools.linear_solver.pywraplp as ort

from GridCal.Engine.basic_structures import ZonalGrouping
from GridCal.Engine.basic_structures import MIPSolvers
from GridCal.Engine.Core.Devices.multi_circuit import MultiCircuit
from GridCal.Engine.Core.DataStructures.numerical_circuit import NumericalCircuit, compile_numerical_circuit_at
from GridCal.Engine.Core.DataStructures.generator_data import GeneratorData
from GridCal.Engine.Core.DataStructures.battery_data import BatteryData
from GridCal.Engine.Core.DataStructures.load_data import LoadData
from GridCal.Engine.Core.DataStructures.branch_data import BranchData
from GridCal.Engine.Core.DataStructures.hvdc_data import HvdcData
from GridCal.Engine.Core.DataStructures.bus_data import BusData
from GridCal.Engine.basic_structures import Logger, Mat, Vec, IntVec, DateVec, BoolVec
import GridCal.ThirdParty.ortools.ortools_extra as pl
from GridCal.Engine.Core.Devices.enumerations import TransformerControlType, HvdcControlType
from GridCal.Engine.Simulations.LinearFactors.linear_analysis import LinearAnalysis, LinearMultiContingency
from GridCal.Engine.Simulations.OPF.linear_opf_ts import BatteryVars, BranchVars, BusVars, GenerationVars, HvdcVars, LoadVars, OpfVars, get_contingency_flow_with_filter

def join(init: str, vals: List[int], sep="_"):
    """
    Generate naming string
    :param init: initial string
    :param vals: concatenation of indices
    :param sep: separator
    :return: naming string
    """
    return init + sep.join([str(x) for x in vals])

def formulate_monitorization_logic(
        branch_data_t: BranchData,
        base_flows: Vec,
        rates:Vec,
        max_alpha: float, # todo: check if consider alpha and alpha_n1 in separate vars
        branch_sensitivity_threshold: float,
        structural_ntc: float,
        ntc_load_rule: float,
        only_sensitive_branches: bool,
        only_ntc_load_rule_branches: bool,
):
    """
    Function to formulate branch monitor status due the given logic
    :param branch_data_t: BranchData
    :param base_flows: branch base flows
    :param rates: array of branch rates
    :param max_alpha: Array of max absolute branch sensitivity to the exchange in n and n-1 condition
    :param branch_sensitivity_threshold: branch sensitivity to the exchange threshold
    :param structural_ntc: Maximum NTC available by thermal interconnection rates.
    :param ntc_load_rule: percentage of loading reserved to exchange flow (Clean Energy Package rule by ACER).
    :param only_sensitive_branches: boolean to apply sensitivity threshold to the monitor logic.
    :param only_ntc_load_rule_branches: boolean to apply ntc load rule to the monitor logic.
    return:
        - monitor: monitor logic final result
        - monitor_loading: monitor logic per branch set by user interface
        - monitor_by_sensitivity: monitor logic due exchange sensibility
        - monitor_by_unrealistic_ntc: monitor logic due unrealistic minimum ntc
        - monitor_by_zero_exchange: monitor logic due zero exchange loading
        - branch_ntc_load_rule: branch minimum ntc to be considered as limiting element
        - branch_zero_exchange_loading: branch loading for zero exchange situation.
    """

    # NTC min for considering as limiting element by CEP rule
    branch_ntc_load_rule = ntc_load_rule * rates / (max_alpha + 1e-20)

    # Branch load without exchange
    branch_zero_exchange_loading = base_flows * (1 - max_alpha) / rates

    # Exclude Branches with not enough sensibility to exchange
    if only_sensitive_branches:
        monitor_by_sensitivity = max_alpha > branch_sensitivity_threshold
    else:
        monitor_by_sensitivity = np.ones(branch_data_t.nelm, dtype=bool)

    # Avoid unrealistic ntc && Exclude Branches with 'interchange zero' flows over CEP rule limit
    if only_ntc_load_rule_branches:
        monitor_by_unrealistic_ntc = branch_ntc_load_rule <= structural_ntc
        monitor_by_zero_exchange = branch_zero_exchange_loading >= (1 - ntc_load_rule)
    else:
        monitor_by_unrealistic_ntc = np.ones(branch_data_t.nelm, dtype=bool)
        monitor_by_zero_exchange = np.ones(branch_data_t.nelm, dtype=bool)

    monitor = branch_data_t.monitor_loading * \
              monitor_by_sensitivity * \
              monitor_by_unrealistic_ntc * \
              monitor_by_zero_exchange

    return monitor, \
           branch_data_t.monitor_loading,\
           monitor_by_sensitivity,\
           monitor_by_unrealistic_ntc, \
           monitor_by_zero_exchange, \
           branch_ntc_load_rule, \
           branch_zero_exchange_loading


def add_ntc_bus_angles_formulation(
        t: Union[int, None],
        bus_data_t: BusData,
        bus_vars: BusVars,
        prob: ort.Solver,
        vd: IntVec,
        logger: Logger = Logger()):
    """
    Add MIP bus angles formulation
    :param t: time step, if None we assume single time step
    :param bus_data_t: BusData structure
    :param bus_vars: BusVars structure
    :param prob: ORTools problem
    :param vd: array of slack bus indices
    :param logger: logger instance
    :return objective function
    """
    f_obj = 0.0

    for k in range(bus_data_t.nelm):

        if bus_data_t.angle_min[k] > bus_data_t.angle_max[k]:
            logger.add_error(
                msg='Theta min > Theta max',
                device=f'Bus {bus_data_t.names[k]}',
                value=f'{bus_data_t.angle_min[k]} > {bus_data_t.angle_max[k]}')

        else:
            bus_vars.theta[t, k] = prob.NumVar(
                lb=bus_data_t.angle_min[k],
                ub=bus_data_t.angle_max[k],
                name=join("bus_theta_", [t, k], "_"))

    for k in vd:
        bus_vars.theta[t, k].SetLb(0)
        bus_vars.theta[t, k].SetUb(0)

    return f_obj


def add_ntc_generation_formulation(
        t: Union[int, None],
        Sbase: float,
        time_array: DateVec,
        gen_data_t: GeneratorData,
        gen_vars: GenerationVars,
        prob: ort.Solver,
        unit_commitment: bool,
        ramp_constraints: bool,
        skip_generation_limits: bool):
    """
    Add MIP generation formulation
    :param t: time step
    :param Sbase: base power (100 MVA)
    :param time_array: complete time array
    :param gen_data_t: GeneratorData structure
    :param gen_vars: GenerationVars structure
    :param prob: ORTools problem
    :param unit_commitment: formulate unit commitment?
    :param ramp_constraints: formulate ramp constraints?
    :param skip_generation_limits: skip the generation limits?
    :return objective function
    """
    f_obj = 0.0
    for k in range(gen_data_t.nelm):

        if gen_data_t.active[k]:

            # declare active power var (limits will be applied later)
            gen_vars.p[t, k] = prob.NumVar(
                lb=-1e20,
                ub=1e20,
                name=join("gen_p_", [t, k], "_"))

            if gen_data_t.dispatchable[k]:

                if unit_commitment:

                    # declare unit commitment vars
                    gen_vars.starting_up[t, k] = prob.IntVar(0, 1, join("gen_starting_up_", [t, k], "_"))
                    gen_vars.producing[t, k] = prob.IntVar(0, 1, join("gen_producing_", [t, k], "_"))
                    gen_vars.shutting_down[t, k] = prob.IntVar(0, 1, join("gen_shutting_down_", [t, k], "_"))

                    # operational cost (linear...)
                    f_obj += gen_data_t.cost_1[k] * gen_vars.p[t, k] + gen_data_t.cost_0[k] * gen_vars.producing[t, k]

                    # start-up cost
                    f_obj += gen_data_t.startup_cost[k] * gen_vars.starting_up[t, k]

                    # power boundaries of the generator
                    if not skip_generation_limits:
                        prob.Add(gen_vars.p[t, k] >= (
                                gen_data_t.availability[k] * gen_data_t.pmin[k] / Sbase * gen_vars.producing[t, k]),
                                 join("gen>=Pmin", [t, k], "_"))
                        prob.Add(gen_vars.p[t, k] <= (
                                gen_data_t.availability[k] * gen_data_t.pmax[k] / Sbase * gen_vars.producing[t, k]),
                                 join("gen<=Pmax", [t, k], "_"))

                    if t is not None:
                        if t == 0:
                            prob.Add(gen_vars.starting_up[t, k] - gen_vars.shutting_down[t, k] ==
                                     gen_vars.producing[t, k] - float(gen_data_t.active[k]),
                                     join("binary_alg1_", [t, k], "_"))
                            prob.Add(gen_vars.starting_up[t, k] + gen_vars.shutting_down[t, k] <= 1,
                                     join("binary_alg2_", [t, k], "_"))
                        else:
                            prob.Add(
                                gen_vars.starting_up[t, k] - gen_vars.shutting_down[t, k] ==
                                gen_vars.producing[t, k] - gen_vars.producing[t - 1, k],
                                join("binary_alg3_", [t, k], "_"))
                            prob.Add(gen_vars.starting_up[t, k] + gen_vars.shutting_down[t, k] <= 1,
                                     join("binary_alg4_", [t, k], "_"))
                else:
                    # No unit commitment

                    # Operational cost (linear...)
                    f_obj += (gen_data_t.cost_1[k] * gen_vars.p[t, k]) + gen_data_t.cost_0[k]

                    if not skip_generation_limits:
                        gen_vars.p[t, k].SetLb(gen_data_t.availability[k] * gen_data_t.pmin[k] / Sbase)
                        gen_vars.p[t, k].SetUb(gen_data_t.availability[k] * gen_data_t.pmax[k] / Sbase)

                # add the ramp constraints
                if ramp_constraints and t is not None:
                    if t > 0:
                        if gen_data_t.ramp_up[k] < gen_data_t.pmax[k] and gen_data_t.ramp_down[k] < gen_data_t.pmax[k]:
                            # if the ramp is actually sufficiently restrictive...
                            dt = (time_array[t] - time_array[t - 1]).seconds / 3600.0  # time increment in hours

                            # - ramp_down · dt <= P(t) - P(t-1) <= ramp_up · dt
                            prob.Add(-gen_data_t.ramp_down[k] / Sbase * dt <= gen_vars.p[t, k] - gen_vars.p[t - 1, k])
                            prob.Add(gen_vars.p[t, k] - gen_vars.p[t - 1, k] <= gen_data_t.ramp_up[k] / Sbase * dt)
            else:

                # it is NOT dispatchable
                p = gen_data_t.p[k] / Sbase

                # Operational cost (linear...)
                f_obj += (gen_data_t.cost_1[k] * gen_vars.p[t, k]) + gen_data_t.cost_0[k]

                # the generator is not dispatchable at time step
                if p > 0:

                    gen_vars.shedding[t, k] = prob.NumVar(0, p, join("gen_shedding_", [t, k], "_"))

                    prob.Add(gen_vars.p[t, k] == gen_data_t.p[k] / Sbase - gen_vars.shedding[t, k],
                             join("gen==PG-PGslack", [t, k], "_"))

                    f_obj += gen_data_t.cost_1[k] * gen_vars.shedding[t, k]

                elif p < 0:
                    # the negative sign is because P is already negative here, to make it positive
                    gen_vars.shedding[t, k] = prob.NumVar(0, -p, join("gen_shedding_", [t, k], "_"))

                    prob.Add(gen_vars.p[t, k] == p + gen_vars.shedding[t, k],
                             join("gen==PG+PGslack", [t, k], "_"))

                    f_obj += gen_data_t.cost_1[k] * gen_vars.shedding[t, k]

                else:
                    # the generation value is exactly zero, pass
                    pass

                gen_vars.producing[t, k] = 1
                gen_vars.shutting_down[t, k] = 0
                gen_vars.starting_up[t, k] = 0

        else:
            # the generator is not available at time step
            gen_vars.p[t, k] = 0.0

    return f_obj


def add_ntc_battery_formulation(
        t: Union[int, None],
        Sbase: float,
        time_array: DateVec,
        batt_data_t: BatteryData,
        batt_vars: BatteryVars,
        prob: ort.Solver,
        unit_commitment: bool,
        ramp_constraints: bool,
        skip_generation_limits: bool,
        energy_0: Vec):
    """
    Add MIP generation formulation
    :param t: time step, if None we assume single time step
    :param Sbase: base power (100 MVA)
    :param time_array: complete time array
    :param batt_data_t: BatteryData structure
    :param batt_vars: BatteryVars structure
    :param prob: ORTools problem
    :param unit_commitment: formulate unit commitment?
    :param ramp_constraints: formulate ramp constraints?
    :param skip_generation_limits: skip the generation limits?
    :param energy_0: initial value of the energy stored
    :return objective function
    """
    f_obj = 0.0
    for k in range(batt_data_t.nelm):

        if batt_data_t.active[k]:

            # declare active power var (limits will be applied later)
            batt_vars.p[t, k] = prob.NumVar(0, 1e20, join("batt_p_", [t, k], "_"))

            if batt_data_t.dispatchable[k]:

                if unit_commitment:

                    # declare unit commitment vars
                    batt_vars.starting_up[t, k] = prob.IntVar(0, 1, join("bat_starting_up_", [t, k], "_"))
                    batt_vars.producing[t, k] = prob.IntVar(0, 1, join("bat_producing_", [t, k], "_"))
                    batt_vars.shutting_down[t, k] = prob.IntVar(0, 1, join("bat_shutting_down_", [t, k], "_"))

                    # operational cost (linear...)
                    f_obj += batt_data_t.cost_1[k] * batt_vars.p[t, k] + batt_data_t.cost_0[k] * batt_vars.producing[t, k]

                    # start-up cost
                    f_obj += batt_data_t.startup_cost[k] * batt_vars.starting_up[t, k]

                    # power boundaries of the generator
                    if not skip_generation_limits:
                        prob.Add(batt_vars.p[t, k] >= (batt_data_t.availability[k] * batt_data_t.pmin[k] / Sbase * batt_vars.producing[t, k]),
                                 join("batt>=Pmin", [t, k], "_"))
                        prob.Add(batt_vars.p[t, k] <= (batt_data_t.availability[k] * batt_data_t.pmax[k] / Sbase * batt_vars.producing[t, k]),
                                 join("batt<=Pmax", [t, k], "_"))

                    if t is not None:
                        if t == 0:
                            prob.Add(batt_vars.starting_up[t, k] - batt_vars.shutting_down[t, k] == batt_vars.producing[t, k] - float(batt_data_t.active[k]),
                                     join("binary_alg1_", [t, k], "_"))
                            prob.Add(batt_vars.starting_up[t, k] + batt_vars.shutting_down[t, k] <= 1,
                                     join("binary_alg2_", [t, k], "_"))
                        else:
                            prob.Add(batt_vars.starting_up[t, k] - batt_vars.shutting_down[t, k] == batt_vars.producing[t, k] - batt_vars.producing[t - 1, k],
                                     join("binary_alg3_", [t, k], "_"))
                            prob.Add(batt_vars.starting_up[t, k] + batt_vars.shutting_down[t, k] <= 1,
                                     join("binary_alg4_", [t, k], "_"))
                else:
                    # No unit commitment

                    # Operational cost (linear...)
                    f_obj += (batt_data_t.cost_1[k] * batt_vars.p[t, k]) + batt_data_t.cost_0[k]

                    # power boundaries of the generator
                    if not skip_generation_limits:
                        batt_vars.p[t, k].SetLb(batt_data_t.availability[k] * batt_data_t.pmin[k] / Sbase)
                        batt_vars.p[t, k].SetUb(batt_data_t.availability[k] * batt_data_t.pmax[k] / Sbase)

                # compute the time increment in hours
                dt = (time_array[t] - time_array[t - 1]).seconds / 3600.0

                if ramp_constraints and t is not None:
                    if t > 0:

                        # add the ramp constraints
                        if batt_data_t.ramp_up[k] < batt_data_t.pmax[k] and \
                                batt_data_t.ramp_down[k] < batt_data_t.pmax[k]:
                            # if the ramp is actually sufficiently restrictive...
                            # - ramp_down · dt <= P(t) - P(t-1) <= ramp_up · dt
                            prob.Add(
                                -batt_data_t.ramp_down[k] / Sbase * dt <= batt_vars.p[t, k] - batt_vars.p[t - 1, k])
                            prob.Add(
                                batt_vars.p[t, k] - batt_vars.p[t - 1, k] <= batt_data_t.ramp_up[k] / Sbase * dt)

                # set the energy  value Et = E(t - 1) + dt * Pb / eff
                batt_vars.e[t, k] = prob.NumVar(batt_data_t.e_min[k] / Sbase, batt_data_t.e_max[k] / Sbase,
                                                join("batt_e_", [t, k], "_"))

                if t > 0:
                    # energy decreases / increases with power · dt
                    prob.Add(
                        batt_vars.e[t, k] == batt_vars.e[t - 1, k] + dt * batt_data_t.efficiency[k] * batt_vars.p[t, k])
                else:
                    # set the initial energy value
                    batt_vars.e[t, k] = energy_0[k] / Sbase

            else:

                # it is NOT dispatchable

                # Operational cost (linear...)
                f_obj += (batt_data_t.cost_1[k] * batt_vars.p[t, k]) + batt_data_t.cost_0[k]

                p = batt_data_t.p[k] / Sbase

                # the generator is not dispatchable at time step
                if p > 0:

                    batt_vars.shedding[t, k] = prob.NumVar(0, p, join("bat_shedding_", [t, k], "_"))

                    prob.Add(batt_vars.p[t, k] == batt_data_t.p[k] / Sbase - batt_vars.shedding[t, k],
                             join("batt==PB-PBslack", [t, k], "_"))

                    f_obj += batt_data_t.cost_1[k] * batt_vars.shedding[t, k]

                elif p < 0:
                    # the negative sign is because P is already negative here
                    batt_vars.shedding[t, k] = prob.NumVar(0, -p, join("bat_shedding_", [t, k], "_"))

                    prob.Add(batt_vars.p[t, k] == batt_data_t.p[k] / Sbase + batt_vars.shedding[t, k],
                             join("batt==PB+PBslack", [t, k], "_"))

                    f_obj += batt_data_t.cost_1[k] * batt_vars.shedding[t, k]

                else:
                    # the generation value is exactly zero, pass
                    pass

                batt_vars.producing[t, k] = 1
                batt_vars.shutting_down[t, k] = 0
                batt_vars.starting_up[t, k] = 0

        else:
            # the generator is not available at time step
            batt_vars.p[t, k] = 0.0

    return f_obj


def add_ntc_load_formulation(
        t: Union[int, None],
        Sbase: float,
        load_data_t: LoadData,
        load_vars: LoadVars,
        prob: ort.Solver):
    """
    Add MIP generation formulation
    :param t: time step, if None we assume single time step
    :param Sbase: base power (100 MVA)
    :param load_data_t: BatteryData structure
    :param load_vars: BatteryVars structure
    :param prob: ORTools problem
    :return objective function
    """
    f_obj = 0.0
    for k in range(load_data_t.nelm):

        if load_data_t.active[k]:

            # store the load
            load_vars.p[t, k] = load_data_t.S[k].real / Sbase

            if load_vars.p[t, k] > 0.0:

                # assign load shedding variable
                load_vars.shedding[t, k] = prob.NumVar(lb=0, ub=load_vars.p[t, k],
                                                       name=join("load_shedding_", [t, k], "_"))

                # minimize the load shedding
                f_obj += load_data_t.cost[k] * load_vars.shedding[t, k]
            else:
                # the load is negative, won't shed?
                load_vars.shedding[t, k] = 0.0

        else:
            # the load is not available at time step
            load_vars.shedding[t, k] = 0.0

    return f_obj


def add_ntc_branches_formulation(
        t: int,
        Sbase: float,
        branch_data_t: BranchData,
        branch_vars: BranchVars,
        gen_vars: GenerationVars,
        gen_data_t: GeneratorData,
        vars_bus: BusVars,
        monitor: BoolVec,
        prob: ort.Solver,
        consider_contingencies: bool,
        multi_contingencies: List[LinearMultiContingency],
        lodf_threshold: float,
        inf=1e20):
    """
    Formulate the branches
    :param t: time index
    :param Sbase: base power (100 MVA)
    :param branch_data_t: BranchData
    :param branch_vars: BranchVars
    :param gen_data_t: GeneratorData
    :param gen_vars: GeneratorVars
    :param vars_bus: BusVars
    :param monitor: monitor logic vector
    :param prob: OR problem
    :param consider_contingencies:
    :param multi_contingencies: List of LinearMultiContingency object
    :param lodf_threshold: LODF threshold to cut-off
    :param inf: number considered infinte
    :return objective function
    """
    f_obj = 0.0

    assert monitor.shape == branch_data_t.monitor_loading.shape

    # for each branch
    for m in range(branch_data_t.nelm):
        fr = branch_data_t.F[m]
        to = branch_data_t.T[m]

        # copy rates
        branch_vars.rates[t, m] = branch_data_t.rates[m]

        if branch_data_t.active[m]:

            # declare the flow LPVar
            branch_vars.flows[t, m] = prob.NumVar(
                lb=-inf,
                ub=inf,
                name=join("flow_", [t, m], "_"))

            # compute the branch susceptance
            if branch_data_t.X[m] == 0.0:
                if branch_data_t.R[m] != 0.0:
                    bk = -1.0 / branch_data_t.R[m]
                else:
                    bk = 1e-20
            else:
                bk = -1.0 / branch_data_t.X[m]

            # compute the flow
            if branch_data_t.control_mode[m] == TransformerControlType.Pt:

                # add angle
                branch_vars.tap_angles[t, m] = prob.NumVar(
                    lb=branch_data_t.tap_angle_min[m],
                    ub=branch_data_t.tap_angle_max[m],
                    name=join("tap_ang_", [t, m], "_"))

                # is a phase shifter device (like phase shifter transformer or VSC with P control)
                flow_ctr = branch_vars.flows[t, m] == bk * (
                        vars_bus.theta[t, fr] -
                        vars_bus.theta[t, to] + # todo: check if + or -
                        branch_vars.tap_angles[t, m])
                prob.Add(
                    constraint=flow_ctr,
                    name=join("Branch_flow_set_with_ps_", [t, m], "_"))

                # power injected and subtracted due to the phase shift
                vars_bus.branch_injections[fr] = -bk * branch_vars.tap_angles[t, m]
                vars_bus.branch_injections[to] = bk * branch_vars.tap_angles[t, m]

            else:  # rest of the branches
                flow_ctr = branch_vars.flows[t, m] == bk * (
                        vars_bus.theta[t, fr] -
                        vars_bus.theta[t, to])
                prob.Add(
                    constraint=flow_ctr,
                    name=join("Branch_flow_set_", [t, m], "_"))

            # add the flow constraint if monitored
            if monitor[m]:
                # TODO: repasar que hacer con los flow slacks
                # branch_vars.flow_slacks_pos[t, m] = prob.NumVar(
                #     lb=0,
                #     ub=inf,
                #     name=join("flow_slack_pos_", [t, m], "_"))
                # branch_vars.flow_slacks_neg[t, m] = prob.NumVar(
                #     lb=0,
                #     ub=inf,
                #     name=join("flow_slack_neg_", [t, m], "_"))

                # add upper rate constraint
                # branch_vars.flow_constraints_ub[t, m] = branch_vars.flows[t, m] + branch_vars.flow_slacks_pos[t, m] - branch_vars.flow_slacks_neg[t, m] <= branch_data_t.rates[m] / Sbase
                # prob.Add(
                #     constraint=branch_vars.flow_constraints_ub[t, m],
                #     name=join("br_flow_upper_lim_", [t, m]))
                #
                # # add lower rate constraint
                # branch_vars.flow_constraints_lb[t, m] = branch_vars.flows[t, m] + branch_vars.flow_slacks_pos[t, m] - branch_vars.flow_slacks_neg[t, m] >= -branch_data_t.rates[m] / Sbase
                # prob.Add(
                #     constraint=branch_vars.flow_constraints_lb[t, m],
                #     name=join("br_flow_lower_lim_", [t, m]))
                #
                # # add to the objective function
                # f_obj += branch_data_t.overload_cost[m] * branch_vars.flow_slacks_pos[t, m]
                # f_obj += branch_data_t.overload_cost[m] * branch_vars.flow_slacks_neg[t, m]

                # todo: de momento para mantener lo del mou yo haría esto otro:
                # Redefine branch limits
                branch_vars.flows[t, m].setUb(branch_data_t.rates[m] / Sbase)
                branch_vars.flows[t, m].setLb(-branch_data_t.rates[m] / Sbase)

                if consider_contingencies:

                    for c, multi_ctg in enumerate(multi_contingencies):

                        # compute contingency injections
                        if len(multi_ctg.bus_indices):
                            injections = pl.lpDot(gen_data_t.C_bus_elm, gen_vars.p[t, :])
                        else:
                            injections = None

                        flow_n1 = get_contingency_flow_with_filter(
                            multi_contingency=multi_ctg,
                            base_flow=branch_vars.flows[t, :],
                            injections=injections,
                            threshold=lodf_threshold,
                            m=m)

                        flow_n1_var = prob.NumVar(
                            ub=inf,
                            lb=-inf,
                            name=join("flown1_", [t, m, c], "_"))

                        prob.Add(
                            constraint=flow_n1_var == flow_n1,
                            name=join("flow_n1_ctr", [t, m, c], "_"))

                        branch_vars.add_contingency_flow(m=m, c=c, flow_var=flow_n1_var)

    return f_obj


def add_ntc_hvdc_formulation(
        t: int,
        Sbase: float,
        hvdc_data_t: HvdcData,
        hvdc_vars: HvdcVars,
        vars_bus: BusVars,
        prob: ort.Solver):
    """

    :param t:
    :param Sbase:
    :param hvdc_data_t:
    :param hvdc_vars:
    :param vars_bus:
    :param prob:
    :return:
    """
    f_obj = 0.0
    for m in range(hvdc_data_t.nelm):

        fr = hvdc_data_t.F[m]
        to = hvdc_data_t.T[m]

        if hvdc_data_t.active[m]:

            # declare the flow var
            hvdc_vars.flows[t, m] = prob.NumVar(
                lb=-hvdc_data_t.rate[m] / Sbase,
                ub=hvdc_data_t.rate[m] / Sbase,
                name=join("hvdc_flow_", [t, m], "_"))

            if hvdc_data_t.control_mode[m] == HvdcControlType.type_0_free:

                # set the flow based on the angular difference
                P0 = hvdc_data_t.Pset[m] / Sbase
                prob.Add(
                    constraint=hvdc_vars.flows[m, t] == P0 + hvdc_data_t.angle_droop[m] * (
                            vars_bus.theta[t, fr] - vars_bus.theta[t, to]),
                    name=join("hvdc_flow_cst_", [t, m], "_"))

                # add the injections matching the flow
                vars_bus.branch_injections[fr] -= hvdc_vars.flows[t, m]
                vars_bus.branch_injections[to] += hvdc_vars.flows[t, m]

            elif hvdc_data_t.control_mode[m] == HvdcControlType.type_1_Pset:

                if hvdc_data_t.dispatchable[m]:

                    # add the injections matching the flow
                    vars_bus.branch_injections[fr] -= hvdc_vars.flows[t, m]
                    vars_bus.branch_injections[to] += hvdc_vars.flows[t, m]

                else:

                    if hvdc_data_t.Pset[m] > hvdc_data_t.rate[m]:
                        P0 = hvdc_data_t.rate[m] / Sbase
                    elif hvdc_data_t.Pset[m] < -hvdc_data_t.rate[m]:
                        P0 = -hvdc_data_t.rate[m] / Sbase
                    else:
                        P0 = hvdc_data_t.Pset[m] / Sbase

                    hvdc_vars.flows[t, m].SetBounds(P0, P0)  # make the flow equal to P0

                    # add the injections matching the flow
                    vars_bus.branch_injections[fr] -= hvdc_vars.flows[t, m]
                    vars_bus.branch_injections[to] += hvdc_vars.flows[t, m]
            else:
                raise Exception('OPF: Unknown HVDC control mode {}'.format(hvdc_data_t.control_mode[m]))
        else:
            # not active, therefore the flow is exactly zero
            hvdc_vars.flows[t, m].SetBounds(0, 0)

    return f_obj


def add_ntc_node_balance(
        t_idx: int,
        Bbus,
        vd: IntVec,
        bus_data: BusData,
        generator_data: GeneratorData,
        battery_data: BatteryData,
        load_data: LoadData,
        bus_vars: BusVars,
        gen_vars: GenerationVars,
        batt_vars: BatteryVars,
        load_vars: LoadVars,
        prob: ort.Solver):
    """
    Add the kirchoff nodal equality
    :param t_idx: time step
    :param Bbus: susceptance matrix (complete)
    :param bus_data: BusData
    :param generator_data: GeneratorData
    :param battery_data: BatteryData
    :param load_data: LoadData
    :param bus_vars: BusVars
    :param gen_vars: GenerationVars
    :param batt_vars: BatteryVars
    :param load_vars: LoadVars
    :param prob: ort.Solver
    """
    B = Bbus.tocsc()

    P_esp = bus_vars.branch_injections[t_idx, :]
    P_esp += pl.lpDot(generator_data.C_bus_elm.tocsc(),
                      gen_vars.p[t_idx, :] - gen_vars.shedding[t_idx, :])
    P_esp += pl.lpDot(battery_data.C_bus_elm.tocsc(),
                      batt_vars.p[t_idx, :] - batt_vars.shedding[t_idx, :])
    P_esp += pl.lpDot(load_data.C_bus_elm.tocsc(),
                      load_vars.shedding[t_idx, :] - load_vars.p[t_idx, :])

    # calculate the linear nodal inyection
    P_calc = pl.lpDot(B, bus_vars.theta[t_idx, :])

    # add the equality restrictions
    for k in range(bus_data.nbus):
        bus_vars.kirchhoff[t_idx, k] = prob.Add(P_calc[k] == P_esp[k], name=join("kirchoff_", [t_idx, k], "_"))

    for i in vd:
        bus_vars.theta[t_idx, i].SetBounds(0, 0)


def run_ntc_opf_ts(grid: MultiCircuit,
                   time_indices: IntVec,
                   solver_type: MIPSolvers = MIPSolvers.CBC,
                   zonal_grouping: ZonalGrouping = ZonalGrouping.NoGrouping,
                   skip_generation_limits: bool = False,
                   consider_contingencies: bool = False,
                   unit_Commitment: bool = False,
                   ramp_constraints: bool = False,
                   lodf_threshold: float = 0.001,
                   maximize_inter_area_flow: bool = False,
                   buses_areas_1=None,
                   buses_areas_2=None,
                   energy_0: Union[Vec, None] = None,
                   logger: Logger = Logger(),
                   progress_text=None,
                   progress_func=None) -> OpfVars:
    """

    :param grid: MultiCircuit instance
    :param time_indices: Time indices (in the general scheme)
    :param solver_type: MIP solver to use
    :param zonal_grouping: Zonal grouping?
    :param skip_generation_limits: Skip the generation limits?
    :param consider_contingencies: Consider the contingencies?
    :param unit_Commitment: Formulate unit commitment?
    :param ramp_constraints: Formulate ramp constraints?
    :param lodf_threshold:
    :param maximize_inter_area_flow:
    :param buses_areas_1:
    :param buses_areas_2:
    :param energy_0:
    :param logger:
    :param progress_text:
    :param progress_func:
    :return:
    """
    bus_dict = {bus: i for i, bus in enumerate(grid.buses)}
    areas_dict = {elm: i for i, elm in enumerate(grid.areas)}

    if time_indices is None:
        time_indices = [None]
    else:
        if len(time_indices) > 0:
            # time indices are ok
            pass
        else:
            time_indices = [None]

    nt = len(time_indices) if len(time_indices) > 0 else 1
    n = grid.get_bus_number()
    nbr = grid.get_branch_number_wo_hvdc()
    ng = grid.get_generators_number()
    nb = grid.get_batteries_number()
    nl = grid.get_calculation_loads_number()
    n_hvdc = grid.get_hvdc_number()

    # declare structures of LP vars
    mip_vars = OpfVars(nt=nt, nbus=n, ng=ng, nb=nb, nl=nl, nbr=nbr, n_hvdc=n_hvdc)

    solver: ort.Solver = ort.Solver.CreateSolver(solver_type.value)

    if solver is None:
        print("{} is not present".format(solver_type.value))
        logger.add_error(msg="Solver is not present", value=solver_type.value)
        return mip_vars.get_values(grid.Sbase)

    # objective function
    f_obj = 0.0

    for t_idx, t in enumerate(time_indices):  # use time_indices = [None] to simulate the snapshot

        # compile the circuit at the master time index ------------------------------------------------------------
        # note: There are very little chances of simplifying this step and experience shows it is not
        #        worth the effort, so compile every time step
        nc: NumericalCircuit = compile_numerical_circuit_at(
            circuit=grid,
            t_idx=t,  # yes, this is not a bug
            bus_dict=bus_dict,
            areas_dict=areas_dict)

        if consider_contingencies:
            # if we want to include contingencies, we'll need the LODF at this time step
            ls = LinearAnalysis(
                numerical_circuit=nc,
                distributed_slack=False,
                correct_values=True)
            ls.run(with_nx=False)
            LODF = ls.LODF
        else:
            LODF = None

        # formulate the bus angles ---------------------------------------------------------------------------------
        for k in range(nc.bus_data.nbus):
            mip_vars.bus_vars.theta[t_idx, k] = solver.NumVar(
                nc.bus_data.angle_min[k],
                nc.bus_data.angle_max[k],
                join("th_", [t_idx, k], "_"))

        # formulate loads ------------------------------------------------------------------------------------------
        f_obj += add_ntc_load_formulation(
            t=t_idx,
            Sbase=nc.Sbase,
            load_data_t=nc.load_data,
            load_vars=mip_vars.load_vars,
            prob=solver)

        # formulate generation -------------------------------------------------------------------------------------
        f_obj += add_ntc_generation_formulation(
            t=t_idx,
            Sbase=nc.Sbase,
            time_array=grid.time_profile,
            gen_data_t=nc.generator_data,
            gen_vars=mip_vars.gen_vars,
            prob=solver,
            unit_commitment=unit_Commitment,
            ramp_constraints=ramp_constraints,
            skip_generation_limits=skip_generation_limits)

        # formulate batteries --------------------------------------------------------------------------------------
        if t_idx == 0 and energy_0 is None:
            # declare the initial energy of the batteries
            energy_0 = nc.battery_data.soc_0 * nc.battery_data.enom  # in MWh here

        f_obj += add_ntc_battery_formulation(
            t=t_idx,
            Sbase=nc.Sbase,
            time_array=grid.time_profile,
            batt_data_t=nc.battery_data,
            batt_vars=mip_vars.batt_vars,
            prob=solver,
            unit_commitment=unit_Commitment,
            ramp_constraints=ramp_constraints,
            skip_generation_limits=skip_generation_limits,
            energy_0=energy_0)

        # formulate hvdc -------------------------------------------------------------------------------------------
        f_obj += add_ntc_hvdc_formulation(
            t=t_idx,
            Sbase=nc.Sbase,
            hvdc_data_t=nc.hvdc_data,
            hvdc_vars=mip_vars.hvdc_vars,
            vars_bus=mip_vars.bus_vars,
            prob=solver)

        if zonal_grouping == ZonalGrouping.NoGrouping:
            # formulate branches -----------------------------------------------------------------------------------
            f_obj += add_ntc_branches_formulation(
                t=t_idx,
                Sbase=nc.Sbase,
                branch_data_t=nc.branch_data,
                branch_vars=mip_vars.branch_vars,
                vars_bus=mip_vars.bus_vars,
                prob=solver,
                consider_contingencies=consider_contingencies,
                LODF=LODF,
                lodf_threshold=lodf_threshold,
                inf=1e20)

            # formulate nodes ---------------------------------------------------------------------------------------

            add_ntc_node_balance(
                t_idx=t_idx,
                Bbus=nc.Bbus,
                vd=nc.vd,
                bus_data=nc.bus_data,
                generator_data=nc.generator_data,
                battery_data=nc.battery_data,
                load_data=nc.load_data,
                bus_vars=mip_vars.bus_vars,
                gen_vars=mip_vars.gen_vars,
                batt_vars=mip_vars.batt_vars,
                load_vars=mip_vars.load_vars,
                prob=solver)

        elif zonal_grouping == ZonalGrouping.All:
            # this is the copper plate approach
            pass

        # production equals demand ---------------------------------------------------------------------------------
        solver.Add(solver.Sum(mip_vars.gen_vars.p[t_idx, :]) + solver.Sum(mip_vars.batt_vars.p[t_idx, :]) >=
                   mip_vars.load_vars.p[t_idx, :].sum() - mip_vars.load_vars.shedding[t_idx].sum(),
                   name="satisfy_demand_at_={0}".format(t_idx))

        if progress_func is not None:
            progress_func((t_idx + 1) / nt * 100.0)

    # set the objective function
    solver.Minimize(f_obj)

    # solve
    if progress_text is not None:
        progress_text("Solving...")

    if progress_func is not None:
        progress_func(0)

    # model_str = solver.ExportModelAsLpFormat(False)
    # with open("lynn5_busopf.lp", "w") as f:
    #     f.write(model_str)

    status = solver.Solve()

    # gather the results
    if status == ort.Solver.OPTIMAL:
        print('Solution:')
        print('Objective value =', solver.Objective().Value())
        mip_vars.acceptable_solution = True
    else:
        print('The problem does not have an optimal solution.')
        mip_vars.acceptable_solution = False
        model_str = solver.ExportModelAsLpFormat(False)
        with open("lp_debug.lp", "w") as f:
            f.write(model_str)

    # extract values from the LP Vars
    vars_v = mip_vars.get_values(grid.Sbase)

    return vars_v
