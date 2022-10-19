# GridCal
# Copyright (C) 2022 Santiago Peñate Vera
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

import numpy as np


def fill_circuit_info(grid: "MultiCircuit"):
    area_dict = {elm: i for i, elm in enumerate(grid.get_areas())}
    bus_dict = grid.get_bus_index_dict()

    area_names = [a.name for a in grid.get_areas()]
    bus_area_indices = np.array([area_dict[b.area] for b in grid.buses])

    branches = grid.get_branches_wo_hvdc()
    F = np.zeros(len(branches), dtype=int)
    T = np.zeros(len(branches), dtype=int)
    for k, elm in enumerate(branches):
        F[k] = bus_dict[elm.bus_from]
        T[k] = bus_dict[elm.bus_to]

    hvdc = grid.get_hvdc()
    hvdc_F = np.zeros(len(hvdc), dtype=int)
    hvdc_T = np.zeros(len(hvdc), dtype=int)

    for k, elm in enumerate(hvdc):
        hvdc_F[k] = bus_dict[elm.bus_from]
        hvdc_T[k] = bus_dict[elm.bus_to]

    return area_names, bus_area_indices, F, T, hvdc_F, hvdc_T


def get_inter_area_flows_snapshot(F, T, Sf, hvdc_F, hvdc_T, hvdc_Pf, area_names, bus_area_indices):
    """
    
    :param F: 
    :param T: 
    :param Sf: 
    :param hvdc_F: 
    :param hvdc_T: 
    :param hvdc_Pf: 
    :param area_names: 
    :param bus_area_indices: 
    :return: 
    """
    na = len(area_names)
    x = np.zeros((na, na), dtype=complex)

    for f, t, flow in zip(F, T, Sf):
        a1 = bus_area_indices[f]
        a2 = bus_area_indices[t]
        if a1 != a2:
            x[a1, a2] += flow
            x[a2, a1] -= flow

    for f, t, flow in zip(hvdc_F, hvdc_T, hvdc_Pf):
        a1 = bus_area_indices[f]
        a2 = bus_area_indices[t]
        if a1 != a2:
            x[a1, a2] += flow
            x[a2, a1] -= flow

    return x


def get_branch_values_per_area_snapshot(branch_values: np.ndarray, F, T, area_names, bus_area_indices):
    """

    :param branch_values:
    :param F:
    :param T:
    :param area_names:
    :param bus_area_indices:
    :return:
    """
    na = len(area_names)
    x = np.zeros((na, na), dtype=branch_values.dtype)

    for f, t, val in zip(F, T, branch_values):
        a1 = bus_area_indices[f]
        a2 = bus_area_indices[t]
        x[a1, a2] += val

    return x


def get_hvdc_values_per_area_snapshot(hvdc_values: np.ndarray, hvdc_F, hvdc_T, area_names, bus_area_indices):
    """

    :param hvdc_values:
    :param hvdc_F:
    :param hvdc_T:
    :param area_names:
    :param bus_area_indices:
    :return:
    """
    na = len(area_names)
    x = np.zeros((na, na), dtype=hvdc_values.dtype)

    for f, t, val in zip(hvdc_F, hvdc_T, hvdc_values):
        a1 = bus_area_indices[f]
        a2 = bus_area_indices[t]
        x[a1, a2] += val

    return x


def get_ordered_area_names(area_names):
    """

    :param area_names:
    :return:
    """
    na = len(area_names)
    x = [''] * (na * na)
    for i, a in enumerate(area_names):
        for j, b in enumerate(area_names):
            x[i * na + j] = a + '->' + b
    return x


def get_inter_area_flows_ts(F, T, Sf, hvdc_F, hvdc_T, hvdc_Pf, area_names, bus_area_indices, time):
    """

    :param F:
    :param T:
    :param Sf:
    :param hvdc_F:
    :param hvdc_T:
    :param hvdc_Pf:
    :param area_names:
    :param bus_area_indices:
    :param time:
    :return:
    """

    na = len(area_names)
    nt = len(time)
    x = np.zeros((nt, na * na), dtype=complex)

    for f, t, flow in zip(F, T, Sf.T):
        a1 = bus_area_indices[f]
        a2 = bus_area_indices[t]
        if a1 != a2:
            x[:, a1 * na + a2] += flow
            x[:, a2 * na + a1] -= flow

    for f, t, flow in zip(hvdc_F, hvdc_T, hvdc_Pf.T):
        a1 = bus_area_indices[f]
        a2 = bus_area_indices[t]
        if a1 != a2:
            x[:, a1 * na + a2] += flow
            x[:, a2 * na + a1] -= flow

    return x


def get_branch_values_per_area_ts(branch_values: np.ndarray, F, T, area_names, bus_area_indices, time):
    """

    :param branch_values:
    :param F:
    :param T:
    :param area_names:
    :param bus_area_indices:
    :param time:
    :return:
    """

    na = len(area_names)
    nt = len(time)
    x = np.zeros((nt, na * na), dtype=branch_values.dtype)

    for f, t, val in zip(F, T, branch_values.T):
        a1 = bus_area_indices[f]
        a2 = bus_area_indices[t]
        x[:, a1 * na + a2] += val

    return x


def get_hvdc_values_per_area_ts(hvdc_values: np.ndarray, hvdc_F, hvdc_T, area_names, bus_area_indices, time):
    """

    :param hvdc_values:
    :param hvdc_F:
    :param hvdc_T:
    :param area_names:
    :param bus_area_indices:
    :param time:
    :return:
    """
    na = len(area_names)
    nt = len(time)
    x = np.zeros((nt, na * na), dtype=hvdc_values.dtype)

    for f, t, val in zip(hvdc_F, hvdc_T, hvdc_values.T):
        a1 = bus_area_indices[f]
        a2 = bus_area_indices[t]
        x[:, a1 * na + a2] += val

    return x
