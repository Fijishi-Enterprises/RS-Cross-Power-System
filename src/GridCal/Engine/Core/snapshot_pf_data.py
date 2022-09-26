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
import numba as nb
import pandas as pd
import scipy.sparse as sp
from typing import List, Tuple

from GridCal.Engine.basic_structures import Logger
from GridCal.Engine.Core.multi_circuit import MultiCircuit
from GridCal.Engine.basic_structures import BranchImpedanceMode
import GridCal.Engine.Core.topology as tp
from GridCal.Engine.Simulations.PowerFlow.NumericalMethods.ac_jacobian import Jacobian
from GridCal.Engine.Simulations.PowerFlow.NumericalMethods.acdc_jacobian import fubm_jacobian
from GridCal.Engine.Core.common_functions import compile_types
from GridCal.Engine.Simulations.sparse_solve import get_sparse_type
import GridCal.Engine.Core.Compilers.circuit_to_data as gc_compiler
import GridCal.Engine.Core.admittance_matrices as ycalc
from GridCal.Engine.Devices.enumerations import TransformerControlType, ConverterControlType

sparse_type = get_sparse_type()


@nb.njit(cache=True)
def compose_generator_voltage_profile(nbus, ntime,
                                      gen_bus_indices, gen_vset, gen_status, gen_is_controlled,
                                      bat_bus_indices, bat_vset, bat_status, bat_is_controlled,
                                      hvdc_bus_f, hvdc_bus_t, hvdc_status, hvdc_vf, hvdc_vt,
                                      iBeqv, iVtma, VfBeqbus, Vtmabus, branch_status, br_vf, br_vt):
    """
    Get the array of voltage set points per bus
    :param nbus: number of buses
    :param ntime: number of time steps
    :param gen_bus_indices: array of bus indices per generator
    :param gen_vset: array of voltage set points (ngen, ntime)
    :param gen_status: array of generator status (ngen, ntime)
    :param gen_is_controlled: array of values indicating if a generator controls the voltage or not (ngen)
    :param bat_bus_indices:  array of bus indices per battery
    :param bat_vset: array of voltage set points (nbat, ntime)
    :param bat_status: array of battery status (nbat, ntime)
    :param bat_is_controlled: array of values indicating if a battery controls the voltage or not (ngen)
    :param hvdc_bus_f:
    :param hvdc_bus_t:
    :param hvdc_status:
    :param hvdc_vf
    :param hvdc_vt
    :param iBeqv: indices of the branches when controlling Vf with Beq
    :param iVtma: indices of the branches when controlling Vt with ma
    :param VfBeqbus: indices of the buses where Vf is controlled by Beq
    :param Vtmabus: indices of the buses where Vt is controlled by ma
    :param branch_status:
    :param br_vf:
    :param br_vt:
    :return: Voltage set points array per bus (nbus, ntime)
    """
    V = np.ones((nbus, ntime), dtype=nb.complex128)
    used = np.zeros((nbus, ntime), dtype=nb.int8)

    # generators
    for i, bus_idx in enumerate(gen_bus_indices):
        if gen_is_controlled[i]:
            for t in range(ntime):
                if used[bus_idx, t] == 0:
                    if gen_status[i, t]:
                        V[bus_idx, t] = complex(gen_vset[i, t], 0)
                        used[bus_idx, t] = 1

    # batteries
    for i, bus_idx in enumerate(bat_bus_indices):
        if bat_is_controlled[i]:
            for t in range(ntime):
                if used[bus_idx, t] == 0:
                    if bat_status[i, t]:
                        V[bus_idx, t] = complex(bat_vset[i, t], 0)
                        used[bus_idx, t] = 1

    # HVDC
    for i in range(hvdc_status.shape[0]):
        from_idx = hvdc_bus_f[i]
        to_idx = hvdc_bus_t[i]
        for t in range(ntime):
            if hvdc_status[i, t] != 0:
                if used[from_idx, t] == 0:
                    V[from_idx, t] = complex(hvdc_vf[i, t], 0)
                    used[from_idx, t] = 1
                if used[to_idx, t] == 0:
                    V[to_idx, t] = complex(hvdc_vt[i, t], 0)
                    used[to_idx, t] = 1

    # branch - from
    for i in iBeqv:  # branches controlling Vf
        from_idx = VfBeqbus[i]
        for t in range(ntime):
            if branch_status[i, t] != 0:
                if used[from_idx, t] == 0:
                    V[from_idx, t] = complex(br_vf[i, t], 0)
                    used[from_idx, t] = 1

    # branch - to
    for i in iVtma:  # branches controlling Vt
        from_idx = Vtmabus[i]
        for t in range(ntime):
            if branch_status[i, t] != 0:
                if used[from_idx, t] == 0:
                    V[from_idx, t] = complex(br_vt[i, t], 0)
                    used[from_idx, t] = 1

    return V


def get_inter_areas_branch(F, T, buses_areas_1, buses_areas_2):
    """
    Get the branches that join two areas
    :param buses_areas_1: Area from
    :param buses_areas_2: Area to
    :return: List of (branch index, branch object, flow sense w.r.t the area exchange)
    """
    nbr = len(F)
    lst: List[Tuple[int, float]] = list()
    for k in range(nbr):
        if F[k] in buses_areas_1 and T[k] in buses_areas_2:
            lst.append((k, 1.0))
        elif F[k] in buses_areas_2 and T[k] in buses_areas_1:
            lst.append((k, -1.0))
    return lst


def get_devices_per_areas(Cdev, buses_in_a1, buses_in_a2):
    """
    Get the generators that belong to the Area 1, Area 2 and the rest of areas
    :param buses_in_a1: List of bus indices of the area 1
    :param buses_in_a2: List of bus indices of the area 2
    :return: Tree lists: (gens_in_a1, gens_in_a2, gens_out) each of the lists contains (bus index, generator index) tuples
    """
    gens_in_a1 = list()
    gens_in_a2 = list()
    gens_out = list()
    for j in range(Cdev.shape[1]):  # for each bus
        for ii in range(Cdev.indptr[j], Cdev.indptr[j + 1]):
            i = Cdev.indices[ii]
            if i in buses_in_a1:
                gens_in_a1.append((i, j))  # i: bus idx, j: gen idx
            elif i in buses_in_a2:
                gens_in_a2.append((i, j))  # i: bus idx, j: gen idx
            else:
                gens_out.append((i, j))  # i: bus idx, j: gen idx

    return gens_in_a1, gens_in_a2, gens_out


class SnapshotData:

    def __init__(self, nbus, nline, ndcline, ntr, nvsc, nupfc, nhvdc, nload, ngen, nbatt, nshunt, nstagen, sbase, ntime=1):
        """

        :param nbus:
        :param nline:
        :param ndcline:
        :param ntr:
        :param nvsc:
        :param nhvdc:
        :param nload:
        :param ngen:
        :param nbatt:
        :param nshunt:
        :param nstagen:
        :param sbase:
        """

        self.nbus = nbus
        self.nline = nline
        self.ndcline = ndcline
        self.ntr = ntr
        self.nvsc = nvsc
        self.nupfc = nupfc
        self.nhvdc = nhvdc

        self.nload = nload
        self.ngen = ngen
        self.nbatt = nbatt
        self.nshunt = nshunt
        self.nstagen = nstagen
        self.ntime = ntime
        self.nbr = nline + ntr + nvsc + ndcline

        self.Sbase = sbase

        self.any_control = False
        self.iPfsh = list()  # indices of the branches controlling Pf flow with theta sh
        self.iQfma = list()  # indices of the branches controlling Qf with ma
        self.iBeqz = list()  # indices of the branches when forcing the Qf flow to zero (aka "the zero condition")
        self.iBeqv = list()  # indices of the branches when controlling Vf with Beq
        self.iVtma = list()  # indices of the branches when controlling Vt with ma
        self.iQtma = list()  # indices of the branches controlling the Qt flow with ma
        self.iPfdp = list()  # indices of the drop-Vm converters controlling the power flow with theta sh
        self.iPfdp_va = list()  # indices of the drop-Va converters controlling the power flow with theta sh
        self.iVscL = list()  # indices of the converters
        self.VfBeqbus = list()  # indices of the buses where Vf is controlled by Beq
        self.Vtmabus = list()  # indices of the buses where Vt is controlled by ma

        # --------------------------------------------------------------------------------------------------------------
        # Data structures
        # --------------------------------------------------------------------------------------------------------------
        self.bus_data = gc_compiler.BusData(nbus=nbus, ntime=ntime)
        self.branch_data = gc_compiler.BranchData(nbr=self.nbr, nbus=nbus, ntime=ntime)
        self.line_data = gc_compiler.LinesData(nline=nline, nbus=nbus, ntime=ntime)
        self.dc_line_data = gc_compiler.DcLinesData(ndcline=ndcline, nbus=nbus, ntime=ntime)
        self.transformer_data = gc_compiler.TransformerData(ntr=ntr, nbus=nbus, ntime=ntime)
        self.vsc_data = gc_compiler.VscData(nvsc=nvsc, nbus=nbus, ntime=ntime)
        self.upfc_data = gc_compiler.UpfcData(nelm=nupfc, nbus=nbus, ntime=ntime)
        self.hvdc_data = gc_compiler.HvdcData(nhvdc=nhvdc, nbus=nbus, ntime=ntime)

        self.load_data = gc_compiler.LoadData(nload=nload, nbus=nbus, ntime=ntime)
        self.static_generator_data = gc_compiler.StaticGeneratorData(nstagen=nstagen, nbus=nbus, ntime=ntime)
        self.battery_data = gc_compiler.BatteryData(nbatt=nbatt, nbus=nbus, ntime=ntime)
        self.generator_data = gc_compiler.GeneratorData(ngen=ngen, nbus=nbus, ntime=ntime)
        self.shunt_data = gc_compiler.ShuntData(nshunt=nshunt, nbus=nbus, ntime=ntime)

        self.original_bus_idx = np.arange(self.nbus)
        self.original_branch_idx = np.arange(self.nbr)
        self.original_line_idx = np.arange(self.nline)
        self.original_tr_idx = np.arange(self.ntr)
        self.original_dc_line_idx = np.arange(self.ndcline)
        self.original_vsc_idx = np.arange(self.nvsc)
        self.original_upfc_idx = np.arange(self.nupfc)
        self.original_hvdc_idx = np.arange(self.nhvdc)
        self.original_gen_idx = np.arange(self.ngen)
        self.original_bat_idx = np.arange(self.nbatt)
        self.original_load_idx = np.arange(self.nload)
        self.original_stagen_idx = np.arange(self.nstagen)
        self.original_shunt_idx = np.arange(self.nshunt)
        self.original_time_idx = np.arange(self.ntime)

        # --------------------------------------------------------------------------------------------------------------
        # Internal variables filled on demand, to be ready to consume once computed
        # --------------------------------------------------------------------------------------------------------------

        self.Cf_ = None
        self.Ct_ = None
        self.A_ = None

        self.Vbus_ = None
        self.Sbus_ = None
        self.Ibus_ = None
        self.YloadBus_ = None
        self.Yshunt_from_devices_ = None

        self.Qmax_bus_ = None
        self.Qmin_bus_ = None
        self.Bmax_bus_ = None
        self.Bmin_bus_ = None

        self.Admittances = None

        # Admittance for HELM / AC linear
        self.Yseries_ = None
        self.Yshunt_ = None

        # Admittances for Fast-Decoupled
        self.B1_ = None
        self.B2_ = None

        # Admittances for Linear
        self.Bbus_ = None
        self.Bf_ = None
        self.Btheta_ = None
        self.Bpqpv_ = None
        self.Bref_ = None

        self.pq_ = None
        self.pv_ = None
        self.vd_ = None
        self.pqpv_ = None
        self.ac_ = None
        self.dc_ = None

        self.available_structures = ['Vbus',
                                     'Sbus',
                                     'Ibus',
                                     'Ybus',  'G', 'B',
                                     'Yf',
                                     'Yt',
                                     'Bbus',
                                     'Bf',
                                     'Cf',
                                     'Ct',
                                     'Yshunt',
                                     'Yseries',
                                     "B'",
                                     "B''",
                                     'Types',
                                     'Jacobian',
                                     'Qmin',
                                     'Qmax',
                                     'pq',
                                     'pv',
                                     'vd',
                                     'pqpv',
                                     'tap_f',
                                     'tap_t',
                                     'original_bus_idx',
                                     'original_branch_idx',
                                     'original_line_idx',
                                     'original_tr_idx',
                                     'original_gen_idx',
                                     'original_bat_idx',
                                     'iPfsh',
                                     'iQfma',
                                     'iBeqz',
                                     'iBeqv',
                                     'iVtma',
                                     'iQtma',
                                     'iPfdp',
                                     'iVscL',
                                     'VfBeqbus',
                                     'Vtmabus'
                                     ]

    def get_injections(self, normalize=True):
        """
        Compute the power
        :return: return the array of power injections in MW if normalized is false, in p.u. otherwise
        """

        # load
        Sbus = self.load_data.get_injections_per_bus()  # MW (negative already)

        # static generators
        Sbus += self.static_generator_data.get_injections_per_bus()

        # generators
        Sbus += self.generator_data.get_injections_per_bus()

        # battery
        Sbus += self.battery_data.get_injections_per_bus()

        # HVDC forced power is not handled here because of the possible islands

        if normalize:
            Sbus /= self.Sbase

        return Sbus

    def consolidate_information(self, use_stored_guess=False):
        """
        Consolidates the information of this object
        :return:
        """

        self.nbus = len(self.bus_data)
        self.nline = len(self.line_data)
        self.ndcline = len(self.dc_line_data)
        self.ntr = len(self.transformer_data)
        self.nvsc = len(self.vsc_data)
        self.nhvdc = len(self.hvdc_data)
        self.nload = len(self.load_data)
        self.ngen = len(self.generator_data)
        self.nbatt = len(self.battery_data)
        self.nshunt = len(self.shunt_data)
        self.nstagen = len(self.static_generator_data)
        self.nupfc = len(self.upfc_data)
        self.nbr = self.nline + self.ntr + self.nvsc + self.ndcline + self.nupfc

        self.original_bus_idx = np.arange(self.nbus)
        self.original_branch_idx = np.arange(self.nbr)
        self.original_line_idx = np.arange(self.nline)
        self.original_tr_idx = np.arange(self.ntr)
        self.original_dc_line_idx = np.arange(self.ndcline)
        self.original_vsc_idx = np.arange(self.nvsc)
        self.original_hvdc_idx = np.arange(self.nhvdc)
        self.original_gen_idx = np.arange(self.ngen)
        self.original_bat_idx = np.arange(self.nbatt)
        self.original_load_idx = np.arange(self.nload)
        self.original_stagen_idx = np.arange(self.nstagen)
        self.original_shunt_idx = np.arange(self.nshunt)

        self.branch_data.C_branch_bus_f = self.branch_data.C_branch_bus_f.tocsc()
        self.branch_data.C_branch_bus_t = self.branch_data.C_branch_bus_t.tocsc()

        self.line_data.C_line_bus = self.line_data.C_line_bus.tocsc()
        self.dc_line_data.C_dc_line_bus = self.dc_line_data.C_dc_line_bus.tocsc()
        self.transformer_data.C_tr_bus = self.transformer_data.C_tr_bus.tocsc()
        self.hvdc_data.C_hvdc_bus_f = self.hvdc_data.C_hvdc_bus_f.tocsc()
        self.hvdc_data.C_hvdc_bus_t = self.hvdc_data.C_hvdc_bus_t.tocsc()
        self.vsc_data.C_vsc_bus = self.vsc_data.C_vsc_bus.tocsc()
        self.upfc_data.C_elm_bus = self.upfc_data.C_elm_bus.tocsc()

        self.load_data.C_bus_load = self.load_data.C_bus_load.tocsr()
        self.battery_data.C_bus_batt = self.battery_data.C_bus_batt.tocsr()
        self.generator_data.C_bus_gen = self.generator_data.C_bus_gen.tocsr()
        self.shunt_data.C_bus_shunt = self.shunt_data.C_bus_shunt.tocsr()
        self.static_generator_data.C_bus_static_generator = self.static_generator_data.C_bus_static_generator.tocsr()

        self.bus_data.bus_installed_power = self.generator_data.get_installed_power_per_bus()
        self.bus_data.bus_installed_power += self.battery_data.get_installed_power_per_bus()

        if not use_stored_guess:
            self.bus_data.Vbus = compose_generator_voltage_profile(nbus=self.nbus,
                                                                   ntime=self.ntime,
                                                                   gen_bus_indices=self.generator_data.get_bus_indices(),
                                                                   gen_vset=self.generator_data.generator_v,
                                                                   gen_status=self.generator_data.generator_active,
                                                                   gen_is_controlled=self.generator_data.generator_controllable,
                                                                   bat_bus_indices=self.battery_data.get_bus_indices(),
                                                                   bat_vset=self.battery_data.battery_v,
                                                                   bat_status=self.battery_data.battery_active,
                                                                   bat_is_controlled=self.battery_data.battery_controllable,
                                                                   hvdc_bus_f=self.hvdc_data.get_bus_indices_f(),
                                                                   hvdc_bus_t=self.hvdc_data.get_bus_indices_t(),
                                                                   hvdc_status=self.hvdc_data.active,
                                                                   hvdc_vf=self.hvdc_data.Vset_f,
                                                                   hvdc_vt=self.hvdc_data.Vset_t,
                                                                   iBeqv=np.array(self.iBeqv, dtype=int),
                                                                   iVtma=np.array(self.iVtma, dtype=int),
                                                                   VfBeqbus=np.array(self.VfBeqbus, dtype=int),
                                                                   Vtmabus=np.array(self.Vtmabus, dtype=int),
                                                                   branch_status=self.branch_data.branch_active,
                                                                   br_vf=self.branch_data.vf_set,
                                                                   br_vt=self.branch_data.vt_set)

        self.determine_control_indices()

    def re_calc_admittance_matrices(self, tap_module, t=0, idx=None):
        """
        Fast admittance recombination
        :param tap_module: transformer taps (if idx is provided, must have the same length as idx,
                           otherwise the length must be the number of branches)
        :param t: time index, 0 by default
        :param idx: Indices of the branches where the tap belongs,
                    if None assumes that the tap sizes is equal to the number of branches
        :return:
        """
        if idx is None:
            Ybus_, Yf_, Yt_ = self.Admittances.modify_taps(self.branch_data.m[:, t], tap_module)
        else:
            Ybus_, Yf_, Yt_ = self.Admittances.modify_taps(self.branch_data.m[np.ix_(idx, t)], tap_module)

        self.Admittances.Ybus = Ybus_
        self.Admittances.Yf = Yf_
        self.Admittances.Yt = Yt_

    def determine_control_indices(self):
        """
        This function fills in the lists of indices to control different magnitudes

        :returns idx_sh, idx_qz, idx_vf, idx_vt, idx_qt, VfBeqbus, Vtmabus

        VSC Control modes:

        in the paper's scheme:
        from -> DC
        to   -> AC

        |   Mode    |   const.1 |   const.2 |   type    |
        -------------------------------------------------
        |   1       |   theta   |   Vac     |   I       |
        |   2       |   Pf      |   Qac     |   I       |
        |   3       |   Pf      |   Vac     |   I       |
        -------------------------------------------------
        |   4       |   Vdc     |   Qac     |   II      |
        |   5       |   Vdc     |   Vac     |   II      |
        -------------------------------------------------
        |   6       | Vdc droop |   Qac     |   III     |
        |   7       | Vdc droop |   Vac     |   III     |
        -------------------------------------------------

        Indices where each control goes:
        mismatch  →  |  ∆Pf	Qf	Qf  Qt	∆Qt
        variable  →  |  Ɵsh	Beq	m	m	Beq
        Indices   →  |  Ish	Iqz	Ivf	Ivt	Iqt
        ------------------------------------
        VSC 1	     |  -	1	-	1	-   |   AC voltage control (voltage “to”)
        VSC 2	     |  1	1	-	-	1   |   Active and reactive power control
        VSC 3	     |  1	1	-	1	-   |   Active power and AC voltage control
        VSC 4	     |  -	-	1	-	1   |   Dc voltage and Reactive power flow control
        VSC 5	     |  -	-	-	1	1   |   Ac and Dc voltage control
        ------------------------------------
        Transformer 0|	-	-	-	-	-   |   Fixed transformer
        Transformer 1|	1	-	-	-	-   |   Phase shifter → controls power
        Transformer 2|	-	-	1	-	-   |   Control the voltage at the “from” side
        Transformer 3|	-	-	-	1	-   |   Control the voltage at the “to” side
        Transformer 4|	1	-	1	-	-   |   Control the power flow and the voltage at the “from” side
        Transformer 5|	1	-	-	1	-   |   Control the power flow and the voltage at the “to” side
        ------------------------------------

        """

        # indices in the global branch scheme
        self.iPfsh = list()  # indices of the branches controlling Pf flow with theta sh
        self.iQfma = list()  # indices of the branches controlling Qf with ma
        self.iBeqz = list()  # indices of the branches when forcing the Qf flow to zero (aka "the zero condition")
        self.iBeqv = list()  # indices of the branches when controlling Vf with Beq
        self.iVtma = list()  # indices of the branches when controlling Vt with ma
        self.iQtma = list()  # indices of the branches controlling the Qt flow with ma
        self.iPfdp = list()  # indices of the drop converters controlling the power flow with theta sh
        self.iVscL = list()  # indices of the converters
        self.iPfdp_va = list()

        self.any_control = False

        for k, tpe in enumerate(self.branch_data.control_mode):

            if tpe == TransformerControlType.fixed:
                pass

            elif tpe == TransformerControlType.Pt:
                self.iPfsh.append(k)
                self.any_control = True

            elif tpe == TransformerControlType.Qt:
                self.iQtma.append(k)
                self.any_control = True

            elif tpe == TransformerControlType.PtQt:
                self.iPfsh.append(k)
                self.iQtma.append(k)
                self.any_control = True

            elif tpe == TransformerControlType.Vt:
                self.iVtma.append(k)
                self.any_control = True

            elif tpe == TransformerControlType.PtVt:
                self.iPfsh.append(k)
                self.iVtma.append(k)
                self.any_control = True

            # VSC ------------------------------------------------------------------------------------------------------
            elif tpe == ConverterControlType.type_0_free:  # 1a:Free
                self.iBeqz.append(k)

                self.iVscL.append(k)
                self.any_control = True

            elif tpe == ConverterControlType.type_I_1:  # 1:Vac
                self.iVtma.append(k)
                self.iBeqz.append(k)

                self.iVscL.append(k)
                self.any_control = True

            elif tpe == ConverterControlType.type_I_2:  # 2:Pdc+Qac

                self.iPfsh.append(k)
                self.iQtma.append(k)
                self.iBeqz.append(k)

                self.iVscL.append(k)
                self.any_control = True

            elif tpe == ConverterControlType.type_I_3:  # 3:Pdc+Vac
                self.iPfsh.append(k)
                self.iVtma.append(k)
                self.iBeqz.append(k)

                self.iVscL.append(k)
                self.any_control = True

            elif tpe == ConverterControlType.type_II_4:  # 4:Vdc+Qac
                self.iBeqv.append(k)
                self.iQtma.append(k)

                self.iVscL.append(k)
                self.any_control = True

            elif tpe == ConverterControlType.type_II_5:  # 5:Vdc+Vac
                self.iBeqv.append(k)
                self.iVtma.append(k)

                self.iVscL.append(k)
                self.any_control = True

            elif tpe == ConverterControlType.type_III_6:  # 6:Droop+Qac
                self.iPfdp.append(k)
                self.iQtma.append(k)

                self.iVscL.append(k)
                self.any_control = True

            elif tpe == ConverterControlType.type_III_7:  # 4a:Droop-slack
                self.iPfdp.append(k)
                self.iVtma.append(k)

                self.iVscL.append(k)
                self.any_control = True

            elif tpe == ConverterControlType.type_IV_I:  # 8:Vdc
                self.iBeqv.append(k)
                self.iVscL.append(k)

                self.any_control = True

            elif tpe == ConverterControlType.type_IV_II:  # 9:Pdc
                self.iPfsh.append(k)
                self.iBeqz.append(k)

                self.any_control = True

            elif tpe == 0:
                pass  # required for the no-control case

            else:
                raise Exception('Unknown control type:' + str(tpe))

        # VfBeqbus_sh = list()
        # for k, is_controlled in enumerate(self.shunt_data.get_controlled_per_bus()):
        #     if is_controlled:
        #         VfBeqbus_sh.append(k)
        #         self.any_control = True

        # FUBM- Saves the "from" bus identifier for Vf controlled by Beq
        #  (Converters type II for Vdc control)
        self.VfBeqbus = self.F[self.iBeqv]

        # FUBM- Saves the "to"   bus identifier for Vt controlled by ma
        #  (Converters and Transformers)
        self.Vtmabus = self.T[self.iVtma]

        self.iPfsh = np.array(self.iPfsh, dtype=np.int)
        self.iQfma = np.array(self.iQfma, dtype=np.int)
        self.iBeqz = np.array(self.iBeqz, dtype=np.int)
        self.iBeqv = np.array(self.iBeqv, dtype=np.int)
        self.iVtma = np.array(self.iVtma, dtype=np.int)
        self.iQtma = np.array(self.iQtma, dtype=np.int)
        self.iPfdp = np.array(self.iPfdp, dtype=np.int)
        self.iPfdp_va = np.array(self.iPfdp_va, dtype=np.int)
        self.iVscL = np.array(self.iVscL, dtype=np.int)

    def get_branch_df(self, t=0):
        return self.branch_data.to_df(t)

    @property
    def line_idx(self):
        return slice(0, self.nline, 1)

    @property
    def transformer_idx(self):
        return slice(self.nline, self.nline + self.ntr, 1)

    @property
    def vsc_idx(self):
        return slice(self.nline + self.ntr, self.nline + self.ntr + self.nvsc, 1)

    @property
    def dc_line_idx(self):
        return slice(self.nline + self.ntr + self.nvsc, self.nline + self.ntr + self.nvsc + self.ndcline, 1)

    @property
    def Vbus(self):

        if self.Vbus_ is None:
            self.Vbus_ = self.bus_data.Vbus

        return self.Vbus_[:, 0]

    @property
    def Sbus(self):
        """
        Returns the power injections in per-unit
        :return: array of power injections (p.u.)
        """

        if self.Sbus_ is None:
            self.Sbus_ = self.get_injections(normalize=True)

        return self.Sbus_[:, 0]

    @property
    def Ibus(self):

        if self.Ibus_ is None:
            self.Ibus_ = self.load_data.get_current_injections_per_bus() / self.Sbase

        return self.Ibus_[:, 0]

    @property
    def YLoadBus(self):

        if self.YloadBus_ is None:
            self.YloadBus_ = self.load_data.get_admittance_injections_per_bus() / self.Sbase

        return self.YloadBus_[:, 0]

    @property
    def Rates(self):
        return self.branch_data.branch_rates[:, 0]

    @property
    def ContingencyRates(self):
        return self.branch_data.branch_contingency_rates[:, 0]

    @property
    def Qmax_bus(self):

        if self.Qmax_bus_ is None:
            self.Qmax_bus_, self.Qmin_bus_ = self.compute_reactive_power_limits()

        return self.Qmax_bus_

    @property
    def Qmin_bus(self):

        if self.Qmin_bus_ is None:
            self.Qmax_bus_, self.Qmin_bus_ = self.compute_reactive_power_limits()

        return self.Qmin_bus_

    @property
    def Bmax_bus(self):

        if self.Bmax_bus_ is None:
            self.Bmax_bus_, self.Bmin_bus_ = self.compute_susceptance_limits()

        return self.Bmax_bus_

    @property
    def Bmin_bus(self):

        if self.Bmin_bus_ is None:
            self.Bmax_bus_, self.Bmin_bus_ = self.compute_susceptance_limits()

        return self.Bmin_bus_

    @property
    def Yshunt_from_devices(self):

        # compute on demand and store
        if self.Yshunt_from_devices_ is None:
            self.Yshunt_from_devices_ = self.shunt_data.get_injections_per_bus() / self.Sbase

        return self.Yshunt_from_devices_

    @property
    def bus_types(self):
        return self.bus_data.bus_types

    @property
    def bus_installed_power(self):
        return self.bus_data.bus_installed_power

    @property
    def bus_names(self):
        return self.bus_data.bus_names

    @property
    def branch_names(self):
        return self.branch_data.branch_names

    @property
    def load_names(self):
        return self.load_data.load_names

    @property
    def generator_names(self):
        return self.generator_data.generator_names

    @property
    def battery_names(self):
        return self.battery_data.battery_names

    @property
    def tr_names(self):
        return self.transformer_data.tr_names

    @property
    def hvdc_names(self):
        return self.hvdc_data.names

    @property
    def tr_tap_position(self):
        return self.transformer_data.tr_tap_position

    @property
    def tr_tap_mod(self):
        return self.transformer_data.tr_tap_mod

    @property
    def tr_bus_to_regulated_idx(self):
        return self.transformer_data.tr_bus_to_regulated_idx

    @property
    def tr_max_tap(self):
        return self.transformer_data.tr_max_tap

    @property
    def tr_min_tap(self):
        return self.transformer_data.tr_min_tap

    @property
    def tr_tap_inc_reg_up(self):
        return self.transformer_data.tr_tap_inc_reg_up

    @property
    def tr_tap_inc_reg_down(self):
        return self.transformer_data.tr_tap_inc_reg_down

    @property
    def tr_vset(self):
        return self.transformer_data.tr_vset

    @property
    def F(self):
        return self.branch_data.F

    @property
    def T(self):
        return self.branch_data.T

    @property
    def branch_rates(self):
        return self.branch_data.branch_rates[:, 0]


    @property
    def ac_indices(self):
        """
        Array of indices of the AC branches
        :return: array of indices
        """
        if self.ac_ is None:
            self.ac_ = self.branch_data.get_ac_indices()

        return self.ac_

    @property
    def dc_indices(self):
        """
        Array of indices of the DC branches
        :return: array of indices
        """
        if self.dc_ is None:
            self.dc_ = self.branch_data.get_dc_indices()

        return self.dc_

    @property
    def Cf(self):
        """
        Connectivity matrix of the "from" nodes
        :return: CSC matrix
        """
        # compute on demand and store
        if self.Cf_ is None:
            self.Cf_, self.Ct_ = ycalc.compute_connectivity(branch_active=self.branch_data.branch_active[:, 0],
                                                            Cf_=self.branch_data.C_branch_bus_f,
                                                            Ct_=self.branch_data.C_branch_bus_t)

        if not isinstance(self.Cf_, sp.csc_matrix):
            self.Cf_ = self.Cf_.tocsc()

        return self.Cf_

    @property
    def Ct(self):
        """
        Connectivity matrix of the "to" nodes
        :return: CSC matrix
        """
        # compute on demand and store
        if self.Ct_ is None:
            self.Cf_, self.Ct_ = ycalc.compute_connectivity(branch_active=self.branch_data.branch_active[:, 0],
                                                            Cf_=self.branch_data.C_branch_bus_f,
                                                            Ct_=self.branch_data.C_branch_bus_t)

        if not isinstance(self.Ct_, sp.csc_matrix):
            self.Ct_ = self.Ct_.tocsc()

        return self.Ct_

    @property
    def A(self):
        """
        Connectivity matrix
        :return: CSC matrix
        """

        if self.A_ is None:
            self.A_ = (self.Cf - self.Ct).tocsc()

        return self.A_

    @property
    def Ybus(self):
        """
        Admittance matrix
        :return: CSC matrix
        """

        # compute admittances on demand
        if self.Admittances is None:

            self.Admittances = ycalc.compute_admittances(R=self.branch_data.R,
                                                         X=self.branch_data.X,
                                                         G=self.branch_data.G,
                                                         B=self.branch_data.B,
                                                         k=self.branch_data.k,
                                                         tap_module=self.branch_data.m[:, 0],
                                                         vtap_f=self.branch_data.tap_f,
                                                         vtap_t=self.branch_data.tap_t,
                                                         tap_angle=self.branch_data.theta[:, 0],
                                                         Beq=self.branch_data.Beq[:, 0],
                                                         Cf=self.Cf,
                                                         Ct=self.Ct,
                                                         G0sw=self.branch_data.G0sw[:, 0],
                                                         If=np.zeros(len(self.branch_data)),
                                                         a=self.branch_data.a,
                                                         b=self.branch_data.b,
                                                         c=self.branch_data.c,
                                                         Yshunt_bus=self.Yshunt_from_devices[:, 0])
        return self.Admittances.Ybus

    @property
    def Yf(self):
        """
        Admittance matrix of the "from" nodes with the branches
        :return: CSC matrix
        """
        if self.Admittances is None:
            x = self.Ybus  # call the constructor of Yf

        return self.Admittances.Yf

    @property
    def Yt(self):
        """
        Admittance matrix of the "to" nodes with the branches
        :return: CSC matrix
        """
        if self.Admittances is None:
            x = self.Ybus  # call the constructor of Yt

        return self.Admittances.Yt

    @property
    def Yseries(self):
        """
        Admittance matrix of the series elements of the pi model of the branches
        :return: CSC matrix
        """
        # compute admittances on demand
        if self.Yseries_ is None:

            self.Yseries_, self.Yshunt_ = ycalc.compute_split_admittances(R=self.branch_data.R,
                                                                          X=self.branch_data.X,
                                                                          G=self.branch_data.G,
                                                                          B=self.branch_data.B,
                                                                          k=self.branch_data.k,
                                                                          m=self.branch_data.m[:, 0],
                                                                          mf=self.branch_data.tap_f,
                                                                          mt=self.branch_data.tap_t,
                                                                          theta=self.branch_data.theta[:, 0],
                                                                          Beq=self.branch_data.Beq[:, 0],
                                                                          Cf=self.Cf,
                                                                          Ct=self.Ct,
                                                                          G0=self.branch_data.G0sw[:, 0],
                                                                          If=np.zeros(len(self.branch_data)),
                                                                          a=self.branch_data.a,
                                                                          b=self.branch_data.b,
                                                                          c=self.branch_data.c,
                                                                          Yshunt_bus=self.Yshunt_from_devices[:, 0])
        return self.Yseries_

    @property
    def Yshunt(self):
        """
        Array of shunt admittances of the pi model of the branches (used in HELM mostly)
        :return: Array of complex values
        """
        if self.Yshunt_ is None:
            x = self.Yseries  # call the constructor of Yshunt

        return self.Yshunt_

    # @property
    # def YshuntHelm(self):
    #     return self.Yshunt_from_devices[:, 0]

    @property
    def B1(self):
        """
        B' matrix of the fast decoupled method
        :return:
        """
        if self.B1_ is None:

            self.B1_, self.B2_ = ycalc.compute_fast_decoupled_admittances(X=self.branch_data.X,
                                                                          B=self.branch_data.B,
                                                                          m=self.branch_data.m[:, 0],
                                                                          mf=self.branch_data.vf_set[:, 0],
                                                                          mt=self.branch_data.vt_set[:, 0],
                                                                          Cf=self.Cf,
                                                                          Ct=self.Ct)
        return self.B1_

    @property
    def B2(self):
        """
        B'' matrix of the fast decoupled method
        :return:
        """
        if self.B2_ is None:
            x = self.B1  # call the constructor of B2

        return self.B2_

    @property
    def Bbus(self):
        """
        Susceptance matrix for the linear methods
        :return:
        """
        if self.Bbus_ is None:
            self.Bbus_, self.Bf_, self.Btheta_ = ycalc.compute_linear_admittances(nbr=self.nbr,
                                                                                  X=self.branch_data.X,
                                                                                  R=self.branch_data.R,
                                                                                  m=self.branch_data.m[:, 0],
                                                                                  active=self.branch_data.branch_active[:, 0],
                                                                                  Cf=self.Cf,
                                                                                  Ct=self.Ct,
                                                                                  ac=self.ac_indices,
                                                                                  dc=self.dc_indices)
            self.Bpqpv_ = self.Bbus_[np.ix_(self.pqpv, self.pqpv)].tocsc()
            self.Bref_ = self.Bbus_[np.ix_(self.pqpv, self.vd)].tocsc()

        return self.Bbus_

    @property
    def Bf(self):
        """
        Susceptance matrix of the "from" nodes to the branches
        :return:
        """
        if self.Bf_ is None:
            x = self.Bbus  # call the constructor of Bf

        return self.Bf_

    @property
    def Btheta(self):

        if self.Bf_ is None:
            x = self.Bbus  # call the constructor of Bf

        return self.Btheta_

    @property
    def Bpqpv(self):

        if self.Bpqpv_ is None:
            x = self.Bbus  # call the constructor of Bpqpv

        return self.Bpqpv_

    @property
    def Bref(self):

        if self.Bref_ is None:
            x = self.Bbus  # call the constructor of Bref

        return self.Bref_

    @property
    def vd(self):

        if self.vd_ is None:
            self.vd_, self.pq_, self.pv_, self.pqpv_ = compile_types(Sbus=self.Sbus,
                                                                     types=self.bus_data.bus_types)

        return self.vd_

    @property
    def pq(self):

        if self.pq_ is None:
            x = self.vd  # call the constructor

        return self.pq_

    @property
    def pv(self):

        if self.pv_ is None:
            x = self.vd  # call the constructor

        return self.pv_

    @property
    def pqpv(self):

        if self.pqpv_ is None:
            x = self.vd  # call the constructor

        return self.pqpv_

    def compute_reactive_power_limits(self):
        """
        compute the reactive power limits in place
        :return: Qmax_bus, Qmin_bus in per unit
        """
        # generators
        Qmax_bus = self.generator_data.get_qmax_per_bus()
        Qmin_bus = self.generator_data.get_qmin_per_bus()

        if self.nbatt > 0:
            # batteries
            Qmax_bus += self.battery_data.get_qmax_per_bus()
            Qmin_bus += self.battery_data.get_qmin_per_bus()

        if self.nshunt > 0:
            # shunts
            Qmax_bus += self.shunt_data.get_b_max_per_bus()
            Qmin_bus += self.shunt_data.get_b_min_per_bus()

        if self.nhvdc > 0:
            # hvdc from
            Qmax_bus += self.hvdc_data.get_qmax_from_per_bus()
            Qmin_bus += self.hvdc_data.get_qmin_from_per_bus()

            # hvdc to
            Qmax_bus += self.hvdc_data.get_qmax_to_per_bus()
            Qmin_bus += self.hvdc_data.get_qmin_to_per_bus()

        # fix zero values
        Qmax_bus[Qmax_bus == 0] = 1e20
        Qmin_bus[Qmin_bus == 0] = -1e20

        return Qmax_bus / self.Sbase, Qmin_bus / self.Sbase

    def compute_susceptance_limits(self):

        Bmin = self.shunt_data.get_b_min_per_bus() / self.Sbase
        Bmax = self.shunt_data.get_b_max_per_bus() / self.Sbase

        return Bmax, Bmin

    def get_inter_areas_branches(self, buses_areas_1, buses_areas_2):
        """
        Get the branches that join two areas
        :param buses_areas_1: Area from
        :param buses_areas_2: Area to
        :return: List of (branch index, branch object, flow sense w.r.t the area exchange)
        """
        return get_inter_areas_branch(self.branch_data.F, self.branch_data.T, buses_areas_1, buses_areas_2)

    def get_inter_areas_hvdc(self, buses_areas_1, buses_areas_2):
        """
        Get the branches that join two areas
        :param buses_areas_1: Area from
        :param buses_areas_2: Area to
        :return: List of (branch index, branch object, flow sense w.r.t the area exchange)
        """
        F = self.hvdc_data.get_bus_indices_f()
        T = self.hvdc_data.get_bus_indices_t()
        return get_inter_areas_branch(F, T, buses_areas_1, buses_areas_2)

    def get_generators_per_areas(self, buses_in_a1, buses_in_a2):
        """
        Get the generators that belong to the Area 1, Area 2 and the rest of areas
        :param buses_in_a1: List of bus indices of the area 1
        :param buses_in_a2: List of bus indices of the area 2
        :return: Tree lists: (gens_in_a1, gens_in_a2, gens_out)
                 each of the lists contains (bus index, generator index) tuples
        """
        if isinstance(self.generator_data.C_bus_gen, sp.csc_matrix):
            Cgen = self.generator_data.C_bus_gen
        else:
            Cgen = self.generator_data.C_bus_gen.tocsc()

        return get_devices_per_areas(Cgen, buses_in_a1, buses_in_a2)

    def get_batteries_per_areas(self, buses_in_a1, buses_in_a2):
        """
        Get the batteries that belong to the Area 1, Area 2 and the rest of areas
        :param buses_in_a1: List of bus indices of the area 1
        :param buses_in_a2: List of bus indices of the area 2
        :return: Tree lists: (batteries_in_a1, batteries_in_a2, batteries_out)
                 each of the lists contains (bus index, generator index) tuples
        """
        if isinstance(self.battery_data.C_bus_batt, sp.csc_matrix):
            Cgen = self.battery_data.C_bus_batt
        else:
            Cgen = self.battery_data.C_bus_batt.tocsc()

        return get_devices_per_areas(Cgen, buses_in_a1, buses_in_a2)

    def get_structure(self, structure_type) -> pd.DataFrame:
        """
        Get a DataFrame with the input.

        Arguments:

            **structure_type** (str): 'Vbus', 'Sbus', 'Ibus', 'Ybus', 'Yshunt', 'Yseries' or 'Types'

        Returns:

            pandas DataFrame

        """

        if structure_type == 'Vbus':

            df = pd.DataFrame(data=self.Vbus, columns=['Voltage (p.u.)'], index=self.bus_data.bus_names)

        elif structure_type == 'Sbus':
            df = pd.DataFrame(data=self.Sbus, columns=['Power (p.u.)'], index=self.bus_data.bus_names)

        elif structure_type == 'Ibus':
            df = pd.DataFrame(data=self.Ibus, columns=['Current (p.u.)'], index=self.bus_data.bus_names)

        elif structure_type == 'tap_f':
            df = pd.DataFrame(data=self.branch_data.tap_f,
                              columns=['Virtual tap from (p.u.)'],
                              index=self.branch_data.branch_names)

        elif structure_type == 'tap_t':
            df = pd.DataFrame(data=self.branch_data.tap_t,
                              columns=['Virtual tap to (p.u.)'],
                              index=self.branch_data.branch_names)

        elif structure_type == 'Ybus':
            df = pd.DataFrame(data=self.Ybus.toarray(),
                              columns=self.bus_data.bus_names,
                              index=self.bus_data.bus_names)

        elif structure_type == 'G':
            df = pd.DataFrame(data=self.Ybus.real.toarray(),
                              columns=self.bus_data.bus_names,
                              index=self.bus_data.bus_names)

        elif structure_type == 'B':
            df = pd.DataFrame(data=self.Ybus.imag.toarray(),
                              columns=self.bus_data.bus_names,
                              index=self.bus_data.bus_names)

        elif structure_type == 'Yf':
            df = pd.DataFrame(data=self.Yf.toarray(),
                              columns=self.bus_data.bus_names,
                              index=self.branch_data.branch_names)

        elif structure_type == 'Yt':
            df = pd.DataFrame(data=self.Yt.toarray(),
                              columns=self.bus_data.bus_names,
                              index=self.branch_data.branch_names)

        elif structure_type == 'Bbus':
            df = pd.DataFrame(data=self.Bbus.toarray(),
                              columns=self.bus_data.bus_names,
                              index=self.bus_data.bus_names)

        elif structure_type == 'Bf':
            df = pd.DataFrame(data=self.Bf.toarray(),
                              columns=self.bus_data.bus_names,
                              index=self.branch_data.branch_names)

        elif structure_type == 'Cf':
            df = pd.DataFrame(data=self.Cf.toarray(),
                              columns=self.bus_data.bus_names,
                              index=self.branch_data.branch_names)

        elif structure_type == 'Ct':
            df = pd.DataFrame(data=self.Ct.toarray(),
                              columns=self.bus_data.bus_names,
                              index=self.branch_data.branch_names)

        elif structure_type == 'Yshunt':
            df = pd.DataFrame(data=self.Yshunt, columns=['Shunt admittance (p.u.)'], index=self.bus_data.bus_names)

        elif structure_type == 'Yseries':
            df = pd.DataFrame(data=self.Yseries.toarray(),
                              columns=self.bus_data.bus_names,
                              index=self.bus_data.bus_names)

        elif structure_type == "B'":
            df = pd.DataFrame(data=self.B1.toarray(), columns=self.bus_data.bus_names, index=self.bus_data.bus_names)

        elif structure_type == "B''":
            df = pd.DataFrame(data=self.B2.toarray(), columns=self.bus_data.bus_names, index=self.bus_data.bus_names)

        elif structure_type == 'Types':
            df = pd.DataFrame(data=self.bus_types, columns=['Bus types'], index=self.bus_data.bus_names)

        elif structure_type == 'Qmin':
            df = pd.DataFrame(data=self.Qmin_bus, columns=['Qmin'], index=self.bus_data.bus_names)

        elif structure_type == 'Qmax':
            df = pd.DataFrame(data=self.Qmax_bus, columns=['Qmax'], index=self.bus_data.bus_names)

        elif structure_type == 'pq':
            df = pd.DataFrame(data=self.pq, columns=['pq'], index=self.bus_data.bus_names[self.pq])

        elif structure_type == 'pv':
            df = pd.DataFrame(data=self.pv, columns=['pv'], index=self.bus_data.bus_names[self.pv])

        elif structure_type == 'vd':
            df = pd.DataFrame(data=self.vd, columns=['vd'], index=self.bus_data.bus_names[self.vd])

        elif structure_type == 'pqpv':
            df = pd.DataFrame(data=self.pqpv, columns=['pqpv'], index=self.bus_data.bus_names[self.pqpv])

        elif structure_type == 'original_bus_idx':
            df = pd.DataFrame(data=self.original_bus_idx, columns=['original_bus_idx'], index=self.bus_data.bus_names)

        elif structure_type == 'original_branch_idx':
            df = pd.DataFrame(data=self.original_branch_idx,
                              columns=['original_branch_idx'],
                              index=self.branch_data.branch_names)

        elif structure_type == 'original_line_idx':
            df = pd.DataFrame(data=self.original_line_idx,
                              columns=['original_line_idx'],
                              index=self.line_data.line_names)

        elif structure_type == 'original_tr_idx':
            df = pd.DataFrame(data=self.original_tr_idx,
                              columns=['original_tr_idx'],
                              index=self.transformer_data.tr_names)

        elif structure_type == 'original_gen_idx':
            df = pd.DataFrame(data=self.original_gen_idx,
                              columns=['original_gen_idx'],
                              index=self.generator_data.generator_names)

        elif structure_type == 'original_bat_idx':
            df = pd.DataFrame(data=self.original_bat_idx,
                              columns=['original_bat_idx'],
                              index=self.battery_data.battery_names)

        elif structure_type == 'Jacobian':

            pvpq = np.r_[self.pv, self.pq]
            i2 = np.r_[self.pq, self.VfBeqbus, self.Vtmabus]
            i4 = np.r_[self.iQfma, self.iBeqz]

            cols = ['1) dVa {0}'.format(i) for i in pvpq]
            cols += ['2) dVm {0}'.format(i) for i in self.pq]
            cols += ['3) dPfsh {0}'.format(i) for i in self.iPfsh]
            cols += ['4) dQfma {0}'.format(i) for i in self.iQfma]
            cols += ['5) dBeqz {0}'.format(i) for i in self.iBeqz]
            cols += ['6) dBeqv {0}'.format(i) for i in self.iBeqv]
            cols += ['7) dVtma {0}'.format(i) for i in self.iVtma]
            cols += ['8) dQtma {0}'.format(i) for i in self.iQtma]
            cols += ['9) dPfdp {0}'.format(i) for i in self.iPfdp]

            rows = ['1) dP {0}'.format(i) for i in pvpq]
            rows += ['2) dQ {0}'.format(i) for i in i2]
            rows += ['3) dQ {0}'.format(i) for i in self.iBeqv]
            rows += ['4) dQ {0}'.format(i) for i in self.iVtma]
            rows += ['5) dPf {0}'.format(i) for i in self.iPfsh]
            rows += ['6) dQf {0}'.format(i) for i in self.iQfma]
            rows += ['7) dQf {0}'.format(i) for i in self.iBeqz]
            rows += ['8) dQt {0}'.format(i) for i in self.iQtma]
            rows += ['9) dPfdp {0}'.format(i) for i in self.iPfdp]

            # compute admittances
            Ys = 1.0 / (self.branch_data.R + 1j * self.branch_data.X)
            Ybus, Yf, Yt, tap = ycalc.compile_y_acdc(Cf=self.Cf, Ct=self.Ct,
                                                     C_bus_shunt=self.shunt_data.C_bus_shunt,
                                                     shunt_admittance=self.shunt_data.shunt_admittance[:, 0],
                                                     shunt_active=self.shunt_data.shunt_active[:, 0],
                                                     ys=Ys,
                                                     B=self.branch_data.B,
                                                     Sbase=self.Sbase,
                                                     m=self.branch_data.m[:, 0],
                                                     theta=self.branch_data.theta[:, 0],
                                                     Beq=self.branch_data.Beq[:, 0],
                                                     Gsw=self.branch_data.G0sw[:, 0],
                                                     mf=self.branch_data.tap_f,
                                                     mt=self.branch_data.tap_t)

            J = fubm_jacobian(self.nbus, self.nbr, self.iPfsh, self.iPfdp,
                              self.iQfma, self.iQtma, self.iVtma, self.iBeqz, self.iBeqv,
                              self.VfBeqbus, self.Vtmabus,
                              self.F, self.T, Ys,
                              self.branch_data.k, tap, self.branch_data.m[:, 0], self.branch_data.B,
                              self.branch_data.Beq[:, 0], self.branch_data.Kdp,
                              self.Vbus, Ybus.tocsc(), Yf.tocsc(), Yt.tocsc(), self.Cf.tocsc(), self.Ct.tocsc(),
                              pvpq, self.pq)

            df = pd.DataFrame(data=J.toarray(), columns=cols, index=rows)


        elif structure_type == 'iPfsh':
            df = pd.DataFrame(data=self.iPfsh, columns=['iPfsh'], index=self.branch_data.branch_names[self.iPfsh])

        elif structure_type == 'iQfma':
            df = pd.DataFrame(data=self.iQfma, columns=['iQfma'], index=self.branch_data.branch_names[self.iQfma])

        elif structure_type == 'iBeqz':
            df = pd.DataFrame(data=self.iBeqz, columns=['iBeqz'], index=self.branch_data.branch_names[self.iBeqz])

        elif structure_type == 'iBeqv':
            df = pd.DataFrame(data=self.iBeqv, columns=['iBeqv'], index=self.branch_data.branch_names[self.iBeqv])

        elif structure_type == 'iVtma':
            df = pd.DataFrame(data=self.iVtma, columns=['iVtma'], index=self.branch_data.branch_names[self.iVtma])

        elif structure_type == 'iQtma':
            df = pd.DataFrame(data=self.iQtma, columns=['iQtma'], index=self.branch_data.branch_names[self.iQtma])

        elif structure_type == 'iPfdp':
            df = pd.DataFrame(data=self.iPfdp, columns=['iPfdp'], index=self.branch_data.branch_names[self.iPfdp])

        elif structure_type == 'iVscL':
            df = pd.DataFrame(data=self.iVscL, columns=['iVscL'], index=self.branch_data.branch_names[self.iVscL])

        elif structure_type == 'VfBeqbus':
            df = pd.DataFrame(data=self.VfBeqbus, columns=['VfBeqbus'], index=self.bus_data.bus_names[self.VfBeqbus])

        elif structure_type == 'Vtmabus':
            df = pd.DataFrame(data=self.Vtmabus, columns=['Vtmabus'], index=self.bus_data.bus_names[self.Vtmabus])
        else:

            raise Exception('PF input: structure type not found' +  str(structure_type))

        return df

    def get_island(self, bus_idx, time_idx=None) -> "SnapshotData":
        """
        Get the island corresponding to the given buses
        :param bus_idx: array of bus indices
        :param time_idx: array of time indices (or None for all time indices)
        :return: SnapshotData
        """

        # if the island is the same as the original bus indices, no slicing is needed
        if len(bus_idx) == len(self.original_bus_idx):
            if np.all(bus_idx == self.original_bus_idx):
                return self

        # find the indices of the devices of the island
        line_idx = self.line_data.get_island(bus_idx)
        dc_line_idx = self.dc_line_data.get_island(bus_idx)
        tr_idx = self.transformer_data.get_island(bus_idx)
        vsc_idx = self.vsc_data.get_island(bus_idx)
        upfc_idx = self.upfc_data.get_island(bus_idx)
        hvdc_idx = self.hvdc_data.get_island(bus_idx)
        br_idx = self.branch_data.get_island(bus_idx)

        load_idx = self.load_data.get_island(bus_idx)
        stagen_idx = self.static_generator_data.get_island(bus_idx)
        gen_idx = self.generator_data.get_island(bus_idx)
        batt_idx = self.battery_data.get_island(bus_idx)
        shunt_idx = self.shunt_data.get_island(bus_idx)

        nc = SnapshotData(nbus=len(bus_idx),
                          nline=len(line_idx),
                          ndcline=len(dc_line_idx),
                          ntr=len(tr_idx),
                          nvsc=len(vsc_idx),
                          nupfc=len(upfc_idx),
                          nhvdc=len(hvdc_idx),
                          nload=len(load_idx),
                          ngen=len(gen_idx),
                          nbatt=len(batt_idx),
                          nshunt=len(shunt_idx),
                          nstagen=len(stagen_idx),
                          sbase=self.Sbase)

        # set the original indices
        nc.original_bus_idx = bus_idx
        nc.original_branch_idx = br_idx
        nc.original_line_idx = line_idx
        nc.original_tr_idx = tr_idx
        nc.original_dc_line_idx = dc_line_idx
        nc.original_vsc_idx = vsc_idx
        nc.original_upfc_idx = upfc_idx
        nc.original_hvdc_idx = hvdc_idx
        nc.original_gen_idx = gen_idx
        nc.original_bat_idx = batt_idx
        nc.original_load_idx = load_idx
        nc.original_stagen_idx = stagen_idx
        nc.original_shunt_idx = shunt_idx

        # slice data
        nc.bus_data = self.bus_data.slice(bus_idx, time_idx)
        nc.branch_data = self.branch_data.slice(br_idx, bus_idx, time_idx)
        nc.line_data = self.line_data.slice(line_idx, bus_idx, time_idx)
        nc.transformer_data = self.transformer_data.slice(tr_idx, bus_idx, time_idx)
        nc.hvdc_data = self.hvdc_data.slice(hvdc_idx, bus_idx, time_idx)
        nc.vsc_data = self.vsc_data.slice(vsc_idx, bus_idx, time_idx)
        nc.dc_line_data = self.dc_line_data.slice(dc_line_idx, bus_idx, time_idx)
        nc.load_data = self.load_data.slice(load_idx, bus_idx, time_idx)
        nc.static_generator_data = self.static_generator_data.slice(stagen_idx, bus_idx, time_idx)
        nc.battery_data = self.battery_data.slice(batt_idx, bus_idx, time_idx)
        nc.generator_data = self.generator_data.slice(gen_idx, bus_idx, time_idx)
        nc.shunt_data = self.shunt_data.slice(shunt_idx, bus_idx, time_idx)

        nc.determine_control_indices()

        return nc

    def split_into_islands(self, ignore_single_node_islands=False) -> List["SnapshotData"]:
        """
        Split circuit into islands
        :param ignore_single_node_islands: ignore islands composed of only one bus
        :return: List[NumericCircuit]
        """

        # compute the adjacency matrix
        A = tp.get_adjacency_matrix(C_branch_bus_f=self.Cf,
                                    C_branch_bus_t=self.Ct,
                                    branch_active=self.branch_data.branch_active[:, 0],
                                    bus_active=self.bus_data.bus_active[:, 0])

        # find the matching islands
        idx_islands = tp.find_islands(A, active=self.bus_data.bus_active[:, 0])

        circuit_islands = list()  # type: List[SnapshotData]

        for bus_idx in idx_islands:

            if ignore_single_node_islands:

                if len(bus_idx) > 1:
                    island = self.get_island(bus_idx)
                    circuit_islands.append(island)

            else:
                island = self.get_island(bus_idx)
                circuit_islands.append(island)

        return circuit_islands


def compile_snapshot_circuit(circuit: MultiCircuit, apply_temperature=False,
                             branch_tolerance_mode=BranchImpedanceMode.Specified,
                             opf_results=None,
                             use_stored_guess=False) -> SnapshotData:
    """

    :param circuit:
    :param apply_temperature:
    :param branch_tolerance_mode:
    :param opf_results:
    :param use_stored_guess:
    :return:
    """

    logger = Logger()

    # declare the numerical circuit
    nc = SnapshotData(nbus=0,
                      nline=0,
                      ndcline=0,
                      ntr=0,
                      nvsc=0,
                      nupfc=0,
                      nhvdc=0,
                      nload=0,
                      ngen=0,
                      nbatt=0,
                      nshunt=0,
                      nstagen=0,
                      sbase=circuit.Sbase)

    bus_dict = {bus: i for i, bus in enumerate(circuit.buses)}

    nc.bus_data = gc_compiler.get_bus_data(circuit=circuit, use_stored_guess=use_stored_guess)

    nc.load_data = gc_compiler.get_load_data(circuit=circuit,
                                             bus_dict=bus_dict,
                                             opf_results=opf_results)

    nc.static_generator_data = gc_compiler.get_static_generator_data(circuit=circuit,
                                                                     bus_dict=bus_dict)

    nc.generator_data = gc_compiler.get_generator_data(circuit=circuit,
                                                       bus_dict=bus_dict,
                                                       Vbus=nc.bus_data.Vbus,
                                                       logger=logger,
                                                       opf_results=opf_results,
                                                       use_stored_guess=use_stored_guess)

    nc.battery_data = gc_compiler.get_battery_data(circuit=circuit,
                                                   bus_dict=bus_dict,
                                                   Vbus=nc.bus_data.Vbus,
                                                   logger=logger,
                                                   opf_results=opf_results,
                                                   use_stored_guess=use_stored_guess)

    nc.shunt_data = gc_compiler.get_shunt_data(circuit=circuit,
                                               bus_dict=bus_dict,
                                               Vbus=nc.bus_data.Vbus,
                                               logger=logger,
                                               use_stored_guess=use_stored_guess)

    nc.line_data = gc_compiler.get_line_data(circuit=circuit,
                                             bus_dict=bus_dict,
                                             apply_temperature=apply_temperature,
                                             branch_tolerance_mode=branch_tolerance_mode)

    nc.transformer_data = gc_compiler.get_transformer_data(circuit=circuit,
                                                           bus_dict=bus_dict)

    nc.vsc_data = gc_compiler.get_vsc_data(circuit=circuit,
                                           bus_dict=bus_dict)

    nc.upfc_data = gc_compiler.get_upfc_data(circuit=circuit,
                                             bus_dict=bus_dict)

    nc.dc_line_data = gc_compiler.get_dc_line_data(circuit=circuit,
                                                   bus_dict=bus_dict,
                                                   apply_temperature=apply_temperature,
                                                   branch_tolerance_mode=branch_tolerance_mode)

    nc.branch_data = gc_compiler.get_branch_data(circuit=circuit,
                                                 bus_dict=bus_dict,
                                                 Vbus=nc.bus_data.Vbus,
                                                 apply_temperature=apply_temperature,
                                                 branch_tolerance_mode=branch_tolerance_mode,
                                                 opf_results=opf_results,
                                                 use_stored_guess=use_stored_guess)

    nc.hvdc_data = gc_compiler.get_hvdc_data(circuit=circuit,
                                             bus_dict=bus_dict,
                                             bus_types=nc.bus_data.bus_types,
                                             opf_results=opf_results)

    nc.consolidate_information(use_stored_guess=use_stored_guess)

    return nc
