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

import numpy as np
import numba as nb
from scipy.sparse import lil_matrix, diags, csc_matrix


import GridCalEngine.Simulations.PowerFlow.NumericalMethods.derivatives as deriv


class AcDcSolSlicer:

    def __init__(self, npq, npv, nVfBeqbus, nVtmabus, nPfsh, nQfma, nBeqz, nQtma, nPfdp):
        """
        Declare the slicing limits in the same order as the Jacobian rows
        :param npq:
        :param npv:
        :param nVfBeqbus:
        :param nVtmabus:
        :param nPfsh:
        :param nQfma:
        :param nBeqz:
        :param nQtma:
        :param nPfdp:
        """
        self.a0 = 0
        self.a1 = self.a0 + npq + npv
        self.a2 = self.a1 + npq
        self.a3 = self.a2 + nBeqz
        self.a4 = self.a3 + nVfBeqbus
        self.a5 = self.a4 + nQfma
        self.a6 = self.a5 + nQtma
        self.a7 = self.a6 + nVtmabus
        self.a8 = self.a7 + nPfsh
        self.a9 = self.a8 + nPfdp

    def split(self, dx):
        """
        Split the linear system solution
        :param dx:
        :return:
        """
        dVa = dx[self.a0:self.a1]
        dVm = dx[self.a1:self.a2]
        dBeq_zero = dx[self.a2:self.a3]
        dBeq_vf = dx[self.a3:self.a4]
        dma_Qf = dx[self.a4:self.a5]
        dma_Qt = dx[self.a5:self.a6]
        dma_Vt = dx[self.a6:self.a7]
        dtheta_Pf = dx[self.a7:self.a8]
        dtheta_Pd = dx[self.a8:self.a9]

        return dVa, dVm, dBeq_vf, dma_Vt, dtheta_Pf, dma_Qf, dBeq_zero, dma_Qt, dtheta_Pd


@nb.njit(cache=True)
def make_lookup(n, arr):
    lookup = np.zeros(n, dtype=np.int32)
    lookup[arr] = np.arange(len(arr), dtype=np.int32)
    return lookup


@nb.njit(cache=True)
def fill_acdc_jacobian_data(Jx, Ji, Jp, Yp, Yi, Ys,
                            dSbus_dVa_x, dSbus_dVm_x,
                            dSf_dVa_x, dSf_dVa_i, dSf_dVa_p,
                            dSf_dVm_x, dSf_dVm_i, dSf_dVm_p,
                            dSt_dVa_x, dSt_dVa_i, dSt_dVa_p,
                            dSt_dVm_x, dSt_dVm_i, dSt_dVm_p,
                            dPfdp_dVm_x, dPfdp_dVm_i, dPfdp_dVm_p,
                            pvpq, pvpq_lookup, npvpq, pq,
                            i2, ni2, i2_lookup, i4, ni4, i4_lookup,
                            i_pf_tau, n_pf_tau, i_pf_tau_lookup,
                            i_qt_m, n_qt_m, i_qt_m_lookup,
                            i_qf_m, n_qf_m,
                            i_vt_m, n_vt_m,
                            i_pf_dp, n_pf_dp, i_pf_dp_lookup,
                            i_zero_beq, n_zero_beq,
                            i_vf_beq, n_vf_beq,
                            F, T, V, tap_modules_m, k2, tap_complex, Bc, b_eq):
    """

    :param Jx:
    :param Ji:
    :param Jp:
    :param Yp:
    :param Yi:
    :param Ys: Array of branch series admittances
    :param dSbus_dVa_x:
    :param dSbus_dVm_x:
    :param dSf_dVa_x:
    :param dSf_dVa_i:
    :param dSf_dVa_p:
    :param dSf_dVm_x:
    :param dSf_dVm_i:
    :param dSf_dVm_p:
    :param dSt_dVa_x:
    :param dSt_dVa_i:
    :param dSt_dVa_p:
    :param dSt_dVm_x:
    :param dSt_dVm_i:
    :param dSt_dVm_p:
    :param dPfdp_dVm_x:
    :param dPfdp_dVm_i:
    :param dPfdp_dVm_p:
    :param pvpq: Array of pv and then pq bus indices (not sorted)
    :param pvpq_lookup:
    :param npvpq:
    :param pvpq: Array of pq bus indices (sorted)
    :param i2:
    :param ni2:
    :param i2_lookup:
    :param i4:
    :param ni4:
    :param i4_lookup:
    :param i_pf_tau: indices of the Pf controlled Branches with tau
    :param n_pf_tau:
    :param i_pf_tau_lookup:
    :param i_qt_m: Indices of the Qt controlled Branches with the tap module "m"
    :param n_qt_m:
    :param i_qt_m_lookup:
    :param i_qf_m: indices of the Qf controlled Branches with the tap module "m"
    :param n_qf_m:
    :param i_vt_m: Indices of the Vt controlled Branches with the tap module "m"
    :param n_vt_m:
    :param i_pf_dp: indices of the droop controlled Branches
    :param n_pf_dp:
    :param i_pf_dp_lookup:
    :param i_zero_beq: Indices of the Qf=0 controlled Branches with b_eq
    :param n_zero_beq:
    :param i_vf_beq: Indices of the Vf Controlled Branches with b_eq
    :param n_vf_beq:
    :param F: Array of "from" bus indices
    :param T: Array of "to" bus indices
    :param V: Array of complex bus voltages
    :param tap_modules_m: Array of tap modules
    :param k2: Array of branch converter losses
    :param tap_complex: Array of complex tap values {remember tap = m * exp(1j * tau) }
    :param Bc: Array of branch full susceptances
    :param b_eq: Array of branch equivalent (variable) susceptances
    :return:
    """

    nnz = 0
    p = 0
    Jp[p] = 0

    # column 1: derivatives w.r.t Va -----------------------------------------------------------------------------------
    for j in pvpq:  # sliced columns

        # J11
        if npvpq:
            for k in range(Yp[j], Yp[j + 1]):  # rows of A[:, j]

                # row index translation to the "rows" space
                i = Yi[k]
                ii = pvpq_lookup[i]

                if pvpq[ii] == i:
                    # entry found
                    Jx[nnz] = dSbus_dVa_x[k].real  # dP/dƟ
                    Ji[nnz] = ii
                    nnz += 1

        # J21 J31 J41
        offset = npvpq
        if ni2:
            for k in range(Yp[j], Yp[j + 1]):  # rows of A[:, j]

                # row index translation to the "rows" space
                i = Yi[k]
                ii = i2_lookup[i]

                if i2[ii] == i:
                    # entry found
                    Jx[nnz] = dSbus_dVa_x[k].imag  # dQ/dƟ
                    Ji[nnz] = ii + offset
                    nnz += 1


        # J51 J61
        offset += n_pf_tau
        if ni4:
            for k in range(dSf_dVa_p[j], dSf_dVa_p[j + 1]):  # rows of A[:, j]

                # row index translation to the "rows" space
                i = dSf_dVa_i[k]
                ii = i4_lookup[i]

                if i4[ii] == i:
                    # entry found
                    Jx[nnz] = dSf_dVa_x[k].imag   # dQf/dƟ
                    Ji[nnz] = ii + offset
                    nnz += 1

        # J71
        offset += ni4
        if n_qt_m:
            for k in range(dSt_dVa_p[j], dSt_dVa_p[j + 1]):  # rows of A[:, j]

                # row index translation to the "rows" space
                i = dSt_dVa_i[k]
                ii = i_qt_m_lookup[i]

                if i_qt_m[ii] == i:
                    # entry found
                    Jx[nnz] = dSt_dVa_x[k].imag   # dQt/dƟ
                    Ji[nnz] = ii + offset
                    nnz += 1

        # J81
        offset += ni2
        if n_pf_tau:
            for k in range(dSf_dVa_p[j], dSf_dVa_p[j + 1]):  # rows of A[:, j]

                # row index translation to the "rows" space
                i = dSf_dVa_i[k]
                ii = i_pf_tau_lookup[i]

                if i_pf_tau[ii] == i:
                    # entry found
                    Jx[nnz] = dSf_dVa_x[k].real  # dPf/dƟ
                    Ji[nnz] = ii + offset
                    nnz += 1

        # J91
        offset += n_qt_m
        if n_pf_dp:
            for k in range(dSf_dVa_p[j], dSf_dVa_p[j + 1]):  # rows of A[:, j]

                # row index translation to the "rows" space
                i = dSf_dVa_i[k]
                ii = i_pf_dp_lookup[i]

                if i_pf_dp[ii] == i:
                    # entry found
                    Jx[nnz] = -dSf_dVa_x[k].real    # dPf/dƟ
                    Ji[nnz] = ii + offset
                    nnz += 1

        # finalize column
        p += 1
        Jp[p] = nnz

    # column 2: derivatives w.r.t Vm -----------------------------------------------------------------------------------
    for j in pq:  # sliced columns

        # J12
        if npvpq:
            for k in range(Yp[j], Yp[j + 1]):  # rows of A[:, j]

                # row index translation to the "rows" space
                i = Yi[k]
                ii = pvpq_lookup[i]

                if pvpq[ii] == i:
                    # entry found
                    Jx[nnz] = dSbus_dVm_x[k].real  # dP/dVm
                    Ji[nnz] = ii
                    nnz += 1

        # J22 J32 J42
        offset = npvpq
        if ni2:
            for k in range(Yp[j], Yp[j + 1]):  # rows of A[:, j]

                # row index translation to the "rows" space
                i = Yi[k]
                ii = i2_lookup[i]

                if i2[ii] == i:
                    # entry found
                    Jx[nnz] = dSbus_dVm_x[k].imag  # dQ/dVm
                    Ji[nnz] = ii + offset
                    nnz += 1

        # J52 J62
        offset += n_pf_tau
        if ni4:
            for k in range(dSf_dVm_p[j], dSf_dVm_p[j + 1]):  # rows of A[:, j]

                # row index translation to the "rows" space
                i = dSf_dVm_i[k]
                ii = i4_lookup[i]

                if i4[ii] == i:
                    # entry found
                    Jx[nnz] = dSf_dVm_x[k].imag  # dQf/dVm
                    Ji[nnz] = ii + offset
                    nnz += 1

        # J72
        offset += ni4
        if n_qt_m:
            for k in range(dSt_dVm_p[j], dSt_dVm_p[j + 1]):  # rows of A[:, j]

                # row index translation to the "rows" space
                i = dSt_dVm_i[k]
                ii = i_qt_m_lookup[i]

                if i_qt_m[ii] == i:
                    # entry found
                    Jx[nnz] = dSt_dVm_x[k].imag  # dQt/dVm
                    Ji[nnz] = ii + offset
                    nnz += 1

        # J82
        offset += ni2
        if n_pf_tau:
            for k in range(dSf_dVm_p[j], dSf_dVm_p[j + 1]):  # rows of A[:, j]

                # row index translation to the "rows" space
                i = dSf_dVm_i[k]
                ii = i_pf_tau_lookup[i]

                if i_pf_tau[ii] == i:
                    # entry found
                    Jx[nnz] = dSf_dVm_x[k].real  # dPf/dVm
                    Ji[nnz] = ii + offset
                    nnz += 1

        # J92
        offset += n_qt_m
        if n_pf_dp:

            for k in range(dPfdp_dVm_p[j], dPfdp_dVm_p[j + 1]):  # rows of A[:, j]

                # row index translation to the "rows" space
                i = dPfdp_dVm_i[k]
                ii = i_pf_dp_lookup[i]

                if i_pf_dp[ii] == i:
                    # entry found
                    Jx[nnz] = dPfdp_dVm_x[k]  # dPfdp/dVm
                    Ji[nnz] = ii + offset
                    nnz += 1

        # finalize column
        p += 1
        Jp[p] = nnz

    # Column 3: derivatives w.r.t Beq for iBeqz + iBeqv ----------------------------------------------------------------
    if n_zero_beq + n_vf_beq:
        indices = np.concatenate((i_zero_beq, i_vf_beq))
        (dSbus_dBeq_data,
         dSbus_dBeq_indices,
         dSbus_dBeq_indptr,
         dSf_dBeqx_data,
         dSf_dBeqx_indices,
         dSf_dBeqx_indptr) = deriv.derivatives_Beq_csc_numba(indices, F, V, tap_modules_m, k2)

        for j in range(n_zero_beq + n_vf_beq):  # sliced columns

            # J13
            if npvpq:
                for k in range(dSbus_dBeq_indptr[j], dSbus_dBeq_indptr[j + 1]):  # rows of A[:, j]

                    # row index translation to the "rows" space
                    i = dSbus_dBeq_indices[k]
                    ii = pvpq_lookup[i]

                    if pvpq[ii] == i:
                        # entry found
                        Jx[nnz] = dSbus_dBeq_data[k].real  # dP/dBeq
                        Ji[nnz] = ii
                        nnz += 1

            # J23 J33 J43
            offset = npvpq
            if ni2:
                for k in range(dSbus_dBeq_indptr[j], dSbus_dBeq_indptr[j + 1]):  # rows of A[:, j]

                    # row index translation to the "rows" space
                    i = dSbus_dBeq_indices[k]
                    ii = i2_lookup[i]

                    if i2[ii] == i:
                        # entry found
                        Jx[nnz] = dSbus_dBeq_data[k].imag  # dQ/dBeq
                        Ji[nnz] = ii + offset
                        nnz += 1


            # J53 J63
            offset += n_pf_tau
            if ni4:
                for k in range(dSf_dBeqx_indptr[j], dSf_dBeqx_indptr[j + 1]):  # rows of A[:, j]

                    # row index translation to the "rows" space
                    i = dSf_dBeqx_indices[k]
                    ii = i4_lookup[i]

                    if i4[ii] == i:
                        # entry found
                        Jx[nnz] = dSf_dBeqx_data[k].imag  # dQf/dBeq
                        Ji[nnz] = ii + offset
                        nnz += 1

            # J73
            offset += ni4
            # if nQtma:  # --> The Jacobian is always zero :|
            #     for k in range(dSt_dBeqz.indptr[j], dSt_dBeqz.indptr[j + 1]):  # rows of A[:, j]
            #
            #         # row index translation to the "rows" space
            #         i = dSt_dBeqz.indices[k]
            #         ii = iQtma_lookup[i]
            #
            #         if iQtma[ii] == i:
            #             # entry found
            #             Jx[nnz] = dSt_dBeqz.data[k].imag
            #             Ji[nnz] = ii + offset
            #             nnz += 1


            # J83
            offset += ni2
            if n_pf_tau:
                for k in range(dSf_dBeqx_indptr[j], dSf_dBeqx_indptr[j + 1]):  # rows of A[:, j]

                    # row index translation to the "rows" space
                    i = dSf_dBeqx_indices[k]
                    ii = i_pf_tau_lookup[i]

                    if i_pf_tau[ii] == i:
                        # entry found
                        Jx[nnz] = dSf_dBeqx_data[k].real  # dPf/dBeq
                        Ji[nnz] = ii + offset
                        nnz += 1

            # J93
            offset += n_qt_m
            if n_pf_dp:
                for k in range(dSf_dBeqx_indptr[j], dSf_dBeqx_indptr[j + 1]):  # rows of A[:, j]

                    # row index translation to the "rows" space
                    i = dSf_dBeqx_indices[k]
                    ii = i_pf_dp_lookup[i]

                    if i_pf_dp[ii] == i:
                        # entry found
                        Jx[nnz] = -dSf_dBeqx_data[k].real  # dPf/dBeq
                        Ji[nnz] = ii + offset
                        nnz += 1

            # finalize column
            p += 1
            Jp[p] = nnz

    # Column 4: derivative w.r.t "m" for iQfma + iQfma + iVtma ---------------------------------------------------------
    if n_qf_m + n_qt_m + n_vt_m:
        indices = np.concatenate((i_qf_m, i_qt_m, i_vt_m))
        (dSbus_dm_data,
         dSbus_dm_indices,
         dSbus_dm_indptr,
         dSf_dm_data,
         dSf_dm_indices,
         dSf_dm_indptr,
         dSt_dm_data,
         dSt_dm_indices,
         dSt_dm_indptr) = deriv.derivatives_ma_csc_numba(indices, F, T, Ys, k2,
                                                         tap_complex, tap_modules_m, Bc, b_eq, V)

        for j in range(n_qf_m + n_qt_m + n_vt_m):  # sliced columns

            # J14
            if npvpq:
                for k in range(dSbus_dm_indptr[j], dSbus_dm_indptr[j + 1]):  # rows of A[:, j]

                    # row index translation to the "rows" space
                    i = dSbus_dm_indices[k]
                    ii = pvpq_lookup[i]

                    if pvpq[ii] == i:
                        # entry found
                        Jx[nnz] = dSbus_dm_data[k].real  # dP/dm
                        Ji[nnz] = ii
                        nnz += 1

            # J24 J34 J44
            offset = npvpq
            if ni2:
                for k in range(dSbus_dm_indptr[j], dSbus_dm_indptr[j + 1]):  # rows of A[:, j]

                    # row index translation to the "rows" space
                    i = dSbus_dm_indices[k]
                    ii = i2_lookup[i]

                    if i2[ii] == i:
                        # entry found
                        Jx[nnz] = dSbus_dm_data[k].imag  # dQ/dm
                        Ji[nnz] = ii + offset
                        nnz += 1

            # J54 J64
            offset += n_pf_tau
            if ni4:
                for k in range(dSf_dm_indptr[j], dSf_dm_indptr[j + 1]):  # rows of A[:, j]

                    # row index translation to the "rows" space
                    i = dSf_dm_indices[k]
                    ii = i4_lookup[i]

                    if i4[ii] == i:
                        # entry found
                        Jx[nnz] = dSf_dm_data[k].imag  # dQf/dm
                        Ji[nnz] = ii + offset
                        nnz += 1

            # J74
            offset += ni4
            if n_qt_m:
                for k in range(dSt_dm_indptr[j], dSt_dm_indptr[j + 1]):  # rows of A[:, j]

                    # row index translation to the "rows" space
                    i = dSt_dm_indices[k]
                    ii = i_qt_m_lookup[i]

                    if i_qt_m[ii] == i:
                        # entry found
                        Jx[nnz] = dSt_dm_data[k].imag  # dQt/dm
                        Ji[nnz] = ii + offset
                        nnz += 1

            # J84
            offset += ni2
            if n_pf_tau:
                for k in range(dSf_dm_indptr[j], dSf_dm_indptr[j + 1]):  # rows of A[:, j]

                    # row index translation to the "rows" space
                    i = dSf_dm_indices[k]
                    ii = i_pf_tau_lookup[i]

                    if i_pf_tau[ii] == i:
                        # entry found
                        Jx[nnz] = dSf_dm_data[k].real  # dPf/dm
                        Ji[nnz] = ii + offset
                        nnz += 1

            # J94
            offset += n_qt_m
            if n_pf_dp:
                for k in range(dSf_dm_indptr[j], dSf_dm_indptr[j + 1]):  # rows of A[:, j]

                    # row index translation to the "rows" space
                    i = dSf_dm_indices[k]
                    ii = i_pf_dp_lookup[i]

                    if i_pf_dp[ii] == i:
                        # entry found
                        Jx[nnz] = -dSf_dm_data[k].real  # dPf/dm
                        Ji[nnz] = ii + offset
                        nnz += 1

            # finalize column
            p += 1
            Jp[p] = nnz

    # Column 5: derivatives w.r.t theta sh for iPfsh + droop -----------------------------------------------------------
    if n_pf_tau + n_pf_dp > 0:

        indices = np.concatenate((i_pf_tau, i_pf_dp))
        (dSbus_dtau_data,
         dSbus_dtau_indices,
         dSbus_dtau_indptr,
         dSf_dtau_data,
         dSf_dtau_indices,
         dSf_dtau_indptr,
         dSt_dtau_data,
         dSt_dtau_indices,
         dSt_dtau_indptr) = deriv.derivatives_tau_csc_numba(indices, F, T, Ys, k2, tap_complex, V)

        for j in range(n_pf_tau + n_pf_dp):  # sliced columns

            # J15
            if npvpq:
                for k in range(dSbus_dtau_indptr[j], dSbus_dtau_indptr[j + 1]):  # rows of A[:, j]

                    # row index translation to the "rows" space
                    i = dSbus_dtau_indices[k]
                    ii = pvpq_lookup[i]

                    if pvpq[ii] == i:
                        # entry found
                        Jx[nnz] = dSbus_dtau_data[k].real  # dP/dtau
                        Ji[nnz] = ii
                        nnz += 1

            # J25 J35 J45
            offset = npvpq
            if ni2:
                for k in range(dSbus_dtau_indptr[j], dSbus_dtau_indptr[j + 1]):  # rows of A[:, j]

                    # row index translation to the "rows" space
                    i = dSbus_dtau_indices[k]
                    ii = i2_lookup[i]

                    if i2[ii] == i:
                        # entry found
                        Jx[nnz] = dSbus_dtau_data[k].imag  # dQ/dtau
                        Ji[nnz] = ii + offset
                        nnz += 1


            # J55 J65
            offset += n_pf_tau
            if ni4:
                for k in range(dSf_dtau_indptr[j], dSf_dtau_indptr[j + 1]):  # rows of A[:, j]

                    # row index translation to the "rows" space
                    i = dSf_dtau_indices[k]
                    ii = i4_lookup[i]

                    if i4[ii] == i:
                        # entry found
                        Jx[nnz] = dSf_dtau_data[k].imag  # dQf/dtau
                        Ji[nnz] = ii + offset
                        nnz += 1

            # J75
            offset += ni4
            if n_qt_m:
                for k in range(dSt_dtau_indptr[j], dSt_dtau_indptr[j + 1]):  # rows of A[:, j]

                    # row index translation to the "rows" space
                    i = dSt_dtau_indices[k]
                    ii = i_qt_m_lookup[i]

                    if i_qt_m[ii] == i:
                        # entry found
                        Jx[nnz] = dSt_dtau_data[k].imag  # dQt/dtau
                        Ji[nnz] = ii + offset
                        nnz += 1

            # J85
            offset += ni2
            if n_pf_tau:
                for k in range(dSf_dtau_indptr[j], dSf_dtau_indptr[j + 1]):  # rows of A[:, j]

                    # row index translation to the "rows" space
                    i = dSf_dtau_indices[k]
                    ii = i_pf_tau_lookup[i]

                    if i_pf_tau[ii] == i:
                        # entry found
                        Jx[nnz] = dSf_dtau_data[k].real  # dPf/dtau
                        Ji[nnz] = ii + offset
                        nnz += 1

            # J95
            offset += n_qt_m
            if n_pf_dp:
                for k in range(dSf_dtau_indptr[j], dSf_dtau_indptr[j + 1]):  # rows of A[:, j]

                    # row index translation to the "rows" space
                    i = dSf_dtau_indices[k]
                    ii = i_pf_dp_lookup[i]

                    if i_pf_dp[ii] == i:
                        # entry found
                        Jx[nnz] = -dSf_dtau_data[k].real  # - dPf/dtau
                        Ji[nnz] = ii + offset
                        nnz += 1

            # finalize column
            p += 1
            Jp[p] = nnz

    # Finalize ----------------------------------------------------------------------------
    #  finalize the Jacobian Pointer
    Jp[p] = nnz

    return nnz, Jx, Ji, Jp


def fubm_jacobian(nb, nl, iPfsh, iPfdp, iQfma, iQtma, iVtma, iBeqz, iBeqv, VfBeqbus, Vtmabus,
                  F, T, Ys, k2, tap, ma, Bc, Beq, Kdp, V, Ybus, Yf, Yt, Cf, Ct, pvpq, pq):
    """
    Compute the FUBM jacobian in a dynamic fashion by only computing the derivatives that are needed
    :param nb: number of buses
    :param nl: Number of lines
    :param iPfsh: indices of the Pf controlled with the shunt susceptance Branches
    :param iPfdp: indices of the Pf-droop controlled Branches
    :param iQfma: indices of the Qf controlled with ma Branches
    :param iQtma: Indices of the Qt controlled with ma Branches
    :param iVtma: Indices of the Vt controlled with ma Branches
    :param iBeqz: Indices of the Qf made zero with the equivalent susceptance Branches
    :param iBeqv: Indices of the Vf Controlled with the equivalent susceptance Branches
    :param F: Array of "from" bus indices
    :param T: Array of "to" bus indices
    :param Ys: Array of branch series admittances
    :param k2: Array of branch converter losses
    :param tap: Array of complex tap values {remember tap = ma * exp(1j * theta) }
    :param ma: Array of tap modules
    :param Bc: Array of branch full susceptances
    :param Beq: Array of branch equivalent (variable) susceptances
    :param Kdp: Array of branch converter droop constants
    :param V: Array of complex bus voltages
    :param Ybus: Admittance matrix
    :param Yf: Admittances matrix of the Branches with the "from" buses
    :param Yt: Admittances matrix of the Branches with the "to" buses
    :param Cf: Connectivity matrix of the Branches with the "from" buses
    :param Ct: Connectivity matrix of the Branches with the "to" buses
    :param pvpq: Array of pv and then pq bus indices (not sorted)
    :param pq: Array of PQ bus indices
    :return: FUBM Jacobian matrix
    """
    nPfsh = len(iPfsh)
    nPfdp = len(iPfdp)
    nQfma = len(iQfma)
    nQtma = len(iQtma)
    nVtma = len(iVtma)
    nBeqz = len(iBeqz)
    nBeqv = len(iBeqv)
    nVfBeqbus = len(VfBeqbus)
    nVtmabus = len(Vtmabus)
    npq = len(pq)
    npvpq = len(pvpq)
    nbus = Ybus.shape[0]
    nbr = Yf.shape[0]

    i2 = np.r_[pq, VfBeqbus, Vtmabus]
    i4 = np.r_[iQfma, iBeqz]
    ni2 = len(i2)
    ni4 = len(i4)
    E = V / np.abs(V)

    # compose the derivatives of the power Injections w.r.t Va and Vm
    dSbus_dVm_x, dSbus_dVa_x = deriv.dSbus_dV_numba_sparse_csc(Ybus.data, Ybus.indptr, Ybus.indices, V, E)

    # compose the derivatives of the branch flow w.r.t Va and Vm
    dSf_dVm, dSf_dVa = deriv.dSf_dV_csc(Yf, V, F, T)

    if nQtma:
        dSt_dVm, dSt_dVa = deriv.dSt_dV_csc(Yt, V, F, T)
    else:
        dSt_dVa = csc_matrix((nl, nb))
        dSt_dVm = csc_matrix((nl, nb))

    if nPfdp:
        # compute the droop derivative
        dVmf_dVm = lil_matrix((nl, nb))
        dVmf_dVm[iPfdp, :] = Cf[iPfdp, :]
        dPfdp_dVm = -dSf_dVm.real + diags(Kdp) * dVmf_dVm
    else:
        dPfdp_dVm = csc_matrix((nPfdp, nb))

    n_cols = npvpq + npq + nBeqz + nBeqv + nQfma + nQtma + nVtma + nPfsh + nPfdp
    n_rows = n_cols

    nnz_estimate = Ybus.nnz * 8
    Jx = np.empty(nnz_estimate, dtype=np.float64)  # data
    Ji = np.empty(nnz_estimate, dtype=np.int32)  # indices
    Jp = np.empty(n_cols + 1, dtype=np.int32)  # pointers

    # generate lookup for the row slicing
    pvpq_lookup = make_lookup(nbus, pvpq)
    i2_lookup = make_lookup(nbus, i2)
    iPfsh_lookup = make_lookup(nbr, iPfsh)
    i4_lookup = make_lookup(nbr, i4)
    iQtma_lookup = make_lookup(nbr, iQtma)
    iPfdp_lookup = make_lookup(nbr, iPfdp)

    # fill the jacobian data with numba
    nnz, Jx, Ji, Jp = fill_acdc_jacobian_data(Jx, Ji, Jp, Ybus.indptr, Ybus.indices, Ys,
                                              dSbus_dVa_x, dSbus_dVm_x,
                                              dSf_dVa.data, dSf_dVa.indices, dSf_dVa.indptr,
                                              dSf_dVm.data, dSf_dVm.indices, dSf_dVm.indptr,
                                              dSt_dVa.data, dSt_dVa.indices, dSt_dVa.indptr,
                                              dSt_dVm.data, dSt_dVm.indices, dSt_dVm.indptr,
                                              dPfdp_dVm.data, dPfdp_dVm.indices, dPfdp_dVm.indptr,
                                              pvpq, pvpq_lookup, npvpq, pq,
                                              i2, ni2, i2_lookup, i4, ni4, i4_lookup,
                                              iPfsh, nPfsh, iPfsh_lookup,
                                              iQtma, nQtma, iQtma_lookup,
                                              iQfma, nQfma,
                                              iVtma, nVtma,
                                              iPfdp, nPfdp, iPfdp_lookup,
                                              iBeqz, nBeqz,
                                              iBeqv, nBeqv,
                                              F, T, V, ma, k2, tap, Bc, Beq)

    Jx = np.resize(Jx, nnz)
    Ji = np.resize(Ji, nnz)

    J = csc_matrix((Jx, Ji, Jp), shape=(n_rows, n_cols))

    return J


