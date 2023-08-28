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
from GridCal.Engine.basic_structures import DateVec, IntVec, StrVec, Vec, Mat, CxVec
from GridCal.Engine.Simulations.result_types import ResultTypes
from GridCal.Engine.Simulations.results_table import ResultsTable
from GridCal.Engine.Simulations.results_template import ResultsTemplate
from GridCal.Engine.Core.DataStructures.numerical_circuit import NumericalCircuit


class LinearAnalysisTimeSeriesResults(ResultsTemplate):

    def __init__(
            self,
            n: int,
            m: int,
            time_array: DateVec,
            bus_names: StrVec,
            bus_types: IntVec,
            branch_names: StrVec,
            clustering_results):
        """
        Constructor
        :param n: number of buses
        :param m: number of Branches
        :param time_array: array of time steps
        :param bus_names: array of bus names
        :param bus_types: array of bus types
        :param branch_names: array of branch names
        """
        ResultsTemplate.__init__(
            self,
            name='Linear Analysis time series',
            available_results=[
                ResultTypes.BusActivePower,
                ResultTypes.BranchActivePowerFrom,
                ResultTypes.BranchLoading
            ],
            data_variables=[
                'bus_names',
                'bus_types',
                'time',
                'branch_names',
                'voltage',
                'S',
                'Sf',
                'loading',
                'losses'
            ],
            time_array=time_array,
            clustering_results=clustering_results
        )

        self.nt: int = len(time_array)
        self.m: int = m
        self.n: int = n
        # self.time_array: DateVec = time_array

        self.bus_names: StrVec = bus_names

        self.bus_types: IntVec = bus_types

        self.branch_names: StrVec = branch_names

        self.voltage: CxVec = np.ones((self.nt, n), dtype=complex)

        self.S: CxVec = np.zeros((self.nt, n), dtype=complex)

        self.Sf: CxVec = np.zeros((self.nt, m), dtype=complex)

        self.loading: Vec = np.zeros((self.nt, m), dtype=float)

        self.losses: CxVec = np.zeros((self.nt, m), dtype=float)

    def apply_new_time_series_rates(self, nc: NumericalCircuit) -> Mat:
        rates = nc.Rates.T
        self.loading = self.Sf / (rates + 1e-9)

    def get_results_dict(self):
        """
        Returns a dictionary with the results sorted in a dictionary
        :return: dictionary of 2D numpy arrays (probably of complex numbers)
        """
        data = {
            'V': self.voltage.tolist(),
            'P': self.S.real.tolist(),
            'Q': self.S.imag.tolist(),
            'Sbr_real': self.Sf.real.tolist(),
            'Sbr_imag': self.Sf.imag.tolist(),
            'loading': np.abs(self.loading).tolist()
        }
        return data

    def mdl(self, result_type: ResultTypes) -> ResultsTable:
        """
        Get ResultsModel instance
        :param result_type:
        :return: ResultsModel instance
        """

        if result_type == ResultTypes.BusActivePower:
            labels = self.bus_names
            data = self.S.real
            y_label = '(MW)'
            title = 'Bus active power '

        elif result_type == ResultTypes.BranchActivePowerFrom:
            labels = self.branch_names
            data = self.Sf.real
            y_label = '(MW)'
            title = 'Branch power '

        elif result_type == ResultTypes.BranchLoading:
            labels = self.branch_names
            data = np.abs(self.loading) * 100
            y_label = '(%)'
            title = 'Branch loading '

        elif result_type == ResultTypes.BranchLosses:
            labels = self.branch_names
            data = self.losses
            y_label = '(MVA)'
            title = 'Branch losses'

        elif result_type == ResultTypes.BusVoltageModule:
            labels = self.bus_names
            data = np.abs(self.voltage)
            y_label = '(p.u.)'
            title = 'Bus voltage'

        else:
            raise Exception('Result type not understood:' + str(result_type))

        if self.time_array is not None:
            index = self.time_array
        else:
            index = list(range(data.shape[0]))

        # assemble model
        return ResultsTable(
            data=data,
            index=index,
            columns=labels,
            title=title,
            ylabel=y_label,
            units=y_label
        )