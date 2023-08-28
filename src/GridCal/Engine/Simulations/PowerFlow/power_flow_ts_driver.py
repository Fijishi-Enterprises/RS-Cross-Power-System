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
import time
from typing import Union
from GridCal.Engine.Simulations.PowerFlow.power_flow_ts_results import PowerFlowTimeSeriesResults
from GridCal.Engine.Core.Devices.multi_circuit import MultiCircuit
from GridCal.Engine.Simulations.PowerFlow.power_flow_options import PowerFlowOptions
from GridCal.Engine.Simulations.driver_types import SimulationTypes
from GridCal.Engine.Simulations.driver_template import TimeSeriesDriverTemplate
from GridCal.Engine.Simulations.Clustering.clustering_results import ClusteringResults
import GridCal.Engine.Simulations.PowerFlow.power_flow_worker as pf_worker
from GridCal.Engine.Core.Compilers.circuit_to_bentayga import bentayga_pf
from GridCal.Engine.Core.Compilers.circuit_to_newton_pa import newton_pa_pf
from GridCal.Engine.Core.Compilers.circuit_to_pgm import pgm_pf
import GridCal.Engine.basic_structures as bs


class PowerFlowTimeSeriesDriver(TimeSeriesDriverTemplate):
    tpe = SimulationTypes.TimeSeries_run
    name = tpe.value

    def __init__(self, grid: MultiCircuit,
                 options: PowerFlowOptions,
                 time_indices: np.ndarray,
                 opf_time_series_results=None,
                 clustering_results: Union[ClusteringResults, None] = None,
                 engine: bs.EngineType = bs.EngineType.GridCal):
        """
        PowerFlowTimeSeries constructor
        :param grid: MultiCircuit instance
        :param options: PowerFlowOptions instance
        :param time_indices: array of time indices to simulate
        :param opf_time_series_results: ClusteringResults instance (optional)
        :param clustering_results: ClusteringResults instance (optional)
        :param engine: Calculation engine to use
        """
        TimeSeriesDriverTemplate.__init__(self,
                                          grid=grid,
                                          time_indices=time_indices,
                                          clustering_results=clustering_results,
                                          engine=engine)

        # reference the grid directly
        # self.grid = grid

        self.options = options

        self.opf_time_series_results = opf_time_series_results

    def run_single_thread(self, time_indices) -> PowerFlowTimeSeriesResults:
        """
        Run single thread time series
        :param time_indices: array of time indices to consider
        :return: TimeSeriesResults instance
        """

        n = self.grid.get_bus_number()
        m = self.grid.get_branch_number_wo_hvdc()

        # initialize the grid time series results we will append the island results with another function
        time_series_results = PowerFlowTimeSeriesResults(n=n,
                                                         m=m,
                                                         n_hvdc=self.grid.get_hvdc_number(),
                                                         bus_names=self.grid.get_bus_names(),
                                                         branch_names=self.grid.get_branch_names_wo_hvdc(),
                                                         hvdc_names=self.grid.get_hvdc_names(),
                                                         bus_types=np.zeros(m),
                                                         time_array=self.grid.time_profile[time_indices],
                                                         clustering_results=self.clustering_results)

        # compile dictionaries once for speed
        bus_dict = {bus: i for i, bus in enumerate(self.grid.buses)}
        areas_dict = {elm: i for i, elm in enumerate(self.grid.areas)}
        self.progress_signal.emit(0.0)
        for it, t in enumerate(time_indices):

            self.progress_text.emit('Time series at ' + str(self.grid.time_profile[t]) + '...')

            progress = ((it + 1) / len(time_indices)) * 100
            self.progress_signal.emit(progress)

            pf_res = pf_worker.multi_island_pf(multi_circuit=self.grid,
                                               t=t,
                                               options=self.options,
                                               opf_results=self.opf_time_series_results,
                                               bus_dict=bus_dict,
                                               areas_dict=areas_dict)

            # gather results
            time_series_results.voltage[it, :] = pf_res.voltage
            time_series_results.S[it, :] = pf_res.Sbus
            time_series_results.Sf[it, :] = pf_res.Sf
            time_series_results.St[it, :] = pf_res.St
            time_series_results.Vbranch[it, :] = pf_res.Vbranch
            time_series_results.loading[it, :] = pf_res.loading
            time_series_results.losses[it, :] = pf_res.losses
            time_series_results.hvdc_losses[it, :] = pf_res.hvdc_losses
            time_series_results.hvdc_Pf[it, :] = pf_res.hvdc_Pf
            time_series_results.hvdc_Pt[it, :] = pf_res.hvdc_Pt
            time_series_results.hvdc_loading[it, :] = pf_res.hvdc_loading
            time_series_results.error_values[it] = pf_res.error
            time_series_results.converged_values[it] = pf_res.converged

            if self.__cancel__:
                return time_series_results

        return time_series_results

    def run_bentayga(self):

        res = bentayga_pf(self.grid, self.options, time_series=True)

        results = PowerFlowTimeSeriesResults(n=self.grid.get_bus_number(),
                                             m=self.grid.get_branch_number_wo_hvdc(),
                                             n_hvdc=self.grid.get_hvdc_number(),
                                             bus_names=res.names,
                                             branch_names=res.names,
                                             hvdc_names=res.hvdc_names,
                                             bus_types=res.bus_types,
                                             time_array=self.grid.time_profile,
                                             clustering_results=self.clustering_results)

        results.voltage = res.V
        results.S = res.S
        results.Sf = res.Sf
        results.St = res.St
        results.loading = res.loading
        results.losses = res.losses
        results.Vbranch = res.Vbranch
        results.If = res.If
        results.It = res.It
        results.Beq = res.Beq
        results.m = res.tap_modules
        results.tap_angle = res.tap_angles

        return results

    def run_newton_pa(self, time_indices=None) -> PowerFlowTimeSeriesResults:
        """
        Run with Newton Power Analytics
        :param time_indices: array of time indices
        :return:
        """
        res = newton_pa_pf(circuit=self.grid,
                           pf_opt=self.options,
                           time_series=True,
                           time_indices=time_indices)

        results = PowerFlowTimeSeriesResults(n=self.grid.get_bus_number(),
                                             m=self.grid.get_branch_number_wo_hvdc(),
                                             n_hvdc=self.grid.get_hvdc_number(),
                                             bus_names=res.bus_names,
                                             branch_names=res.branch_names,
                                             hvdc_names=res.hvdc_names,
                                             bus_types=res.bus_types,
                                             time_array=self.grid.time_profile[time_indices],
                                             clustering_results=self.clustering_results)

        results.voltage = res.voltage
        results.S = res.Scalc
        results.Sf = res.Sf
        results.St = res.St
        results.loading = res.Loading
        results.losses = res.Losses
        # results.Vbranch = res.Vbranch
        # results.If = res.If
        # results.It = res.It
        results.Beq = res.Beq
        results.tap_module = res.tap_module
        results.tap_angle = res.tap_angle
        results.F = res.F
        results.T = res.T
        results.hvdc_F = res.hvdc_F
        results.hvdc_T = res.hvdc_T
        results.hvdc_Pf = res.hvdc_Pf
        results.hvdc_Pt = res.hvdc_Pt
        results.hvdc_loading = res.hvdc_loading
        results.hvdc_losses = res.hvdc_losses
        results.error_values = res.error

        return results

    def run(self):
        """
        Run the time series simulation
        @return:
        """

        self.tic()

        if self.engine == bs.EngineType.GridCal:
            self.results = self.run_single_thread(time_indices=self.time_indices)

        elif self.engine == bs.EngineType.Bentayga:
            self.progress_text.emit('Running Bentayga... ')
            self.results = self.run_bentayga()

        elif self.engine == bs.EngineType.NewtonPA:
            self.progress_text.emit('Running Newton power analytics... ')
            self.results = self.run_newton_pa(time_indices=self.time_indices)

        elif self.engine == bs.EngineType.PGM:
            self.progress_text.emit('Running Power Grid Model... ')
            self.results = pgm_pf(self.grid, self.options, logger=self.logger, time_series=True)
            self.results.area_names = [a.name for a in self.grid.areas]

        else:
            raise Exception('Engine not implemented for Time Series:' + self.engine.value)

        # fill F, T, Areas, etc...
        self.results.fill_circuit_info(self.grid)

        self.toc()
