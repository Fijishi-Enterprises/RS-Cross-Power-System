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

import time
from typing import Union
from GridCal.Engine.Core.Devices.multi_circuit import MultiCircuit
from GridCal.Engine.Simulations.driver_types import SimulationTypes
from GridCal.Engine.Simulations.driver_template import DriverTemplate
from GridCal.Engine.Simulations.Clustering.clustering_results import ClusteringResults
from GridCal.Engine.Simulations.Clustering.clustering_options import ClusteringAnalysisOptions
from GridCal.Engine.Simulations.Clustering.clustering import kmeans_sampling


class ClusteringDriver(DriverTemplate):
    name = 'Clustering analysis'
    tpe = SimulationTypes.ClusteringAnalysis_run

    def __init__(self, grid: MultiCircuit, options: ClusteringAnalysisOptions):
        """
        Clustering analysis driver constructor
        :param grid: Multicircuit instance
        :param options: ClusteringAnalysisOptions
        """
        DriverTemplate.__init__(self, grid=grid)

        self.options: ClusteringAnalysisOptions = options

        self.results: Union[ClusteringResults, None] = None

    def run(self):
        """
        Run thread
        """
        start = time.time()
        self.progress_text.emit('Clustering')
        self.progress_signal.emit(0)

        time_indices, sampled_probabilities, sample_idx = kmeans_sampling(x_input=self.grid.get_Pbus_prof(),
                                                                          n_points=self.options.n_points)
        self.results = ClusteringResults(
            time_indices=time_indices,
            sampled_probabilities=sampled_probabilities,
            time_array=self.grid.time_profile,
            original_sample_idx=sample_idx
        )

        end = time.time()
        self.elapsed = end - start

    def get_steps(self):
        """
        Get variations list of strings
        """
        if self.results is not None:
            return [self.grid.time_profile[i].strftime('%d-%m-%Y %H:%M') for i in self.results.time_indices]
        else:
            return list()
