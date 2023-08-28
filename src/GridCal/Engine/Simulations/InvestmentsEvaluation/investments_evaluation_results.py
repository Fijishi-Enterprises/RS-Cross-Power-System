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
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.colors as plt_colors
from GridCal.Engine.Simulations.results_template import ResultsTemplate
from GridCal.Engine.Simulations.result_types import ResultTypes
from GridCal.Engine.Simulations.results_table import ResultsTable
from GridCal.Engine.Core.Devices.multi_circuit import MultiCircuit
from GridCal.Engine.basic_structures import Vec, IntVec, StrVec


class InvestmentsEvaluationResults(ResultsTemplate):
    tpe = 'Investments Evaluation Results'

    def __init__(self, grid: MultiCircuit, max_eval: int):
        """
        Construct the analysis
        :param grid: MultiCircuit
        :param max_eval: maximum number of evaluations
        """
        available_results = {
            ResultTypes.ReportsResults: [ResultTypes.InvestmentsReportResults, ],
            ResultTypes.SpecialPlots: [ResultTypes.InvestmentsParetoPlot,
                                       ResultTypes.InvestmentsIterationsPlot]
        }

        ResultsTemplate.__init__(self,
                                 name='Investments Evaluation',
                                 available_results=available_results,
                                 data_variables=[],
                                 time_array=None,
                                 clustering_results=None)

        self.grid = grid
        self.n_groups = len(grid.investments_groups)
        self.max_eval = max_eval

        self._combinations: IntVec = np.zeros((max_eval, self.n_groups), dtype=int)
        self._capex: Vec = np.zeros(max_eval, dtype=float)
        self._opex: Vec = np.zeros(max_eval, dtype=float)
        self._losses: Vec = np.zeros(max_eval, dtype=float)
        self._overload_score: Vec = np.zeros(max_eval, dtype=float)
        self._voltage_score: Vec = np.zeros(max_eval, dtype=float)
        self._losses: Vec = np.zeros(max_eval, dtype=float)
        self._f_obj: Vec = np.zeros(max_eval, dtype=float)
        self._index_names: Vec = np.zeros(max_eval, dtype=object)

    def get_index(self) -> StrVec:
        """
        Return index names
        """
        return self._index_names

    def set_at(self, eval_idx, capex, opex, losses, overload_score, voltage_score, objective_function,
               combination: IntVec, index_name: str) -> None:
        """
        Set the results at an investment group
        :param eval_idx: evaluation index
        :param capex:
        :param opex:
        :param losses:
        :param overload_score:
        :param voltage_score:
        :param objective_function:
        :param combination: vector of size (n_investment_groups) with ones in those investments used
        :param index_name: Name of the evaluation
        """
        self._capex[eval_idx] = capex
        self._opex[eval_idx] = opex
        self._losses[eval_idx] = losses
        self._overload_score[eval_idx] = overload_score
        self._voltage_score[eval_idx] = voltage_score
        self._f_obj[eval_idx] = objective_function
        self._combinations[eval_idx, :] = combination
        self._index_names[eval_idx] = index_name

    def mdl(self, result_type) -> "ResultsTable":
        """
        Plot the results
        :param result_type: type of results (string)
        :return: DataFrame of the results
                (or None if the result was not understood)
        """

        if result_type == ResultTypes.InvestmentsReportResults:
            labels = self._index_names
            columns = ["CAPEX (M€)",
                       "OPEX (M€/yr)",
                       "Losses (MW)",
                       "Overload cost (€)",
                       "Voltage deviations cost (€)",
                       "Objective function"] + self.grid.get_investment_groups_names()
            data = np.c_[self._capex,
                         self._opex,
                         self._losses,
                         self._overload_score,
                         self._voltage_score,
                         self._f_obj,
                         self._combinations]
            y_label = ''
            title = ''

        elif result_type == ResultTypes.InvestmentsParetoPlot:
            labels = self._index_names
            columns = ["CAPEX (M€) + OPEX (M€)", "Objective function"]
            x = self._capex + self._opex
            y = self._f_obj
            data = np.c_[x, y]
            y_label = ''
            title = ''

            plt.ion()
            color_norm = plt_colors.LogNorm()
            fig = plt.figure(figsize=(8, 6))
            ax3 = plt.subplot(1, 1, 1)
            sc3 = ax3.scatter(x, y, c=y, norm=color_norm)
            ax3.set_xlabel('Investment cost (M€)')
            ax3.set_ylabel('Total cost of losses (M€)')
            plt.colorbar(sc3, fraction=0.05, label='Objective function')
            fig.suptitle(result_type.value[0])
            plt.tight_layout()
            plt.show()

        elif result_type == ResultTypes.InvestmentsIterationsPlot:
            labels = self._index_names
            columns = ["CAPEX (M€) + OPEX (M€)", "Objective function"]
            x = np.arange(self.max_eval)
            y = self._f_obj
            data = np.c_[x, y]
            y_label = ''
            title = ''

            plt.ion()
            fig = plt.figure(figsize=(8, 6))
            ax3 = plt.subplot(1, 1, 1)
            ax3.plot(x, y, '.')
            # plt.plot(iters, self.best_y[0:self.iter], 'r')
            ax3.set_xlabel('Iteration')
            ax3.set_ylabel('Objective')
            fig.suptitle(result_type.value[0])
            plt.grid()
            plt.show()

        else:
            columns = []
            labels = []
            data = np.zeros(0)
            y_label = '(MW)'
            title = ''

        mdl = ResultsTable(data=data,
                           index=np.array(labels),
                           columns=np.array(columns),
                           title=title,
                           ylabel=y_label,
                           xlabel='',
                           units=y_label)
        return mdl
