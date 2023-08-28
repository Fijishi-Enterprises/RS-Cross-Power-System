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
from PySide6 import QtWidgets
from matplotlib import pyplot as plt

import GridCal.Engine.Simulations as sim
import GridCal.Gui.GuiFunctions as gf
from GridCal.Gui.messages import error_msg, warning_msg
from GridCal.Gui.Main.gui_sub_classes.simulations import SimulationsMain
from GridCal.Gui.Session.session import ResultsModel


class ResultsMain(SimulationsMain):
    """
    Diagrams Main
    """

    def __init__(self, parent=None):
        """

        @param parent:
        """

        # create main window
        SimulationsMain.__init__(self, parent)

        self.results_mdl: sim.ResultsTable = sim.ResultsTable(data=np.zeros((0, 0)),
                                                              columns=np.zeros(0),
                                                              index=np.zeros(0))

        # --------------------------------------------------------------------------------------------------------------
        self.ui.actionSet_OPF_generation_to_profiles.triggered.connect(self.copy_opf_to_profiles)

        # Buttons
        self.ui.saveResultsButton.clicked.connect(self.save_results_df)
        self.ui.copy_results_pushButton.clicked.connect(self.copy_results_data)
        self.ui.copy_numpy_button.clicked.connect(self.copy_results_data_as_numpy)
        self.ui.plot_data_pushButton.clicked.connect(self.plot_results)
        self.ui.search_results_Button.clicked.connect(self.search_in_results)
        self.ui.deleteDriverButton.clicked.connect(self.delete_results_driver)

        # tree-click
        self.ui.results_treeView.clicked.connect(self.results_tree_view_click)

        # line edit enter
        self.ui.sear_results_lineEdit.returnPressed.connect(self.search_in_results)

    def results_tree_view_click(self, index):
        """
        Display the simulation results on the results table
        """
        tree_mdl = self.ui.results_treeView.model()
        item = tree_mdl.itemFromIndex(index)
        path = gf.get_tree_item_path(item)

        if len(path) > 1:

            if len(path) == 2:
                study_name = path[0]
                result_name = path[1]
            elif len(path) == 3:
                study_name = path[0]
                result_name = path[2]
            else:
                raise Exception('Path len ' + str(len(path)) + ' not supported')

            if study_name in self.available_results_dict.keys():
                if result_name in self.available_results_dict[study_name].keys():

                    study_type = self.available_results_dict[study_name][result_name]

                    self.results_mdl = None

                    self.results_mdl = self.session.get_results_model_by_name(study_name=study_name,
                                                                              study_type=study_type)

                    if self.results_mdl is not None:

                        if self.ui.results_as_abs_checkBox.isChecked():
                            self.results_mdl.convert_to_abs()

                        if self.ui.results_as_cdf_checkBox.isChecked():
                            self.results_mdl.convert_to_cdf()

                        # set the table model
                        self.ui.resultsTableView.setModel(self.results_mdl)
                        self.ui.units_label.setText(self.results_mdl.units)
                    else:
                        self.ui.resultsTableView.setModel(None)
                        self.ui.units_label.setText("")

                else:
                    self.ui.resultsTableView.setModel(None)
                    self.ui.units_label.setText("")
            else:
                self.ui.resultsTableView.setModel(None)
                self.ui.units_label.setText("")

    def plot_results(self):
        """
        Plot the results
        """
        mdl: ResultsModel = self.ui.resultsTableView.model()

        if mdl is not None:

            plt.rcParams["date.autoformatter.minute"] = "%Y-%m-%d %H:%M:%S"

            # get the selected element
            obj_idx = self.ui.resultsTableView.selectedIndexes()

            # create figure to plot
            fig = plt.figure(figsize=(12, 8))
            ax = fig.add_subplot(111)

            if len(obj_idx):

                # get the unique columns in the selected cells
                cols = np.zeros(len(obj_idx), dtype=int)
                rows = np.zeros(len(obj_idx), dtype=int)

                for i in range(len(obj_idx)):
                    cols[i] = obj_idx[i].column()
                    rows[i] = obj_idx[i].row()

                cols = np.unique(cols)
                rows = np.unique(rows)

            else:
                # plot all
                cols = None
                rows = None

            # none selected, plot all
            mdl.plot(ax=ax, selected_col_idx=cols, selected_rows=rows, stacked=False)

            plt.show()

    def save_results_df(self):
        """
        Save the data displayed at the results as excel
        """
        mdl: ResultsModel = self.ui.resultsTableView.model()

        if mdl is not None:
            file, filter_ = QtWidgets.QFileDialog.getSaveFileName(self, "Export results", '',
                                                                  filter="CSV (*.csv);;Excel files (*.xlsx)")

            if file != '':
                if 'xlsx' in filter_:
                    f = file
                    if not f.endswith('.xlsx'):
                        f += '.xlsx'
                    mdl.save_to_excel(f)
                    print('Saved!')
                if 'csv' in filter_:
                    f = file
                    if not f.endswith('.csv'):
                        f += '.csv'
                    mdl.save_to_csv(f)
                    print('Saved!')
                else:
                    error_msg(file[0] + ' is not valid :(')
        else:
            warning_msg('There is no profile displayed, please display one', 'Copy profile to clipboard')

    def copy_results_data(self):
        """
        Copy the current displayed profiles to the clipboard
        """
        mdl = self.ui.resultsTableView.model()
        if mdl is not None:
            mdl.copy_to_clipboard()
            print('Copied!')
        else:
            warning_msg('There is no profile displayed, please display one', 'Copy profile to clipboard')

    def copy_results_data_as_numpy(self):
        """
        Copy the current displayed profiles to the clipboard
        """
        mdl = self.ui.resultsTableView.model()
        if mdl is not None:
            mdl.copy_numpy_to_clipboard()
            print('Copied!')
        else:
            warning_msg('There is no profile displayed, please display one', 'Copy profile to clipboard')

    def search_in_results(self):
        """
        Search in the results model
        """

        if self.results_mdl is not None:
            text = self.ui.sear_results_lineEdit.text().strip()

            if text != '':
                mdl = self.results_mdl.search(text)
            else:
                mdl = None

            self.ui.resultsTableView.setModel(mdl)

    def delete_results_driver(self):
        """
        Delete the driver
        :return:
        """
        idx = self.ui.results_treeView.selectedIndexes()
        if len(idx) > 0:
            tree_mdl = self.ui.results_treeView.model()
            item = tree_mdl.itemFromIndex(idx[0])
            path = gf.get_tree_item_path(item)

            if len(path) > 0:
                study_name = path[0]
                study_type = self.available_results_dict[study_name]

                quit_msg = "Do you want to delete the results driver " + study_name + "?"
                reply = QtWidgets.QMessageBox.question(self, 'Message',
                                                       quit_msg,
                                                       QtWidgets.QMessageBox.StandardButton.Yes,
                                                       QtWidgets.QMessageBox.StandardButton.No)

                if reply == QtWidgets.QMessageBox.StandardButton.Yes.value:
                    self.session.delete_driver_by_name(study_name)
                    self.update_available_results()

    def copy_opf_to_profiles(self):
        """
        Copy the results from the OPF time series to the profiles
        """
        drv, results = self.session.get_driver_results(sim.SimulationTypes.OPFTimeSeries_run)

        if results is not None:

            reply = QtWidgets.QMessageBox.question(self, 'Message',
                                                   'Are you sure that you want to overwrite '
                                                   'the generation profiles with the OPF results?',
                                                   QtWidgets.QMessageBox.StandardButton.Yes,
                                                   QtWidgets.QMessageBox.StandardButton.No)

            if reply == QtWidgets.QMessageBox.StandardButton.Yes.value:
                for i, gen in enumerate(self.circuit.get_generators()):
                    gen.P_prof = results.generator_power[:, i]

        else:
            warning_msg('The OPF time series has no results :(')
