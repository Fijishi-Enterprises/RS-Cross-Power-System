

# from PySide2 import QtCore, QtGui, QtWidgets

from GridCal.Gui.ProfilesInput.excel_sheet_selection import *
import sys
import xlrd


class ExcelDialog(QtWidgets.QDialog):

    def __init__(self, excel_file=None):

        QtWidgets.QDialog.__init__(self)
        self.ui = Ui_ExcelSelectionDialog()
        self.ui.setupUi(self)

        # click
        self.ui.acceptButton.clicked.connect(self.accepted)
        self.ui.cancelButton.clicked.connect(self.rejected)
        self.ui.sheets_list.doubleClicked.connect(self.accepted)

        self.excel_sheet = None

        self.sheet_names = list()

        # self.fill_from_file(excel_file=excel_file)

    def fill_from_file(self, excel_file):
        if excel_file is not None:
            xls = xlrd.open_workbook(excel_file, on_demand=True)
            self.sheet_names = xls.sheet_names()
            self.ui.sheets_list.addItems(self.sheet_names)

            if len(self.sheet_names) > 0:
                self.excel_sheet = 0

    def accepted(self):
        """

        :return:
        """
        if len(self.ui.sheets_list.selectedIndexes()):
            self.excel_sheet = self.ui.sheets_list.selectedIndexes()[0].row()
        print('Accepted: self.excel_sheet: ', self.excel_sheet)

        self.close()

    def rejected(self):
        """

        :return:
        """
        print('Rejected: self.excel_sheet: ', self.excel_sheet)
        self.close()


def excel_dialogue_show(excel_file):

    # app = QtWidgets.QApplication()
    ExcelSelectionDialog = QtWidgets.QDialog()
    ui = Ui_ExcelSelectionDialog()
    ui.setupUi(ExcelSelectionDialog)

    excel_sheet = None

    if excel_file is not None:
        xls = xlrd.open_workbook(excel_file, on_demand=True)
        sheet_names = xls.sheet_names()
        ui.sheets_list.addItems(sheet_names)

        if len(sheet_names) > 0:
            excel_sheet = 0

    def accepted():
        if len(ui.sheets_list.selectedIndexes()):
            excel_sheet = ui.sheets_list.selectedIndexes()[0].row()
        print('Accepted: self.excel_sheet: ', excel_sheet)

        app.quit()

    def rejected():
        print('Rejected: self.excel_sheet: ', excel_sheet)
        app.quit()

    # make connections
    ui.acceptButton.clicked.connect(accepted)
    ui.cancelButton.clicked.connect(rejected)
    ui.sheets_list.doubleClicked.connect(accepted)

    # show
    ExcelSelectionDialog.show()

    return excel_sheet


if __name__ == "__main__":
    excel_file = '/home/santi/Documentos/GitHub/GridCal/Grids_and_profiles/profiles/Total_profiles_1W_1H.xlsx'

    app = QtWidgets.QApplication(sys.argv)
    window = ExcelDialog(excel_file)
    window.show()
    # excel_dialogue_show(excel_file)
    sys.exit(app.exec_())

