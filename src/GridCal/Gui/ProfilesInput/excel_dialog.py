
import os
import sys
import xlrd
from PySide2 import QtCore, QtGui, QtWidgets
# from PyQt5 import QtCore, QtGui, QtWidgets


class ExcelDialog(QtWidgets.QDialog):

    def __init__(self, excel_file=None, items=()):
        """
        Constructor
        :param excel_file: excel file to list
        :param items: items to show if the file is none
        """

        QtWidgets.QDialog.__init__(self)

        self.setObjectName("ExcelSelectionDialog")
        self.resize(272, 229)
        self.setMaximumSize(QtCore.QSize(272, 229))
        self.setModal(True)
        self.verticalLayout = QtWidgets.QVBoxLayout(self)
        self.verticalLayout.setContentsMargins(1, 1, 1, 1)
        self.verticalLayout.setObjectName("verticalLayout")
        self.sheets_list = QtWidgets.QListWidget(self)
        self.sheets_list.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.sheets_list.setObjectName("sheets_list")
        self.verticalLayout.addWidget(self.sheets_list)
        self.frame = QtWidgets.QFrame(self)
        self.frame.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame.setObjectName("frame")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.frame)
        self.horizontalLayout.setContentsMargins(1, 1, 1, 1)
        self.horizontalLayout.setObjectName("horizontalLayout")
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.cancelButton = QtWidgets.QPushButton(self.frame)
        self.cancelButton.setObjectName("cancelButton")
        self.horizontalLayout.addWidget(self.cancelButton)
        self.acceptButton = QtWidgets.QPushButton(self.frame)
        self.acceptButton.setObjectName("acceptButton")
        self.horizontalLayout.addWidget(self.acceptButton)
        self.verticalLayout.addWidget(self.frame)

        self.retranslateUi(self)
        QtCore.QMetaObject.connectSlotsByName(self)

        # click
        self.acceptButton.clicked.connect(self.accepted)
        self.cancelButton.clicked.connect(self.rejected)
        self.sheets_list.doubleClicked.connect(self.accepted)

        self.excel_sheet = None

        self.sheet_names = list()
        if excel_file is not None:
            if os.path.exists(excel_file):
                self.fill_from_file(excel_file=excel_file)
            else:
                self.sheets_list.addItems(items)
        else:
            self.sheets_list.addItems(items)

    def fill_from_file(self, excel_file):
        """

        :param excel_file:
        :return:
        """
        if excel_file is not None:
            xls = xlrd.open_workbook(excel_file, on_demand=True)
            self.sheet_names = xls.sheet_names()
            self.sheets_list.addItems(self.sheet_names)

            if len(self.sheet_names) > 0:
                self.excel_sheet = 0

    def accepted(self):
        """

        :return:
        """
        if len(self.sheets_list.selectedIndexes()):
            self.excel_sheet = self.sheets_list.selectedIndexes()[0].row()
        print('Accepted: self.excel_sheet: ', self.excel_sheet)

        self.close()

    def rejected(self):
        """

        :return:
        """
        print('Rejected: self.excel_sheet: ', self.excel_sheet)
        self.close()

    def retranslateUi(self, ExcelSelectionDialog):
        """

        :param ExcelSelectionDialog:
        :return:
        """
        ExcelSelectionDialog.setWindowTitle(QtWidgets.QApplication.translate("ExcelSelectionDialog", "Excel sheet selection", None, -1))
        self.cancelButton.setText(QtWidgets.QApplication.translate("ExcelSelectionDialog", "Cancel", None, -1))
        self.acceptButton.setText(QtWidgets.QApplication.translate("ExcelSelectionDialog", "Accept", None, -1))


if __name__ == "__main__":
    excel_file = None
    app = QtWidgets.QApplication(sys.argv)
    window = ExcelDialog(excel_file, items=['A', 'B', 'C'])
    window.show()
    sys.exit(app.exec_())

