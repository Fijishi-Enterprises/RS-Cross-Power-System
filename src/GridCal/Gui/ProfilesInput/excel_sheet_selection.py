# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'excel_sheet_selection.ui',
# licensing of 'excel_sheet_selection.ui' applies.
#
# Created: Tue May 28 20:57:33 2019
#      by: pyside2-uic  running on PySide2 5.12.3
#
# WARNING! All changes made in this file will be lost!

from PySide2 import QtCore, QtGui, QtWidgets


class Ui_ExcelSelectionDialog(object):

    def setupUi(self, ExcelSelectionDialog):

        ExcelSelectionDialog.setObjectName("ExcelSelectionDialog")
        ExcelSelectionDialog.resize(272, 229)
        ExcelSelectionDialog.setMaximumSize(QtCore.QSize(272, 229))
        ExcelSelectionDialog.setModal(True)
        self.verticalLayout = QtWidgets.QVBoxLayout(ExcelSelectionDialog)
        self.verticalLayout.setContentsMargins(1, 1, 1, 1)
        self.verticalLayout.setObjectName("verticalLayout")
        self.sheets_list = QtWidgets.QListWidget(ExcelSelectionDialog)
        self.sheets_list.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.sheets_list.setObjectName("sheets_list")
        self.verticalLayout.addWidget(self.sheets_list)
        self.frame = QtWidgets.QFrame(ExcelSelectionDialog)
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

        self.retranslateUi(ExcelSelectionDialog)
        QtCore.QMetaObject.connectSlotsByName(ExcelSelectionDialog)

    def retranslateUi(self, ExcelSelectionDialog):
        ExcelSelectionDialog.setWindowTitle(QtWidgets.QApplication.translate("ExcelSelectionDialog", "Excel sheet selection", None, -1))
        self.cancelButton.setText(QtWidgets.QApplication.translate("ExcelSelectionDialog", "Cancel", None, -1))
        self.acceptButton.setText(QtWidgets.QApplication.translate("ExcelSelectionDialog", "Accept", None, -1))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    ExcelSelectionDialog = QtWidgets.QDialog()
    ui = Ui_ExcelSelectionDialog()
    ui.setupUi(ExcelSelectionDialog)

    excel_sheet = None

    ui.sheets_list.addItems(['A', 'B', 'C'])

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
    sys.exit(app.exec_())

