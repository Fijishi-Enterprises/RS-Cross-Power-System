# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'gui.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(949, 537)
        self.verticalLayout = QtWidgets.QVBoxLayout(Dialog)
        self.verticalLayout.setContentsMargins(1, 1, 1, 1)
        self.verticalLayout.setObjectName("verticalLayout")
        self.tabWidget = QtWidgets.QTabWidget(Dialog)
        self.tabWidget.setObjectName("tabWidget")
        self.tab_2 = QtWidgets.QWidget()
        self.tab_2.setObjectName("tab_2")
        self.verticalLayout_6 = QtWidgets.QVBoxLayout(self.tab_2)
        self.verticalLayout_6.setObjectName("verticalLayout_6")
        self.splitter_3 = QtWidgets.QSplitter(self.tab_2)
        self.splitter_3.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_3.setObjectName("splitter_3")
        self.frame_8 = QtWidgets.QFrame(self.splitter_3)
        self.frame_8.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_8.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_8.setObjectName("frame_8")
        self.verticalLayout_5 = QtWidgets.QVBoxLayout(self.frame_8)
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.tabWidget_2 = QtWidgets.QTabWidget(self.frame_8)
        self.tabWidget_2.setObjectName("tabWidget_2")
        self.tab_4 = QtWidgets.QWidget()
        self.tab_4.setObjectName("tab_4")
        self.verticalLayout_8 = QtWidgets.QVBoxLayout(self.tab_4)
        self.verticalLayout_8.setObjectName("verticalLayout_8")
        self.tower_tableView = QtWidgets.QTableView(self.tab_4)
        self.tower_tableView.setObjectName("tower_tableView")
        self.verticalLayout_8.addWidget(self.tower_tableView)
        self.frame = QtWidgets.QFrame(self.tab_4)
        self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame.setObjectName("frame")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.frame)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.add_to_tower_pushButton = QtWidgets.QPushButton(self.frame)
        self.add_to_tower_pushButton.setObjectName("add_to_tower_pushButton")
        self.horizontalLayout_2.addWidget(self.add_to_tower_pushButton)
        self.delete_from_tower_pushButton = QtWidgets.QPushButton(self.frame)
        self.delete_from_tower_pushButton.setObjectName("delete_from_tower_pushButton")
        self.horizontalLayout_2.addWidget(self.delete_from_tower_pushButton)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem)
        self.compute_pushButton = QtWidgets.QPushButton(self.frame)
        self.compute_pushButton.setObjectName("compute_pushButton")
        self.horizontalLayout_2.addWidget(self.compute_pushButton)
        self.verticalLayout_8.addWidget(self.frame)
        self.tabWidget_2.addTab(self.tab_4, "")
        self.tab_5 = QtWidgets.QWidget()
        self.tab_5.setObjectName("tab_5")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.tab_5)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.wires_tableView = QtWidgets.QTableView(self.tab_5)
        self.wires_tableView.setObjectName("wires_tableView")
        self.verticalLayout_4.addWidget(self.wires_tableView)
        self.frame_4 = QtWidgets.QFrame(self.tab_5)
        self.frame_4.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_4.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_4.setObjectName("frame_4")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.frame_4)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.add_wire_pushButton = QtWidgets.QPushButton(self.frame_4)
        self.add_wire_pushButton.setObjectName("add_wire_pushButton")
        self.horizontalLayout.addWidget(self.add_wire_pushButton)
        self.delete_wire_pushButton = QtWidgets.QPushButton(self.frame_4)
        self.delete_wire_pushButton.setObjectName("delete_wire_pushButton")
        self.horizontalLayout.addWidget(self.delete_wire_pushButton)
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem1)
        self.verticalLayout_4.addWidget(self.frame_4)
        self.tabWidget_2.addTab(self.tab_5, "")
        self.verticalLayout_5.addWidget(self.tabWidget_2)
        self.PlotFrame = QtWidgets.QFrame(self.splitter_3)
        self.PlotFrame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.PlotFrame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.PlotFrame.setObjectName("PlotFrame")
        self.verticalLayout_7 = QtWidgets.QVBoxLayout(self.PlotFrame)
        self.verticalLayout_7.setContentsMargins(9, 9, 9, 9)
        self.verticalLayout_7.setObjectName("verticalLayout_7")
        self.label_4 = QtWidgets.QLabel(self.PlotFrame)
        self.label_4.setObjectName("label_4")
        self.verticalLayout_7.addWidget(self.label_4)
        self.plotwidget = MatplotlibWidget(self.PlotFrame)
        self.plotwidget.setObjectName("plotwidget")
        self.verticalLayout_7.addWidget(self.plotwidget)
        self.verticalLayout_6.addWidget(self.splitter_3)
        self.tabWidget.addTab(self.tab_2, "")
        self.tab = QtWidgets.QWidget()
        self.tab.setObjectName("tab")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.tab)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.label = QtWidgets.QLabel(self.tab)
        self.label.setObjectName("label")
        self.verticalLayout_3.addWidget(self.label)
        self.z_tableView_abcn = QtWidgets.QTableView(self.tab)
        self.z_tableView_abcn.setObjectName("z_tableView_abcn")
        self.verticalLayout_3.addWidget(self.z_tableView_abcn)
        self.label_6 = QtWidgets.QLabel(self.tab)
        self.label_6.setObjectName("label_6")
        self.verticalLayout_3.addWidget(self.label_6)
        self.z_tableView_abc = QtWidgets.QTableView(self.tab)
        self.z_tableView_abc.setObjectName("z_tableView_abc")
        self.verticalLayout_3.addWidget(self.z_tableView_abc)
        self.label_7 = QtWidgets.QLabel(self.tab)
        self.label_7.setObjectName("label_7")
        self.verticalLayout_3.addWidget(self.label_7)
        self.z_tableView_seq = QtWidgets.QTableView(self.tab)
        self.z_tableView_seq.setObjectName("z_tableView_seq")
        self.verticalLayout_3.addWidget(self.z_tableView_seq)
        self.tabWidget.addTab(self.tab, "")
        self.tab_3 = QtWidgets.QWidget()
        self.tab_3.setObjectName("tab_3")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.tab_3)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.label_2 = QtWidgets.QLabel(self.tab_3)
        self.label_2.setObjectName("label_2")
        self.verticalLayout_2.addWidget(self.label_2)
        self.y_tableView_abcn = QtWidgets.QTableView(self.tab_3)
        self.y_tableView_abcn.setObjectName("y_tableView_abcn")
        self.verticalLayout_2.addWidget(self.y_tableView_abcn)
        self.label_3 = QtWidgets.QLabel(self.tab_3)
        self.label_3.setObjectName("label_3")
        self.verticalLayout_2.addWidget(self.label_3)
        self.y_tableView_abc = QtWidgets.QTableView(self.tab_3)
        self.y_tableView_abc.setObjectName("y_tableView_abc")
        self.verticalLayout_2.addWidget(self.y_tableView_abc)
        self.label_5 = QtWidgets.QLabel(self.tab_3)
        self.label_5.setObjectName("label_5")
        self.verticalLayout_2.addWidget(self.label_5)
        self.y_tableView_seq = QtWidgets.QTableView(self.tab_3)
        self.y_tableView_seq.setObjectName("y_tableView_seq")
        self.verticalLayout_2.addWidget(self.y_tableView_seq)
        self.tabWidget.addTab(self.tab_3, "")
        self.verticalLayout.addWidget(self.tabWidget)

        self.retranslateUi(Dialog)
        self.tabWidget.setCurrentIndex(0)
        self.tabWidget_2.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Tower creation"))
        self.add_to_tower_pushButton.setText(_translate("Dialog", "Add"))
        self.delete_from_tower_pushButton.setText(_translate("Dialog", "Delete"))
        self.compute_pushButton.setText(_translate("Dialog", "Compute"))
        self.tabWidget_2.setTabText(self.tabWidget_2.indexOf(self.tab_4), _translate("Dialog", "Wire composition (Tower)"))
        self.add_wire_pushButton.setText(_translate("Dialog", "Add"))
        self.delete_wire_pushButton.setText(_translate("Dialog", "Delete"))
        self.tabWidget_2.setTabText(self.tabWidget_2.indexOf(self.tab_5), _translate("Dialog", "Wires catalogue"))
        self.label_4.setText(_translate("Dialog", "Tower"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), _translate("Dialog", "Tower"))
        self.label.setText(_translate("Dialog", "   Z series (Ohm / km) fpr ABCN"))
        self.label_6.setText(_translate("Dialog", "   Z series (Ohm / km) for ABC"))
        self.label_7.setText(_translate("Dialog", "   Z series (Ohm / km) in sequence components"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("Dialog", "Z"))
        self.label_2.setText(_translate("Dialog", "   Y shunt (uS / km) for ABCN"))
        self.label_3.setText(_translate("Dialog", "   Y shunt (uS / km) for ABC"))
        self.label_5.setText(_translate("Dialog", "   Y shunt (uS / km) for the sequence components"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_3), _translate("Dialog", "Y"))

from .matplotlibwidget import MatplotlibWidget

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = QtWidgets.QDialog()
    ui = Ui_Dialog()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec_())

