# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'button group.ui'
#
# Created by: PyQt4 UI code generator 4.11.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName(_fromUtf8("Form"))
        Form.resize(874, 212)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Form.sizePolicy().hasHeightForWidth())
        Form.setSizePolicy(sizePolicy)
        self.horizontalLayout = QtGui.QHBoxLayout(Form)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.beamPlotXAxisCombo = QtGui.QListWidget(Form)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.beamPlotXAxisCombo.sizePolicy().hasHeightForWidth())
        self.beamPlotXAxisCombo.setSizePolicy(sizePolicy)
        self.beamPlotXAxisCombo.setMinimumSize(QtCore.QSize(0, 100))
        self.beamPlotXAxisCombo.setObjectName(_fromUtf8("beamPlotXAxisCombo"))
        self.horizontalLayout.addWidget(self.beamPlotXAxisCombo)
        self.beamPlotYAxisCombo = QtGui.QListWidget(Form)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.beamPlotYAxisCombo.sizePolicy().hasHeightForWidth())
        self.beamPlotYAxisCombo.setSizePolicy(sizePolicy)
        self.beamPlotYAxisCombo.setMinimumSize(QtCore.QSize(0, 100))
        self.beamPlotYAxisCombo.setObjectName(_fromUtf8("beamPlotYAxisCombo"))
        self.horizontalLayout.addWidget(self.beamPlotYAxisCombo)
        self.groupBox = QtGui.QGroupBox(Form)
        self.groupBox.setObjectName(_fromUtf8("groupBox"))
        self.formLayout = QtGui.QFormLayout(self.groupBox)
        self.formLayout.setFieldGrowthPolicy(QtGui.QFormLayout.AllNonFixedFieldsGrow)
        self.formLayout.setObjectName(_fromUtf8("formLayout"))
        self.chooseSliceWidth = QtGui.QRadioButton(self.groupBox)
        self.chooseSliceWidth.setObjectName(_fromUtf8("chooseSliceWidth"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.LabelRole, self.chooseSliceWidth)
        self.sliceWidth = QtGui.QDoubleSpinBox(self.groupBox)
        self.sliceWidth.setObjectName(_fromUtf8("sliceWidth"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.FieldRole, self.sliceWidth)
        self.chooseSliceNumber = QtGui.QRadioButton(self.groupBox)
        self.chooseSliceNumber.setObjectName(_fromUtf8("chooseSliceNumber"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.LabelRole, self.chooseSliceNumber)
        self.sliceNumber = QtGui.QDoubleSpinBox(self.groupBox)
        self.sliceNumber.setObjectName(_fromUtf8("sliceNumber"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.FieldRole, self.sliceNumber)
        self.horizontalLayout.addWidget(self.groupBox)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(_translate("Form", "Form", None))
        self.groupBox.setTitle(_translate("Form", "Slices", None))
        self.chooseSliceWidth.setText(_translate("Form", "Slice Width", None))
        self.chooseSliceNumber.setText(_translate("Form", "Number of Slices", None))

