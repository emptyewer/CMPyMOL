# -*- coding: utf-8 -*-

# The MIT License (MIT)
#
# Copyright (c) 2015 Venkatramanan Krishnamani
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


import sys
from PyQt4 import QtCore, QtGui
from functions.mitlicense import License_Dialog

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

class About_Dialog(QtGui.QDialog):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        self.setupUi(self)

    def setupUi(self, Dialog):
        Dialog.setObjectName(_fromUtf8("Dialog"))
        Dialog.resize(329, 358)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Dialog.sizePolicy().hasHeightForWidth())
        Dialog.setSizePolicy(sizePolicy)
        Dialog.setMinimumSize(QtCore.QSize(329, 358))
        Dialog.setMaximumSize(QtCore.QSize(329, 358))
        self.horizontalLayout = QtGui.QHBoxLayout(Dialog)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.groupBox = QtGui.QGroupBox(Dialog)
        self.groupBox.setObjectName(_fromUtf8("groupBox"))
        self.verticalLayout = QtGui.QVBoxLayout(self.groupBox)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        spacerItem = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem)
        self.label_3 = QtGui.QLabel(self.groupBox)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.verticalLayout.addWidget(self.label_3)
        self.label_2 = QtGui.QLabel(self.groupBox)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_2.sizePolicy().hasHeightForWidth())
        self.label_2.setSizePolicy(sizePolicy)
        self.label_2.setMinimumSize(QtCore.QSize(256, 100))
        self.label_2.setMaximumSize(QtCore.QSize(256, 100))
        self.label_2.setText(_fromUtf8(""))
        if sys.platform == 'darwin':
            self.label_2.setPixmap(QtGui.QPixmap(_fromUtf8("icon/Icon_256x256.png")))
        elif sys.platform == 'linux2' or sys.platform == 'linux':
            self.label_2.setPixmap(QtGui.QPixmap(_fromUtf8("icon/Icon_256x256.png")))
        elif sys.platform == 'win32':
            self.label_2.setPixmap(QtGui.QPixmap(_fromUtf8("icon/Icon_256x256.png")))
        self.label_2.setAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignTop)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.verticalLayout.addWidget(self.label_2, QtCore.Qt.AlignHCenter)
        self.label = QtGui.QLabel(self.groupBox)
        self.label.setObjectName(_fromUtf8("label"))
        self.verticalLayout.addWidget(self.label)
        self.textBrowser = QtGui.QTextBrowser(self.groupBox)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.textBrowser.sizePolicy().hasHeightForWidth())
        self.textBrowser.setSizePolicy(sizePolicy)
        self.textBrowser.setMaximumSize(QtCore.QSize(256, 80))
        self.textBrowser.setObjectName(_fromUtf8("textBrowser"))
        self.verticalLayout.addWidget(self.textBrowser)
        self.pushButton_2 = QtGui.QPushButton(self.groupBox)
        self.pushButton_2.setObjectName(_fromUtf8("pushButton_2"))
        self.verticalLayout.addWidget(self.pushButton_2)
        self.horizontalLayout.addWidget(self.groupBox)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)


    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(_translate("Dialog", "About", None))
        self.label_3.setText(_translate("Dialog", "<html><head/><body><p align=\"center\"><span style=\" "
                                                  "font-size:14pt; "
                                                  "font-weight:600;\">CMPyMOL 2.0</span></p></body></html>", None))
        self.label.setText(_translate("Dialog", "<html><head/><body><p align=\"center\">(c) 2016 Venkatramanan "
                                                "Krishnamani</p><p align=\"center\">University of Iowa</p></body></html>", None))
        self.textBrowser.setHtml(_translate("Dialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'.Helvetica Neue DeskInterface\'; font-size:13pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Liberation Sans\'; font-size:11pt;\">This software is licensed under the</span><span style=\" font-family:\'Liberation Sans\'; font-size:11pt; font-weight:600;\"> MIT License. </span><span style=\" font-family:\'Liberation Sans\'; font-size:11pt;\">It is a free software license that permits resue within proprietary software provided all copies of this software include a copy of the </span><span style=\" font-family:\'Liberation Sans\'; font-size:11pt; font-weight:600;\">MIT License</span><span style=\" font-family:\'Liberation Sans\'; font-size:11pt;\"> terms and copyright notice.</span></p></body></html>", None))
        self.pushButton_2.setText(_translate("Dialog", "STRIDE Copyright Notice", None))
        self.pushButton_2.clicked.connect(self.showLicense)

    def showLicense(self):
        lic_dialog = License_Dialog()
        lic_dialog.exec_()

