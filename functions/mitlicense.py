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

class License_Dialog(QtGui.QDialog):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        self.setupUi(self)

    def setupUi(self, Dialog):
        Dialog.setObjectName(_fromUtf8("Dialog"))
        Dialog.resize(713, 451)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Dialog.sizePolicy().hasHeightForWidth())
        Dialog.setSizePolicy(sizePolicy)
        Dialog.setMinimumSize(QtCore.QSize(713, 451))
        Dialog.setMaximumSize(QtCore.QSize(713, 451))
        self.horizontalLayout = QtGui.QHBoxLayout(Dialog)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.textBrowser = QtGui.QTextBrowser(Dialog)
        self.textBrowser.setObjectName(_fromUtf8("textBrowser"))
        self.horizontalLayout.addWidget(self.textBrowser)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(_translate("License", "License", None))
        self.textBrowser.setHtml(_translate("Dialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Liberation Sans\'; font-size:11pt; font-weight:400; font-style:normal;\">\n"
"<p>STRIDE COPYRIGHT NOTICE</p>"
"<p>All rights reserved, whether the whole or part of the program is concerned. Permission to use, copy,"
"and modify this software and its documentation is granted for academic use, provided that:</p>"
"<p> i. this copyright notice appears in all copies of the software and  related documentation;</p>"
"<p> ii. the reference given below (Frishman and Argos, 1995) must be  cited in any publication of scientific results based in part or  completely on the use of the program;</p>"
"<p> iii. bugs will be reported to the authors.</p>"
"<p> The use of the software in commercial activities is not allowed  without a prior written commercial license agreement.</p>"
"<p> WARNING: STRIDE is provided &quot;as-is&quot; and without warranty of any kind,  express, implied or otherwise,"
"including without limitation any  warranty of merchantability or fitness for a particular purpose. In no  event will "
"the authors be liable for any special, incidental, indirect  or consequential damages of any kind, or any damages "
"whatsoever  resulting from loss of data or profits, whether or not advised of the  possibility of damage, and on"
" any theory of liability, arising out of  or in connection with the use or performance of this software.</p>"
"<p> For calculation of the residue solvent accessible area the program NSC  [3,4] is used and was kindly provided by"
"Dr. F.Eisenhaber  (EISENHABER@EMBL-HEIDELBERG.DE). Please direct to him all questions "
"concerning specifically accessibility calculations.</p>"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>", None))

