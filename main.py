import numpy as np
import pyqtgraph as pg
from PyQt5 import uic
from PyQt5.QtWidgets import QApplication
from pyqtgraph import PlotWidget

import os
import solving as s

def button_clicked():
    с = 1.35
    l = 5
    k_thermal_cond=0.065
    alpha = 0.008
    T=150
    beta = 0.1
    R = 3
    u_0 = 0
    p = 60
    a = 0.1*R
    i = k = 1
    I = K = 1000

    if (form.lineEdit.text()!=''):
        l = (float)(form.lineEdit.text())
    if (form.lineEdit_2.text() != ''):
        R = (float)(form.lineEdit_2.text())
    if (form.lineEdit_3.text() != ''):
        beta = (float)(form.lineEdit_3.text())
    if (form.lineEdit_4.text() != ''):
        k_thermal_cond = (float)(form.lineEdit_4.text())
    if (form.lineEdit_5.text() != ''):
        c = (float)(form.lineEdit_5.text())
    if (form.lineEdit_6.text() != ''):
        alpha = (float)(form.lineEdit_6.text())
    if (form.lineEdit_7.text() != ''):
        T = (float)(form.lineEdit_7.text())
    if (form.lineEdit_8.text() != ''):
        u_0 = (float)(form.lineEdit_8.text())
    if (form.lineEdit_11.text() != ''):
        i = (float)(form.lineEdit_11.text())
    if (form.lineEdit_12.text() != ''):
        k = (float)(form.lineEdit_12.text())
    if (form.lineEdit_9.text() != ''):
        I = (float)(form.lineEdit_9.text())
    if (form.lineEdit_10.text() != ''):
        K = (float)(form.lineEdit_10.text())
    if (form.lineEdit_13.text() != ''):
        p = (float)(form.lineEdit_13.text())
    if (form.lineEdit_14.text() != ''):
        a = (float)(form.lineEdit_14.text())*R
    x = np.arange(1000)
    y = np.random.normal(size=(3, 1000))
    plotWidget = pg.plot(title="График для фиксированного i="+(str)(i))
    plotWidget2 = pg.plot(title="График для фиксированного k="+(str)(k))
    plotWidget.plot(x, y[i])  ## setting pen=(i,3) automaticaly creates three different-colored pens
    plotWidget2.plot(x, y[i])  ## setting pen=(i,3) automaticaly creates three different-colored pens

def button_clicked_2():
    с = 1.35
    l = 5
    k_thermal_cond = 0.065
    alpha = 0.008
    T = 150
    beta = 0.1
    R = 3
    u_0 = 0
    p = 60
    a = 0.1 * R
    i = k = 1
    I = K = 1000

    if (form.lineEdit.text() != ''):
        l = (float)(form.lineEdit.text())
    if (form.lineEdit_2.text() != ''):
        R = (float)(form.lineEdit_2.text())
    if (form.lineEdit_3.text() != ''):
        beta = (float)(form.lineEdit_3.text())
    if (form.lineEdit_4.text() != ''):
        k_thermal_cond = (float)(form.lineEdit_4.text())
    if (form.lineEdit_5.text() != ''):
        c = (float)(form.lineEdit_5.text())
    if (form.lineEdit_6.text() != ''):
        alpha = (float)(form.lineEdit_6.text())
    if (form.lineEdit_7.text() != ''):
        T = (float)(form.lineEdit_7.text())
    if (form.lineEdit_8.text() != ''):
        u_0 = (float)(form.lineEdit_8.text())
    if (form.lineEdit_11.text() != ''):
        i = (float)(form.lineEdit_11.text())
    if (form.lineEdit_12.text() != ''):
        k = (float)(form.lineEdit_12.text())
    if (form.lineEdit_9.text() != ''):
        I = (float)(form.lineEdit_9.text())
    if (form.lineEdit_10.text() != ''):
        K = (float)(form.lineEdit_10.text())
    if (form.lineEdit_13.text() != ''):
        p = (float)(form.lineEdit_13.text())
    if (form.lineEdit_14.text() != ''):
        a = (float)(form.lineEdit_14.text()) * R
    x = np.arange(1000)
    y = np.random.normal(size=(3, 1000))
    plotWidget = pg.plot(title="График для фиксированного i=" + (str)(i))
    plotWidget2 = pg.plot(title="График для фиксированного k=" + (str)(k))
    plotWidget.plot(x, y[i])  ## setting pen=(i,3) automaticaly creates three different-colored pens
    plotWidget2.plot(x, y[i])  ## setting pen=(i,3) automaticaly creates three different-colored pens

Form, Window = uic.loadUiType('D:\PythonProjects\\nmmp\interface.ui')
app = QApplication([])
window = Window()
form = Form()
form.setupUi(window)
form.pushButton.clicked.connect(button_clicked)
form.pushButton_2.clicked.connect(button_clicked_2)
window.show()
app.exec_()

