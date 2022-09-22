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

def get_a_implicit(k, c, alpha, I, l):
    a = np.arange(0, I, dtype=float)
    h_z=l/I
    for i in range(0, I-1):
        a[i]=-k/(c*h_z**2)
    a[I-1]=-k/(alpha*h_z**2)
    return a

def get_b_implicit(k, c, alpha, I, l):
    b = np.arange(0, I, dtype=float)
    h_z=l/I
    for i in range(1, I):
        b[i]=-k/(c*h_z**2)
    b[0]=-k/(alpha*h_z**2)
    return b

def get_c_implicit(k, c, alpha, I, l, T, K):
    c = np.arange(0, I+1, dtype=float)
    h_z=l/I
    h_t=T/K
    c[0]=1-k/(alpha*h_z)
    for i in range(1, I):
        c[i]=1/h_t+2*k/(c*h_z**2)
    c[I]=1+k/(alpha*h_z)
    return c

def get_f_implicit(u_0, c, R, beta, p, a, u_k_minus_1):
    f = np.arange(0, I+1, dtype=float)
    h_t=T/K
    f[0]=u_0
    phi = get_phi(c, R, beta, p, a, I, l)
    for i in range(1, I):
        f[i]=phi[i]+u_k_minus_1/h_t
    c[I]=1+k/(alpha*h_z)
    return c

def get_phi(c, R, beta, p, a, I, l):
    integral = 0
    z = np.arange(0, l, l/I)
    for i in range(1001):
        r = i/1000*R
        integral+=p*np.exp(-(r/a)**2)/2/a**2*r
    return 2*np.exp(-beta*z)/c/R**2

def compute(I, k_thermal_cond, alpha, c, i, k, l, T, K, u_0):
    a = get_a_implicit(k, c, alpha, I, l)
    b = get_b_implicit(k, c, alpha, I, l)
    c = get_c_implicit(k, c, alpha, I, l, T, K)
    f = get_f_implicit(u_0, c, R, beta, p, a, u_k_minus_1)


Form, Window = uic.loadUiType('D:\PythonProjects\\nmmp\interface.ui')
app = QApplication([])
window = Window()
form = Form()
form.setupUi(window)
form.pushButton.clicked.connect(button_clicked)
form.pushButton_2.clicked.connect(button_clicked_2)
window.show()
app.exec_()

