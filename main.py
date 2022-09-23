import numpy as np
import pyqtgraph as pg
from PyQt5 import uic
from PyQt5.QtWidgets import QApplication
from pyqtgraph import PlotWidget

import os
import solving as s

def get_a_implicit(k, c, alpha, I, l):
    a = np.arange(0, I, dtype=float)
    h_z=l/I
    for i in range(0, I-1):
        a[i]=-k/(c*h_z**2)
    a[I-1]=-k/(alpha*h_z**2)
    return a

def get_a_krank(k, c, alpha, I, l):
    a = np.arange(0, I, dtype=float)
    h_z=l/I
    for i in range(0, I-1):
        a[i]=-k/(2*c*h_z)
    a[I-1]=-k/(alpha*h_z)
    return a

def get_b_implicit(k, c, alpha, I, l):
    b = np.arange(0, I, dtype=float)
    h_z=l/I
    for i in range(1, I):
        b[i]=-k/(c*h_z**2)
    b[0]=-k/(alpha*h_z**2)
    return b

def get_b_krank(k, c, alpha, I, l):
    b = np.arange(0, I, dtype=float)
    h_z=l/I
    for i in range(1, I):
        b[i]=-k/(2*c*h_z)
    b[0]=-k/(alpha*h_z)
    return b

def get_c_implicit(k, c, alpha, I, l, T, K):
    cm = np.arange(0, I+1, dtype=float)
    h_z=l/I
    h_t=T/K
    cm[0]=1-k/(alpha*h_z)
    for i in range(1, I):
        cm[i]=1/h_t+2*k/(c*h_z**2)
    cm[I]=1+k/(alpha*h_z)
    return cm

def get_c_krank(k, c, alpha, I, l, T, K):
    cm = np.arange(0, I+1, dtype=float)
    h_z=l/I
    h_t=T/K
    cm[0]=1+k/alpha/h_z+c*h_z/2/alpha/h_t
    for i in range(1, I):
        cm[i]=1/h_t+1/(c*h_z**2)
    cm[I]=1+k/(alpha*h_z)+c*h_z/2/alpha/h_t
    return cm

def get_f_implicit(u_0, c, R, beta, p, a, u_k_minus_1, alpha, I, T, K, l):
    f = np.arange(0, I+1, dtype=float)
    h_t=T/K
    phi = get_phi(c, R, beta, p, a, I, l)
    f[0]=u_0
    for i in range(1, I):
        f[i]=phi[i]+u_k_minus_1[i]/h_t
    f[I]=u_0
    return f

def get_f_krank(u_0, c, R, beta, p, a, u_k_minus_1, alpha, I, T, K, l, k):
    f = np.arange(0, I+1, dtype=float)
    h_t=T/K
    h_z=l/I
    phi = get_phi(c, R, beta, p, a, I, l)
    f[0]=c*h_z/(2*alpha*h_t)*u_k_minus_1[0]+c*h_z/2/alpha*phi[0]+u_0
    for i in range(1, I):
        f[i]=phi[i]+k/(2*c*h_z**2)*u_k_minus_1[i+1]+k/(2*c*h_z**2)*u_k_minus_1[i-1]+(1/h_t-k/(c*h_z**2))*u_k_minus_1[i]
    f[I]=c*h_z/(2*alpha*h_t)*u_k_minus_1[I]+c*h_z/2/alpha*phi[I]+u_0
    return f

def get_phi(c, R, beta, p, a, I, l):
    integral = 0
    z = np.arange(0, l, l/(I+1))
    for i in range(1, 1001):
        r = i/1000*R
        r_1 = (i-1)/1000*R
        integral+=p*np.exp(-(r/a)**2)/2/a**2*r+p*np.exp(-(r_1/a)**2)/2/a**2*r_1
    integral*=1/2000
    return 2*np.exp(-beta*z)/c/R**2

def get_phi_krank(c, R, beta, p, a, I, l):
    integral = 0
    z = np.arange(0, l, l/(I+1))
    for i in range(1, 1001):
        r = i/1000*R
        r_1 = (i-1)/1000*R
        integral+=p*np.exp(-(r/a)**2)/2/a**2*r+p*np.exp(-(r_1/a)**2)/2/a**2*r_1
    integral*=1/2000
    return 2*np.exp(-beta*z)/c/R**2

def compute_implicit(I, k_thermal_cond, alpha, c, i, k, l, T, K, u_0, R, beta, p, a):
    am = get_a_implicit(k, c, alpha, I, l)
    bm = get_b_implicit(k, c, alpha, I, l)
    cm = get_c_implicit(k, c, alpha, I, l, T, K)
    fm = get_f_implicit(u_0, c, R, beta, p, a, np.full(I, u_0), alpha, I, T, K, l)
    u = np.zeros((K, I+1))
    u[0, :] = s.thomas_method(am, bm, cm, fm)
    for i in range(1, K):
        fm = get_f_implicit(u_0, c, R, beta, p, a, u[i - 1, :], alpha, I, T, K, l)
        u[i, :]=s.thomas_method(am, bm, cm, fm)
    return u[i, :], u[:, k]

def compute_krank(I, k_thermal_cond, alpha, c, i, k, l, T, K, u_0, R, beta, p, a):
    am = get_a_krank(k, c, alpha, I, l)
    bm = get_b_krank(k, c, alpha, I, l)
    cm = get_c_krank(k, c, alpha, I, l, T, K)
    fm = get_f_krank(u_0, c, R, beta, p, a, np.full(I+1, u_0), alpha, I, T, K, l, k)
    u = np.zeros((K, I+1))
    u[0, :] = s.thomas_method(am, bm, cm, fm)
    for i in range(1, K):
        fm = get_f_krank(u_0, c, R, beta, p, a, u[i - 1, :], alpha, I, T, K, l, k)
        u[i, :]=s.thomas_method(am, bm, cm, fm)
    return u[i, :], u[:, k]


def button_clicked():
    c = 1.35
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
        I = (int)(form.lineEdit_9.text())
    if (form.lineEdit_10.text() != ''):
        K = (int)(form.lineEdit_10.text())
    if (form.lineEdit_13.text() != ''):
        p = (float)(form.lineEdit_13.text())
    if (form.lineEdit_14.text() != ''):
        a = (float)(form.lineEdit_14.text())*R
    u_i, u_k = compute_implicit(I, k_thermal_cond, alpha, c, i, k, l, T, K, u_0,R, beta, p ,a)
    zs = np.arange(0, l, l/(I+1))
    ts = np.arange(0, T, T / (K))
    plotWidget = pg.plot(title="График для фиксированного i="+(str)(i))
    plotWidget2 = pg.plot(title="График для фиксированного k="+(str)(k))
    plotWidget.plot(zs, u_i)
    plotWidget2.plot(ts, u_k)

def button_clicked_2():
    c = 1.35
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
        I = (int)(form.lineEdit_9.text())
    if (form.lineEdit_10.text() != ''):
        K = (int)(form.lineEdit_10.text())
    if (form.lineEdit_13.text() != ''):
        p = (float)(form.lineEdit_13.text())
    if (form.lineEdit_14.text() != ''):
        a = (float)(form.lineEdit_14.text()) * R
    u_i, u_k = compute_krank(I, k_thermal_cond, alpha, c, i, k, l, T, K, u_0, R, beta, p, a)
    zs = np.arange(0, l, l / (I+1))
    ts = np.arange(0, T, T / (K))
    plotWidget = pg.plot(title="График для фиксированного i=" + (str)(i))
    plotWidget2 = pg.plot(title="График для фиксированного k=" + (str)(k))
    plotWidget.plot(zs, u_i)
    plotWidget2.plot(ts, u_k)





Form, Window = uic.loadUiType('D:\PythonProjects\\nmmp\interface.ui')
app = QApplication([])
window = Window()
form = Form()
form.setupUi(window)
form.pushButton.clicked.connect(button_clicked)
form.pushButton_2.clicked.connect(button_clicked_2)
window.show()
app.exec_()

