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


def get_a_krank(k, c, alpha, I, l, K, T):
    ht = T/K
    hz = l/I
    gamma = ht/hz/hz
    a = np.ones(I)*(-k*gamma/2/c)
    a[-1] *= 2
    return a

def get_b_implicit(k, c, alpha, I, l):
    b = np.arange(0, I, dtype=float)
    h_z=l/I
    for i in range(1, I):
        b[i]=-k/(c*h_z**2)
    b[0]=-k/(alpha*h_z**2)
    return b


def get_b_krank(k, c, alpha, I, l, K, T):
    ht = T / K
    hz = l / I
    gamma = ht / hz / hz
    b = np.ones(I) * (-k * gamma / 2 / c)
    b[0] *= 2
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
    ht = T / K
    hz = l / I
    gamma = ht / hz / hz
    cm = np.ones(I + 1) * (1 + k * gamma / c)
    cm[0] += alpha * gamma * hz / c
    cm[-1] += alpha * gamma * hz / c
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
    ht = T/K
    hz = l/I
    gamma = ht / hz / hz
    phi = get_phi(c, R, beta, p, a, I, l)
    f[0] = (1 - k*gamma/c - alpha*gamma*hz/c)*u_k_minus_1[0] + k*gamma/c*u_k_minus_1[1] + 2*alpha*gamma*hz/c*u_0 + ht*phi[0]
    f[1:I] = k*gamma/2/c*u_k_minus_1[:I-1]+(1-k*gamma/c)*u_k_minus_1[1:I]+k*gamma/2/c*u_k_minus_1[2:]+ht*phi[1:I]
    f[I] = (1 - k*gamma/c - alpha*gamma*hz/c)*u_k_minus_1[I] + k*gamma/c*u_k_minus_1[I-1] + 2*alpha*gamma*hz/c*u_0 + ht*phi[I]
    return f

def get_phi(c, R, beta, p, a, I, l):
    integral = 0
    z = np.arange(0, l, l/(I+1))
    for i in range(1, 1001):
        r = i/1000*R
        r_1 = (i-1)/1000*R
        integral+=p*np.exp(-(r/a)**2)/2/a**2*r+p*np.exp(-(r_1/a)**2)/2/a**2*r_1
    integral*=1/2000
    return beta*2*np.exp(-beta*z)/c/R**2


def get_phi_krank(c, R, beta, p, a, I, l):
    integral = 0
    z = np.arange(0, l, l/(I+1))
    for i in range(1, 1001):
        r = i/1000*R
        r_1 = (i-1)/1000*R
        integral+=p*np.exp(-(r/a)**2)/2/a**2*r+p*np.exp(-(r_1/a)**2)/2/a**2*r_1
    integral*=1/2000
    return 2*beta*np.exp(-beta*z)/c/R**2

def compute_implicit(I, k_thermal_cond, alpha, c, i, k, l, T, K, u_0, R, beta, p, a):
    am = get_a_implicit(k, c, alpha, I, l)
    bm = get_b_implicit(k, c, alpha, I, l)
    cm = get_c_implicit(k, c, alpha, I, l, T, K)
    fm = get_f_implicit(u_0, c, R, beta, p, a, np.full(I, u_0), alpha, I, T, K, l)
    u = np.zeros((K, I+1))

    #u[0, :] = s.thomas_method(am, bm, cm, fm)
    u[0, :] = u_0
    for i in range(1, K):
        fm = get_f_implicit(u_0, c, R, beta, p, a, u[i - 1, :], alpha, I, T, K, l)
        u[i, :]=s.thomas_method(am, bm, cm, fm)
    return u[i, :], u[:, k]


def compute_krank(I, k_thermal_cond, alpha, c, i, k, l, T, K, u_0, R, beta, p, a):
    am = get_a_krank(k, c, alpha, I, l, K, T)
    bm = get_b_krank(k, c, alpha, I, l, K, T)
    cm = get_c_krank(k, c, alpha, I, l, T, K)
    fm = get_f_krank(u_0, c, R, beta, p, a, np.full(I+1, u_0), alpha, I, T, K, l, k)
    u = np.zeros((K+1, I+1))
    u[0, :] = u_0
    for n in range(1, K+1):
        fm = get_f_krank(u_0, c, R, beta, p, a, u[n - 1, :], alpha, I, T, K, l, k)
        u[n, :] = s.thomas_method(am, bm, cm, fm)
    return u[:, i], u[k, :]


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
    ts = np.arange(0, T, T / (K+1))
    plotWidget = pg.plot(title="График для фиксированного i="+(str)(i))
    plotWidget2 = pg.plot(title="График для фиксированного k="+(str)(k))
    plotWidget.plot(ts, u_i)
    plotWidget2.plot(zs, u_k)


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

    if form.lineEdit.text() != '':
        l = float(form.lineEdit.text())
    if form.lineEdit_2.text() != '':
        R = float(form.lineEdit_2.text())
    if form.lineEdit_3.text() != '':
        beta = float(form.lineEdit_3.text())
    if form.lineEdit_4.text() != '':
        k_thermal_cond = float(form.lineEdit_4.text())
    if form.lineEdit_5.text() != '':
        c = float(form.lineEdit_5.text())
    if form.lineEdit_6.text() != '':
        alpha = float(form.lineEdit_6.text())
    if form.lineEdit_7.text() != '':
        T = float(form.lineEdit_7.text())
    if form.lineEdit_8.text() != '':
        u_0 = float(form.lineEdit_8.text())
    if form.lineEdit_11.text() != '':
        i = int(form.lineEdit_11.text())
    if form.lineEdit_12.text() != '':
        k = int(form.lineEdit_12.text())
    if form.lineEdit_9.text() != '':
        I = int(form.lineEdit_9.text())
    if form.lineEdit_10.text() != '':
        K = int(form.lineEdit_10.text())
    if form.lineEdit_13.text() != '':
        p = float(form.lineEdit_13.text())
    if form.lineEdit_14.text() != '':
        a = float(form.lineEdit_14.text()) * R
    u_i, u_k = compute_krank(I, k_thermal_cond, alpha, c, i, k, l, T, K, u_0, R, beta, p, a)
    zs = np.arange(0, l, l / (I+1))
    ts = np.arange(0, T, T / (K+1))
    plot_widget = pg.plot(title="График для фиксированного i=" + str(i))
    plot_widget2 = pg.plot(title="График для фиксированного k=" + str(k))
    plot_widget.plot(ts, u_i)
    plot_widget2.plot(zs, u_k)


Form, Window = uic.loadUiType('interface.ui')
app = QApplication([])
window = Window()
form = Form()
form.setupUi(window)
form.pushButton.clicked.connect(button_clicked)
form.pushButton_2.clicked.connect(button_clicked_2)
window.show()
app.exec_()

