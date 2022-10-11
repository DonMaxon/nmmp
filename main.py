import numpy as np
import pyqtgraph as pg
from PyQt5 import uic
from PyQt5.QtWidgets import QApplication
from pyqtgraph import PlotWidget
from matplotlib import pyplot as plt

import os
import solving as s


def get_a_implicit(k, c, alpha, I, l):
    h_z = l / I
    a = -k/(c*h_z**2)
    a1 = -k/(alpha*h_z)
    return a, a1


def get_a_krank(k, c, I, l, K, T):
    ht = T/K
    hz = l/I
    gamma = ht/hz/hz
    a = -k*gamma/2/c
    a1 = 2*a
    return a, a1


def get_c_implicit(k, c, alpha, I, l, T, K):
    h_z=l/I
    h_t=T/K
    cm = 1/h_t+2*k/(c*h_z**2)
    cm1=1+k/(alpha*h_z)
    return cm, cm1


def get_c_krank(k, c, alpha, I, l, T, K):
    ht = T / K
    hz = l / I
    gamma = ht / hz / hz
    cm = 1 + k * gamma / c
    cm1 = cm + alpha * gamma * hz / c
    return cm, cm1


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
    ht = T / K
    hz = l / I
    gamma = ht / hz / hz
    phi = get_phi(c, R, beta, p, a, I, l)
    f[0] = (1 - k * gamma / c - alpha * gamma * hz / c) * u_k_minus_1[0] + k * gamma / c * u_k_minus_1[1] + 2 * alpha * gamma * hz / c * u_0 + ht * phi[0]
    f[1:I] = k * gamma / 2 / c * u_k_minus_1[:I - 1] + (1 - k * gamma / c) * u_k_minus_1[1:I] + k * gamma / 2 / c * u_k_minus_1[2:] + ht * phi[1:I]
    f[I] = (1 - k * gamma / c - alpha * gamma * hz / c) * u_k_minus_1[I] + k * gamma / c * u_k_minus_1[I - 1] + 2 * alpha * gamma * hz / c * u_0 + ht * phi[I]
    return f



def get_phi(c, R, beta, p, a, I, l):
    z = np.arange(0, l, l/(I+1))
    return beta*2*np.exp(-beta*z)/c/R**2*p/4*(1-np.exp(-R**2/a**2))


def compute_implicit(I, k_thermal_cond, alpha, c, i, k, l, T, K, u_0, R, beta, p, a):
    am, am1 = get_a_implicit(k_thermal_cond, c, alpha, I, l)
    cm, cm1 = get_c_implicit(k_thermal_cond, c, alpha, I, l, T, K)
    u = np.zeros((K+1, I+1))
    u[0, :] = u_0
    for j in range(1, K+1):
        fm = get_f_implicit(u_0, c, R, beta, p, a, u[j - 1, :], alpha, I, T, K, l)
        u[j, :]=s.thomas_method(am, cm, am1, cm1, fm)
    return u


def compute_krank(I, k_thermal_cond, alpha, c, l, T, K, u_0, R, beta, p, a):
    am, am1 = get_a_krank(k_thermal_cond, c, I, l, K, T)
    cm, cm1 = get_c_krank(k_thermal_cond, c, alpha, I, l, T, K)
    u = np.zeros((K+1, I+1))
    u[0, :] = u_0
    for j in range(1, K+1):
        fm = get_f_krank(u_0, c, R, beta, p, a, u[j - 1, :], alpha, I, T, K, l, k_thermal_cond)
        u[j, :] = s.thomas_method(am, cm, am1, cm1, fm)
    return u


def button_clicked():
    c = 1.35
    l = 5
    k_thermal_cond=0.065
    alpha = 0.008
    T = 150
    beta = 0.1
    R = 3
    u_0 = 0
    p = 60
    a = 0.1*R
    i = k = 1
    I = K = 1000
    mult_graph_k = False
    mult_graph_i = False
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
    else:
        mult_graph_i = True
    if form.lineEdit_12.text() != '':
        k = int(form.lineEdit_12.text())
    else:
        mult_graph_k = True
    if form.lineEdit_9.text() != '':
        I = int(form.lineEdit_9.text())
    if form.lineEdit_10.text() != '':
        K = int(form.lineEdit_10.text())
    if form.lineEdit_13.text() != '':
        p = float(form.lineEdit_13.text())
    if form.lineEdit_14.text() != '':
        a = float(form.lineEdit_14.text()) * R
    plot_widget = pg.plot(title="График для i")
    plot_widget.showGrid(x=True, y=True)
    plot_widget.setBackground('w')
    plot_widget2 = pg.plot(title="График для k")
    plot_widget2.setBackground('w')
    plot_widget2.showGrid(x=True, y=True)
    if not mult_graph_i and not mult_graph_k:
        u = compute_implicit(I, k_thermal_cond, alpha, c, i, k, l, T, K, u_0, R, beta, p, a)
        zs = np.arange(0, l, l / (I + 1))
        ts = np.arange(0, T, T / (K + 1))
        plot_widget.plot(ts, u[:, i], pen=pg.mkPen(color=(0, 0, 0)))
        plot_widget2.plot(zs, u[k, :], pen=pg.mkPen(color=(0, 0, 0)))
    if not mult_graph_i and mult_graph_k:
        u = compute_implicit(I, k_thermal_cond, alpha, c, i, k, l, T, K, u_0, R, beta, p, a)
        ts = np.arange(0, T, T / (K + 1))
        zs = np.arange(0, l, l / (I + 1))
        plot_widget.plot(ts, u[:, i], pen=pg.mkPen(color=(0, 0, 0)))
        plot_widget2.plot(zs, u[K - 1, :], pen=pg.mkPen(color=(0, 0, 0)))
        for j in range(4):
            plot_widget2.plot(zs, u[int(j * K / 4) + 1, :], pen=pg.mkPen(color=(0, 0, 0)))
    if mult_graph_i and not mult_graph_k:
        u = compute_implicit(I, k_thermal_cond, alpha, c, I - 1, k, l, T, K, u_0, R, beta, p, a)
        ts = np.arange(0, T, T / (K + 1))
        zs = np.arange(0, l, l / (I + 1))
        plot_widget.plot(ts, u[:, I - 1], pen=pg.mkPen(color=(0, 0, 0)))
        plot_widget2.plot(zs, u[k, :], pen=pg.mkPen(color=(0, 0, 0)))
        for j in range(4):
            plot_widget.plot(ts, u[:, int(j * I / 4) + 1], pen=pg.mkPen(color=(0, 0, 0)))
    if mult_graph_i and mult_graph_k:
        u = compute_implicit(I, k_thermal_cond, alpha, c, I - 1, K - 1, l, T, K, u_0, R, beta, p, a)
        ts = np.arange(0, T, T / (K + 1))
        zs = np.arange(0, l, l / (I + 1))
        plot_widget.plot(ts, u[:, I - 1], pen=pg.mkPen(color=(0, 0, 0)))
        plot_widget2.plot(zs, u[K - 1, :], pen=pg.mkPen(color=(0, 0, 0)))
        for j in range(4):
            plot_widget2.plot(zs, u[int(j * K / 4) + 1, :], pen=pg.mkPen(color=(0, 0, 0)))
        for j in range(4):
            plot_widget.plot(ts, u[:, int(j * I / 4) + 1], pen=pg.mkPen(color=(0, 0, 0)))


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
    mult_graph_k = False
    mult_graph_i = False
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
    else:
        mult_graph_i = True
    if form.lineEdit_12.text() != '':
        k = int(form.lineEdit_12.text())
    else:
        mult_graph_k = True
    if form.lineEdit_9.text() != '':
        I = int(form.lineEdit_9.text())
    if form.lineEdit_10.text() != '':
        K = int(form.lineEdit_10.text())
    if form.lineEdit_13.text() != '':
        p = float(form.lineEdit_13.text())
    if form.lineEdit_14.text() != '':
        a = float(form.lineEdit_14.text()) * R
    plot_widget = pg.plot(title="График для i")
    plot_widget.showGrid(x=True, y=True)
    plot_widget.setBackground('w')
    plot_widget2 = pg.plot(title="График для k")
    plot_widget2.setBackground('w')
    plot_widget2.showGrid(x=True, y=True)
    if not mult_graph_i and not mult_graph_k:
        u = compute_krank(I, k_thermal_cond, alpha, c, l, T, K, u_0, R, beta, p, a)
        zs = np.arange(0, l, l / (I + 1))
        ts = np.arange(0, T, T / (K + 1))
        plot_widget.plot(ts, u[:, i], pen=pg.mkPen(color=(0, 0, 0)))
        plot_widget2.plot(zs, u[k, :], pen=pg.mkPen(color=(0, 0, 0)))
    if not mult_graph_i and mult_graph_k:
        u = compute_krank(I, k_thermal_cond, alpha, c, i, T, K, u_0, R, beta, p, a)
        ts = np.arange(0, T, T / (K + 1))
        zs = np.arange(0, l, l / (I + 1))
        plot_widget.plot(ts, u[:, i], pen=pg.mkPen(color=(0, 0, 0)))
        plot_widget2.plot(zs, u[K-1, :], pen=pg.mkPen(color=(0, 0, 0)))
        for j in range(4):
            plot_widget2.plot(zs, u[int(j * K / 4) + 1, :], pen=pg.mkPen(color=(0, 0, 0)))
    if mult_graph_i and not mult_graph_k:
        u = compute_krank(I, k_thermal_cond, alpha, c, l, T, K, u_0, R, beta, p, a)
        ts = np.arange(0, T, T / (K + 1))
        zs = np.arange(0, l, l / (I + 1))
        plot_widget.plot(ts, u[:, I-1], pen=pg.mkPen(color=(0, 0, 0)))
        plot_widget2.plot(zs, u[k, :], pen=pg.mkPen(color=(0, 0, 0)))
        for j in range(4):
            plot_widget.plot(ts, u[:, int(j * I / 4) + 1], pen=pg.mkPen(color=(0, 0, 0)))
    if mult_graph_i and mult_graph_k:
        u = compute_krank(I, k_thermal_cond, alpha, c, l, T, K, u_0, R, beta, p, a)
        ts = np.arange(0, T, T / (K + 1))
        zs = np.arange(0, l, l / (I + 1))
        plot_widget.plot(ts,  u[:, I-1], pen=pg.mkPen(color=(0, 0, 0)))
        plot_widget2.plot(zs, u[K-1, :], pen=pg.mkPen(color=(0, 0, 0)))
        for j in range(4):
            plot_widget2.plot(zs, u[int(j * K / 4) + 1, :], pen=pg.mkPen(color=(0, 0, 0)))
        for j in range(4):
            plot_widget.plot(ts, u[:, int(j * I / 4) + 1], pen=pg.mkPen(color=(0, 0, 0)))


def test_cranck_button():
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
    if form.lineEdit_13.text() != '':
        p = float(form.lineEdit_13.text())
    if form.lineEdit_14.text() != '':
        a = float(form.lineEdit_14.text()) * R
    I = 5
    K = 10
    fig, axes = plt.subplots(nrows=1, ncols=2)
    for i in range(10):
        u = compute_krank(I, k_thermal_cond, alpha, c, l, T, K, u_0, R, beta, p, a)
        for_i_xs = np.linspace(0, T, K+1)
        for_i_us = u[:, I // 2]
        for_n_xs = np.linspace(0, l, I+1)
        for_n_us = u[K // 2, :]
        axes[0].plot(for_i_xs, for_i_us, label='I = ' + str(I) + ', K = ' + str(K))
        axes[1].plot(for_n_xs, for_n_us, label='I = ' + str(I) + ', K = ' + str(K))
        I *= 2
        K *= 2
    name_0 = "Распределение температуры в точке z={:4.2f} см".format(l/2)
    name_1 = "Распределение температуры в момент времени t={:4.2f} с".format(T / 2)
    axes[0].set(title=name_0, xlabel="Время, с", ylabel="Темпераура, град")
    axes[1].set(title=name_1, xlabel="Расстояние, см", ylabel="Темпераура, град")
    axes[0].grid()
    axes[1].grid()
    axes[0].legend()
    axes[1].legend()
    fig.show()


def shod_cranck_button():
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
    if form.lineEdit_13.text() != '':
        p = float(form.lineEdit_13.text())
    if form.lineEdit_14.text() != '':
        a = float(form.lineEdit_14.text()) * R
    I = 100
    K = 10
    fig, axes = plt.subplots(nrows=1, ncols=1)
    us = []
    for i in range(13):
        u = compute_krank(I, k_thermal_cond, alpha, c, l, T, K, u_0, R, beta, p, a)
        for_i_xs = np.linspace(0, l, I + 1)
        for_i_us = u[K//2, :]
        axes.plot(for_i_xs, for_i_us, label='K = ' + str(K))
        us.append(u[K // 2, I // 2])
        print("\nh_t = " + str(T/K))
        print("u_ht = " + str(us[-1]) + "\n")
        K *= 2
    name_0 = "Распределение температуры в точке z={:4.2f} см".format(l / 2)
    axes.set(title=name_0, xlabel="Время, с", ylabel="Темпераура, град")
    axes.grid()
    axes.legend()
    fig.show()





Form, Window = uic.loadUiType('interface.ui')
app = QApplication([])
window = Window()
form = Form()
form.setupUi(window)
form.pushButton.clicked.connect(button_clicked)
form.pushButton_2.clicked.connect(button_clicked_2)
form.pushButton_4.clicked.connect(test_cranck_button)
form.pushButton_6.clicked.connect(shod_cranck_button)
window.show()
app.exec_()