import numpy as np


def thomas_method(a, c, a1, c1, f):
	"""
	На вход поступают элементы a, c, a1, c1 матрицы системы и вектор правой части f
	На основе этих массивов находятся прогоночные коэффициенты
	С помощью прогоночных коэффициентов формируется массив, являющийся решением системы
	Функция возвращает массив-решение системы
	"""
	length = len(f)
	u = np.zeros(length)
	p = np.array([-a1/c1])
	q = np.array([f[0]/c1])
	print(p)
	print(q)
	for i in range(1, length):
		p = np.append(p, -a/(c+a*p[i-1]))
		q = np.append(q, (f[i]-a*q[i-1])/(c+a*p[i-1]))
	print(p)
	print(q)
	q = np.append(q, (f[length-1]-a*q[length-2])/(c+a*p[length-2]))
	u[length-1] = q[length-1]
	for i in range(2, length+1, 1):
		u[length-i] = q[length-i]+p[length-i]*u[length-i+1]
	return u
