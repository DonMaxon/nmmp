import numpy as np


def thomas_method(a, b, c, f):
	"""
	На вход поступают массивы с элементами {a_i}, {b_i}, {c_i} и {f_i}
	{a_i} - массив длины I, содержащий элементы, находящиеся под главной диагональю
	{b_i} - массив длины I, содержащий элементы, находящиеся над главной диагональю
	{c_i} - массив длины I+1, содержащий элементы, находящиеся на главной диагонали
	{f_i} - массив длины I+1, содержащий элементы правой части системы
	На основе этих массивов находятся прогоночные коэффициенты
	С помощью прогоночных коэффициентов формируется массив, являющийся решением системы
	Функция возвращает массив-решение системы
	"""
	if len(c) != len(f) or len(a) != len(b) or len(c) - len(a) != 1:
		raise Exception(f"""Неверные размеры массивов входных данных!\n
						len(a) = {len(a)}\n
						len(b) = {len(b)}\n
						len(c) = {len(c)}\n
						len(f) = {len(f)}""")
	length = len(a)
	u = np.zeros(length+1)
	alphas = np.array([b[0]/c[0]])
	betas = np.array([f[0]/c[0]])
	for i in range(1, length):
		alphas = np.append(alphas, b[i]/(c[i]-a[i]*alphas[i-1]))
		betas = np.append(betas, (f[i]-a[i]*betas[i-1])/(c[i]-a[i]*alphas[i-1]))
	betas = np.append(betas, (f[-1]-a[-1]*betas[-1])/(c[-1]-a[-1]*alphas[-1]))
	u[-1] = betas[-1]
	for i in range(length-1, -1, -1):
		u[i] = betas[i]-alphas[i]*u[i+1]
	return u
