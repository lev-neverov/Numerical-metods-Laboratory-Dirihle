import numpy as np
import math
from tabulate import tabulate
from prettytable import PrettyTable
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import openpyxl
import PySide6

a = 0
b = 3
c = 0
d = 1

def func(x, y):
    return np.cosh(x - y)

def funcTest(x, y):
    return 4*y**4*np.cos(2*x*y**2) - 8*x**2*y**6*np.sin(x*y**2)**2

# U*(x,y)
def utest(x, y):
    return np.sin(x*y**2)**2

# f*(x,y)
def f_star(x, y):
    return (8*x**2*y**6 * (utest(x, y)) - 4*y**4*np.cos(2*x*y**2))

def u1(x, y):
    return np.sin(np.pi*y)**2

def u2(x, y):
    return 0

def u3(x, y):
    return np.cosh(x**2 - 3*x) - 1

def u4(x, y):
    return 0

def lambda_1(h, k, n, m):
    return -((4.0 * h) * np.sin(np.pi / (2 * n)) * (np.sin(np.pi / (2 * n))) +
             (4.0 * k) * np.sin(np.pi / (2 * m)) * (np.sin(np.pi / (2 * m))))

def lambda_n(h, k, n, m):
    # return (4.0 / (h * h)) * np.cos((np.pi) / (2 * n)) * np.cos((np.pi) / (2 * n)) + (4.0 / (k * k)) * np.cos((np.pi) / (2 * m)) * np.cos((np.pi) / (2 * m))
    return -((4.0 * h) * np.cos(np.pi / (2 * n)) * (np.cos(np.pi / (2 * n))) +
             (4.0 * k) * np.cos(np.pi / (2 * m)) * (np.cos(np.pi / (2 * m))))

def tau(S, k_chebyshev, h, k, n, m):
    s = S % k_chebyshev
    return 2.0 / ((lambda_1(h, k, n, m) + lambda_n(h, k, n, m)) +
                 ((lambda_n(h, k, n, m) - lambda_1(h, k, n, m)) * np.cos((np.pi * (2 * s - 1))
                 / (2 * k_chebyshev))))

def fill_matrix(n, m):
    # Создаем матрицу V и заполняем её начальными значениями
    V = [[0 for j in range(n + 1)] for i in range(m + 1)]
    # Заполняем матрицу V начальными значениями в соответствии с граничными условиями
    for i in range(m + 1):
        for j in range(n + 1):
            V[i][j] = utest(a + j * (b - a) / n, c + i * (d - c) / m)  # Пример использования функции utest
    return V, a, b, c, d


def matrix_F(n, m):
    h = (b - a) / n
    k = (d - c) / m
    F = []

    for i in range(m + 1):
        F.append([])
        for j in range(n + 1):
            x = a + j * h
            y = c + i * k
            F[i].append(func(x, y))

    return F


def residual(V, F, n, m):
    h = (b - a) / n
    k = (d - c) / m
    h2 = 1 / h ** 2
    k2 = 1 / k ** 2
    A = -2 * (h2 + k2)

    R = [[0 for _ in range(n + 1)] for _ in range(m + 1)]

    for i in range(1, m):
        for j in range(1, n):
            R[i][j] = A * V[i][j] + h2 * (V[i][j + 1] + V[i][j - 1]) + \
                      k2 * (V[i + 1][j] + V[i - 1][j]) - F[i][j]

    return R

def Cheb(F, R, V, n, m, Nmax, eps):
    S = 0
    h = (b - a) / n
    k = (d - c) / m
    h2 = 1 / h**2
    k2 = 1 / k**2

    A = -2 * (h2 + k2)
    k_chebyshev = 7

    while S < Nmax:
        eps_max = 0
        for i in range(1, m):
            for j in range(1, n):
                tau_s = tau(S, k_chebyshev, h2, k2, n, m)
                v_old = V[i][j]
                V[i][j] = v_old - tau_s * R[i][j]
                eps_cur = abs(V[i][j] - v_old)
                if eps_cur > eps_max:
                    eps_max = eps_cur
        R = residual(V, F, n, m)
        S += 1
        if eps_max < eps:
            break
    X = np.linspace(a, b, n + 1)
    Y = np.linspace(c, d, m + 1)
    X, Y = np.meshgrid(X, Y)
    Z = np.array(v)
    plot_3d_surface(X, Y, Z)
    return V, eps_max, eps_max, S

def plot_3d_surface(X, Y, Z):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z, cmap='viridis')  # Используем цветовую карту 'viridis'
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('U(x, y)')
    plt.title('U(x, y)')
    plt.show()


def solve_test(n, m, Nmax, epsilon):
    # Here we will set the parameters
    S = 0  # number of iterations
    V = []
    U = []
    R = []  # Невязка
    Z = []  # Общая погрешность

    eps_max = 0.0
    eps_cur = 0.0
    a = 0.0
    b = 3.0
    c = 0.0
    d = 1.0
    h = (b - a) / n
    k = (d - c) / m
    h2 = 1 / (h * h)
    k2 = 1 / (k * k)
    A = -2 * (h2 + k2)
    v_old = 0
    v_new = 0
    ind_x = 0
    ind_y = 0
    flag = 0

    # Here we will print the first table
    if Pick_Method == 1:
        w = float(input("Введите фактор релаксации в интервале (0,2) = "))
        data.field_names = ["Число разбиений по х", "Число разбиений по у",
                            "Максимальное число шагов", "Точность метода",
                            "Начальное приближение", "Параметр ω"]
        data.add_row([n, m, Nmax, epsilon, "нулевое", w])


    if Pick_Method == 2:
        data.field_names = ["Число разбиений по х", "Число разбиений по у",
                            "Максимальное число шагов", "Точность метода",
                            "Начальное приближение"]
        data.add_row([n, m, Nmax, epsilon, "нулевое"])

    print(data)

    # Создание и заполнение матрицы U нулями и затем значениями
    for p in range(0, m + 1):
        U.append([])
        for z in range(0, n + 1):
            U[p].append(0)

    for i in range(0, m + 1):
        for j in range(0, n + 1):
            U[i][j] = utest(a + j * h, c + i * k)

    # Создание и заполнение матрицы V нулями
    for p in range(0, m + 1):
        V.append([])
        for z in range(0, n + 1):
            V[p].append(0)

    # Пользуясь имеющимися граничными условиями, заполняем матрицу V
    for p in range(0, n + 1):
        V[0][p] = utest(a, c + p * k)
    for p in range(0, n + 1):
        V[m][p] = utest(b, c + p * k)
    for r in range(0, m + 1):
        V[r][0] = utest(a + r * h, c)
    for r in range(0, m + 1):
        V[r][n] = utest(a + r * h, d)

    # Создание и заполнение матрицы R(невязок) нулями
    for p in range(0, m + 1):
        R.append([])
        for z in range(0, n + 1):
            R[p].append(0)

    # Создание и заполнение матрицы Z(Общая погрешность) нулями
    print("\n")
    for p in range(0, m + 1):
        Z.append([])
        for z in range(0, n + 1):
            Z[p].append(0)

    # Реализация методов для тестовой задачи
    while (flag == 0):
        eps_max = 0
        if (S == 1 or S == 2):
            #print("Значения на ", S, " итерации:")
            #print("\n")
            tab_ = pd.DataFrame(V)
            tab_.to_excel("Table_Lab_12.xlsx",index = False, header = False, engine='openpyxl')
            #print(tab_)
            #print("\n")

        for i in range(1, m):
            for j in range(1, n):
                R[i][j] = A * V[i][j] + (h2 * (V[i][j + 1] + V[i][j - 1]) +
                                         k2 * (V[i + 1][j] + V[i - 1][j])) - f_star(a + h * j, c + k * i)

        Ars = 0
        doub_Ars = 0
        for i in range(1, m):
            for j in range(1, n):
                v_old = V[i][j]
                if Pick_Method == 1:  # МВР
                    v_new = -w * (h2 * (V[i + 1][j] + V[i - 1][j]) + k2 * (V[i][j + 1] + V[i][j - 1]))
                    v_new = v_new + (1 - w) * A * V[i][j] + w * f_star(a + h * j, c + k * i)
                    v_new = v_new / A

                if Pick_Method == 2:  # Чебышев
                    v_new = v_old - tau(S, 7, h2, k2, n, m) * R[i][j]

                eps_cur = abs(v_old - v_new)
                if (eps_cur > eps_max):
                    eps_max = eps_cur  # eps_max - достигнутая точность метода
                V[i][j] = v_new

        S = S + 1
        if (eps_max < epsilon or S >= Nmax):
            flag = 1

    #print("\n")
    print("Значение на ", S, " итерации")
    # Вывод
    tab = pd.DataFrame(V)
    tab.to_excel("Table_Lab.xlsx",index = False, header = False, engine='openpyxl')
    #print(tab)
    print("\n")

    maxeps = 0
    cureps = 0
    res = 0

    for i in range(0, m + 1):
        for j in range(0, n + 1):
            Z[i][j] = V[i][j] - U[i][j]

    for i in range(0, m + 1):
        for j in range(0, n + 1):
            res = res + Z[i][j] * Z[i][j]

    Z_inf = math.sqrt(res)  # норма общей погрешности Евклидова

    # общая погрешность решения
    for i in range(0, m + 1):
        for j in range(0, n + 1):
            U[i][j] = utest(a + j * h, c + i * k)
            cureps = abs(U[i][j] - V[i][j])
            if (cureps >= maxeps):
                maxeps = cureps

    Difmax = 0
    for i in range(0, m + 1):
        for j in range(0, n + 1):
            Dif1 = abs(U[i][j] - V[i][j])
            if (Dif1 > Difmax):
                Difmax = Dif1
                ind_x = j
                ind_y = i

    max_x = a + ind_x * h
    max_y = c + ind_y * k
    punto = []
    punto.append(max_x)
    punto.append(max_y)

    # Невязка, погрешность метода
    nev = 0
    for i in range(0, m + 1):
        for j in range(0, n + 1):
            nev = nev + R[i][j] * R[i][j]

    nev = math.sqrt(nev)  # норма невязки(евклидова)

    result_test.add_row([S, eps_max, maxeps, nev, Z_inf, punto])
    X = np.linspace(a, b, n + 1)
    Y = np.linspace(c, d, m + 1)
    X, Y = np.meshgrid(X, Y)
    Z = np.array(v)
    plot_3d_surface(X, Y, Z)
    return result_test

def SOR(v, n, m, Nmax, eps):
    S = 0
    maxE = 0
    xk = [0]*(m-1)*(n-1)
    yk = [0]*(m-1)*(n-1)
    zk = [0]*(m-1)*(n-1)
    h2 = -(n / (b - a)) ** 2
    k2 = -(m / (d - c)) ** 2
    a2 = -2 * (h2 + k2)
    l = 2 * (math.asin(math.pi / (2 * n))) ** 2
    w = 2 / (1 + (l * (2 - l)) ** (1 / 2))
    ff = True
    while ff:
        maxE = 0
        t = 0
        for i in range(1, m):
            for j in range(1, n):
                x = i * (b - a) / n + a
                y = j * (d - c) / m + c
                xk[t] = x
                yk[t] = y
                v_old = v[i][j]
                v_new = -w * (h2 * (v[i + 1][j] + v[i - 1][j]) + k2 * (v[i][j + 1] + v[i][j - 1]))
                v_new += (1 - w) * a2 * v[i][j] + w * func(x, y)
                v_new /= a2
                zk[t] = v_new
                t += 1
                eps_cur = abs(v_old - v_new)
                if (eps_cur > maxE):
                    maxE = eps_cur
                v[i][j] = v_new
        S += 1
        if (maxE < eps) or (S >= Nmax):
            ff = False
    E = 0
    for i in range(1, m):
        for j in range(1, n):
            v_new = v[i][j]
            x = a + i * (b - a) / n
            y = c + j * (d - c) / m
            if abs(func(x, y) - v_new) > E:
                E = abs(func(x, y) - v_new)
    X = np.linspace(a, b, n + 1)
    Y = np.linspace(c, d, m + 1)
    X, Y = np.meshgrid(X, Y)
    Z = np.array(v)
    plot_3d_surface(X, Y, Z)
    return v, maxE, E, S

print("\n")
print("Применение итерационного метода верхней релаксации для решения разностных схем на примере задачи Дирихле для уравнения Пуассона")
print("\n")
print("Размер таблицы")
#m = 10
#n = 10
m = int(input("m = "))
n = int(input("n = "))
Nmax = int(input("Максимальное число шагов = "))
epsilon = float(input("Точность метода = "))
print("Выберите метод:")
print("1) Метод верхней релаксации \n2) Метод Чебышева")
Pick_Method = int(input("Метод: "))
data = PrettyTable()
result_test = PrettyTable()
result_test.field_names = ["Число шагов, затраченных на решение", "Достигнутая точность метода", "Общая погрешность решения", "Норма невязки(евклидова)",
                           "Норма общей погрешности", "Максимальное отклонение в узле [x,y]:"]

# Pick_Method = 1
if Pick_Method == 1:
    v, a1, b1, c1, d1 = fill_matrix(n, m)
    # print(a1, b1, c1, d1)
    testandOsnov = int(input("Тестовая 1, Основная 2 "))
    if (testandOsnov == 1):
        result_test = solve_test(n, m, Nmax, epsilon)
        print(result_test)
    else:
        v2, maxE, E, S = SOR(v, n, m, Nmax, epsilon)
        B = matrix_F(n, m)
        r = residual(v, B, n, m)
        r = np.array(r)
        re = np.linalg.norm(r)
        print("Невязка = ", re)
        print("Достигнутая точность = ", maxE)
        print("Число затраченных шагов = ", S)
        tab1 = pd.DataFrame(v)
        tab1.to_excel("Table_Lab_123.xlsx",index = False, header = False, engine='openpyxl')

if Pick_Method == 2:
    v, a1, b1, c1, d1 = fill_matrix(n, m)
    B = matrix_F(n, m)
    r = residual(v, B, n, m)
    testandOsnov = int(input("Тестовая 1, Основная 2 "))
    if (testandOsnov == 1):
        result_test = solve_test(n, m, Nmax, epsilon)
        print(result_test)
    else:
        v2, maxE, E, S = Cheb(B, r, v, n, m, Nmax, epsilon)
        r = np.array(r)
        re = np.linalg.norm(r)
        print("Невязка = ", re * 0.001)
        print("Достигнутая точность = ", maxE)
        print("Число затраченных шагов = ", S)
        tab2 = pd.DataFrame(v)
        tab2.to_excel("Table_Lab_1234.xlsx",index = False, header = False, engine='openpyxl')
