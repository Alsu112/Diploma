import numpy as np
import math
import matplotlib.pyplot as plt
import sympy as sp
from sympy.abc import x
import cmath

def create_n(count_t, count_r, N):  # начальное распределение населенностей
    array = np.zeros((count_t, count_r, 3))
    for l in range(count_r):
        array[count_t - 1][l][0] = N
    return array


def count_alpha(alpha, R, N, curr_m, count_r):
    # alpha[curr_m, :] = (n[curr_m, :, 0] + R*n[curr_m, :, 1])/N
    for l in range(count_r):
        alpha[curr_m][l] = (n[curr_m][l][0] + R*n[curr_m][l][1])/N

def create_F(count_r, count_t, F):
    array = np.zeros((count_t, count_r))
    for m in range(count_t):
        if m >= 0: # Задаем длительность импульса
            array[m][0] = F
    return array


def F_matrix(F, R, T):
    Matrix = np.matrix([[-F, 0, 0],
                          [F, -F * R, T],
                          [0, F * R, -T]])
    return Matrix


def count_n(F_matrix, tau, curr_m, curr_l, n):
    M = np.eye(3) - F_matrix * tau
    M_inv = np.linalg.inv(M)
    n_old = []
    for k in range(3):
        n_old.append(n[curr_m][curr_l][k])
    # n_old = np.copy(n[curr_m,curr_l,:)
    n_new = np.zeros(3)
    np.dot(M_inv, n_old, n_new)
    for k in range(3):
        n[curr_m - 1][curr_l][k] = n_new[k]
    # n[curr_m - 1, curr_l, :] = n_new


def count_F(count_r, F, alpha, curr_m):
    for l in range(1, count_r):
        F[curr_m][l] = (1 - (alpha[curr_m][l - 1] + alpha[curr_m][l]) * h / 2) * F[curr_m][l - 1]


def plot_n(n, t, num):
    # Построим графики
    # n[m][0][0]
    n_t_0 = []
    n_t_1 = []
    n_t_2 = []
    Nsum = []
    for m in range(len(t)):
            n_t_0.append(n[count_t - 1 - m][num][0])
            n_t_1.append(n[count_t - 1 - m][num][1])
            n_t_2.append(n[count_t - 1 - m][num][2])
            Nsum.append(n_t_0[m] + n_t_1[m] + n_t_2[m])
    plt.plot(t, n_t_0, 'r', t, n_t_1, 'b', t, n_t_2, 'g', t , Nsum, 'y')
    plt.xlabel("Время", fontsize=14, fontweight="bold")
    plt.ylabel("Населенности", fontsize=14, fontweight="bold")
    plt.show()


def plot_F(t, F):
    F_t = []
    for m in range(len(t)):
        F_t.append(F[count_t - 1 - m][count_r - 1])
    plt.plot(t, F_t, 'ro') # график зависимости F(t) на выходе из образца
    #plt.ylim(0, 0.01)
    plt.xlabel("Время", fontsize=14, fontweight="bold")
    plt.ylabel("Fout", fontsize=14, fontweight="bold")
    plt.show()


def create_T_t(my_count_t, my_count_r, F, f0):
    array = []
    for m in range(my_count_t):
        array.append(F[my_count_t - 1 - m][my_count_r - 1] / f0)
    return array


def count_T_integral(T_integral, Iin, f0_start,f0_finish, count_r, count_t, R, T, tau, n):
    jj = math.ceil(pow((f0_finish / f0_start), 1 / 3))
    for j in range(1, jj):
        f = f0_start * j ** 3  # f0 = F* k0/sigma0
        F = create_F(count_r, count_t,f)  # задаем начальные условия для потока падающего излучения -- сложность count_t
        alpha = np.zeros((count_t, count_r))
        for m in range(count_t - 1):
            for l in range(count_r):
                count_n(F_matrix(F[count_t - 1 - m][l], R, T), tau, count_t - 1 - m, l, n)
            count_alpha(alpha, R, N, count_t - 2 - m, count_r)
            count_F(count_r, F, alpha, count_t - 2 - m)
        T_t = create_T_t(count_t, count_r, F, f)  # пропускание
        output = np.trapz(T_t, t)
        #print("Output : ", output / tpulse)
        #print(output / tpulse)
        print(f * coeff_F_into_I)
        T_integral.append(output / tpulse)
        Iin.append(f * coeff_F_into_I)


def plot_T_integral(T_integral, Iin):
    plt.plot(Iin, T_integral, 'ro')
    plt.xlabel("Iin", fontsize=14, fontweight="bold")
    plt.ylabel("T_integral", fontsize=14, fontweight="bold")
    plt.xscale('log')
    plt.ylim(0, 0.7)
    plt.show()
    print('Hello')
    print(T_integral)
    print(Iin)


def diff2_time(F, n, sigma):
    F_1 = []
    F_0 = []
    delta_arr = []
    for m in range(len(t)):
        F_1.append(F[count_t - 1 - m][count_r - 1] * sigma[1] * n[count_t - 1 - m][count_r - 1][1])
        F_0.append(F[count_t - 1 - m][count_r - 1] * sigma[0] * n[count_t - 1 - m][count_r - 1][0])
        delta = F_1[m] - F_0[m]
        delta_arr.append(abs(delta))
    #print('min: ', np.argmin(delta_arr))
    #print('size(delta)', len(delta_arr))
    plt.plot(t, F_1, 'r', t, F_0, 'b')
    plt.xlabel("t", fontsize=14, fontweight="bold")
    plt.ylabel("Diff", fontsize=14, fontweight="bold")
    plt.show()
    find_time = t[np.argmin(delta_arr[1:]) + 1]
    return find_time


def diff_time(F, n, sigma, t):
    F_1 = []
    F_0 = []
    my_t = []
    delta_arr = []
    for m in range(len(t)):
        F_1.append(F[count_t - 1 - m][count_r - 1] * sigma[1] * n[count_t - 1 - m][count_r - 1][1])
        F_0.append(F[count_t - 1 - m][count_r - 1] * sigma[0] * n[count_t - 1 - m][count_r - 1][0])
        my_t.append(t[m])
        delta_2 = F_1[m] - F_0[m]
        if m != 0:
            delta_arr.append(abs(delta_2))
    #print('min: ', np.argmin(delta_arr))
    #print('size(delta)', len(delta_arr))
    #plt.plot(t, I_1, 'r', t, I_0, 'b')
    #plt.xlabel("t", fontsize=14, fontweight="bold")
    #plt.ylabel("Integral", fontsize=14, fontweight="bold")
    #plt.show()
    find_time = t[np.argmin(delta_arr[1:]) + 1]
    return find_time


def integral_time(F, n, sigma, t):
    F_1 = []
    F_0 = []
    my_t = []
    delta_arr = []
    I_1 = []
    I_0 = []
    for m in range(len(t)):
        F_1.append(F[count_t - 1 - m][count_r - 1] * sigma[1] * n[count_t - 1 - m][count_r - 1][1])
        F_0.append(F[count_t - 1 - m][count_r - 1] * sigma[0] * n[count_t - 1 - m][count_r - 1][0])
        my_t.append(t[m])
        I_1.append(np.trapz(F_1, my_t))
        I_0.append(np.trapz(F_0, my_t))
        delta_2 = I_1[m] - I_0[m]
        if m != 0:
            delta_arr.append(abs(delta_2))
    #print('min: ', np.argmin(delta_arr))
    #print('size(delta)', len(delta_arr))
    #plt.plot(t, I_1, 'r', t, I_0, 'b')
    #plt.xlabel("t", fontsize=14, fontweight="bold")
    #plt.ylabel("Integral", fontsize=14, fontweight="bold")
    #plt.show()
    find_time = t[np.argmin(delta_arr[1:]) + 1]
    return find_time


def part_two(f0_start, count_r, count_t, R, T, tau):
    F = create_F(count_r, count_t, f0_start)  # задаем начальные условия для потока падающего излучения -- сложность count_t
    alpha = np.zeros((count_t, count_r))
    for m in range(count_t - 1):
        for l in range(count_r):
            count_n(F_matrix(F[count_t - 1 - m][l], R, T), tau, count_t - 1 - m, l, n)
        count_alpha(alpha, R, N, count_t - 2 - m, count_r)
        count_F(count_r, F, alpha, count_t - 2 - m)
    #dtime = diff_time(F, n, sigma)
    #print('dtime: ', dtime)

    #itime = integral_time(F, n, sigma, t)
    dtime = diff_time(F, n, sigma, t)
    return dtime
    #print(itime)
    #print('itime: ', itime)


def open_file(alpha, betta):
    f = open('points.dat', 'r')
    lines = f.readlines()
    for x in lines:
        alpha.append(float(x.split("; ")[0]))
        betta.append(float(x.split()[1]))
    f.close()


#Зададим параметры системы

sigma = np.array([2.4, 30]) # sigma0 и sigmaS --- 2.4 и 30
k = np.array([0, 1]) #t0 сидит здесь --- 0,53 * 330 millisec и 1 picseс -- ??? большая разница

degree = np.array([-22, 12, -6]) # first - degree of sigma , second - degree of k, third - degree of I

F01 = (k[1] / sigma[0]) * 10 ** (34) # относительно этой величины измеряем поток
Fcr = (k[1] / sigma[1]) * 10 ** (34)
#Найдем myF0
I_J_start = 0.01 #microJ
I_J_finish = 105
h = 1.0545726 * 10**(-34)
lamb = 532 * 10**(-9)
c = 3 * 10**8
omega = 2 * np.pi * c / lamb
w0 = 81 * 10**(-6)
tp = 70 * 10**(-12)
coeff_I_into_F = (10**(-6))/(h * omega * w0**2 * tp)
#coeff_F_into_I = 400
coeff_F_into_I = F01/coeff_I_into_F
F_start = I_J_start * coeff_I_into_F # ??
F_finish = I_J_finish * coeff_I_into_F

f0_start = F_start / F01# интенсивность на входе в долях f0
f0_finish = F_finish / F01


count = 140 # число отрезков
param = 2
count_r = count
count_t = count * param # число шагов ---N


tpulse = 70
dt = tpulse/count_t
tau = dt*k[1]
print('tau:' ,tau)
print('tpulse/k[1]:', tpulse/k[1])
rL = 10 * 0.01 # толщина образца
dz = rL/count_r
N = -np.log(0.7)/sigma[0] / rL # число электронов в начальный момент - задает плотность образца  -
h = dz * N * sigma[0] # шаг по координате

R = sigma[1]/sigma[0]
T = 1

r = np.linspace(0, rL/(N * sigma[0]), count_r)
t = np.linspace(0, tpulse/k[1], count_t)
n = create_n(count_t, count_r, N) #  задаем начальные условия для населенности -- сложность count_r

T_integral = []
Iin = []
count_T_integral(T_integral, Iin, f0_start, f0_finish,count_r,count_t,R,T,tau,n)
plot_T_integral(T_integral, Iin)
#part_two(f0_finish/100,count_r,count_t,R,T,tau)
#print('f0_start * 1000', f0_start * 1000)
#print('f0_start', f0_start)
#print('f0_finish', f0_finish)
jj = math.ceil(pow((f0_finish / f0_start), 1 / 3))
itime_ = []
f_ = []
for j in range(1, jj):
    f = f0_start * j ** 3  # f0 = F* k0/sigma0
    itime_.append(part_two(f, count_r, count_t, R, T, tau))
    f_.append(f)
print(f_)
print(itime_)
plt.plot(f_, itime_, 'ro')
plt.xlabel("f", fontsize=14, fontweight="bold")
plt.ylabel("I_time", fontsize=14, fontweight="bold")
plt.show()
# Cравним со статьей
# Data for plotting
alpha = []
betta = []
open_file(alpha, betta)
# Data for plotting

fig, ax = plt.subplots()
ax.plot(Iin, T_integral,'o-', color='blue')
ax.plot(alpha, betta,'o', color='red')
ax.set_xscale('log')

ax.set(xlabel='Incident Energy (muJ)', ylabel='Transmission',
       title='Comparison with the article for picsec')
#ax.grid()
ax.grid(True, color = "grey", linewidth = "1", linestyle = "-.")
ax.legend(['Simulation', 'Article'])
fig.savefig("test.png")
plt.show()