import numpy as np
import math
import matplotlib.pyplot as plt

def create_n(count_t, count_r, N):  # начальное распределение населенностей
    array = np.zeros((count_t, count_r, 5))
    for l in range(count_r):
        array[count_t - 1][l][0] = N
    return array


def count_alpha(alpha, R, R2, N, curr_m, count_r):
    for l in range(count_r):
        alpha[curr_m][l] = (n[curr_m][l][0] + R*n[curr_m][l][1] + R2 * n[curr_m][l][3])/N


def create_F(count_r, count_t, F):
    array = np.zeros((count_t, count_r))
    for m in range(count_t):
        if m >= 0: # Задаем длительность импульса
            array[m][0] = F
    return array


def F_matrix(F, R, R_t, T, T_2, T_st, T_ts):
    Matrix = np.matrix([[-F, 1, 0, T_ts, 0],
                          [F, -1 - F*R - T_st, T, 0, 0],
                          [0, F * R, -T, 0, 0],
                          [0, T_st, 0, -F*R_t - T_ts, T_2],
                          [0, 0, 0, F * R_t, -T_2]])
    return Matrix


def count_n(F_matrix, tau, curr_m, curr_l, n):
    M = np.eye(5) - F_matrix * tau
    M_inv = np.linalg.inv(M)
    n_old = []
    for k in range(5):
        n_old.append(n[curr_m][curr_l][k])
    n_new = np.zeros(5)
    np.dot(M_inv, n_old, n_new)
    for k in range(5):
        n[curr_m - 1][curr_l][k] = n_new[k]


def count_F(count_r, F, alpha, curr_m):
    for l in range(1, count_r):
        F[curr_m][l] = (1 - (alpha[curr_m][l - 1] + alpha[curr_m][l]) * h / 2) * F[curr_m][l - 1]


def plot_n(n, t, num):
    # Построим графики
    # n[m][0][0]
    n_t_0 = []
    n_t_1 = []
    n_t_2 = []
    n_t_3 = []
    n_t_4 = []
    Nsum = []
    for m in range(len(t)):
            n_t_0.append(n[count_t - 1 - m][num][0])
            n_t_1.append(n[count_t - 1 - m][num][1])
            n_t_2.append(n[count_t - 1 - m][num][2])
            n_t_3.append(n[count_t - 1 - m][num][3])
            n_t_4.append(n[count_t - 1 - m][num][4])
            Nsum.append(n_t_0[m] + n_t_1[m] + n_t_2[m])
    plt.plot(t, n_t_0, 'r', t, n_t_1, 'b', t, n_t_2, 'g', t, n_t_3, 'm',  t, n_t_4, 'c', t, Nsum, 'y')
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


def count_T_integral(T_integral, Iin, f0_start,f0_finish, count_r, count_t,R, R_t, T, T_2, T_st, T_ts, tau, n, coeff_F_into_I):
    jj = math.ceil(pow((f0_finish / f0_start), 1 / 3))
    for j in range(1, jj):
        f = f0_start * j ** 3  # f0 = F* k0/sigma0
        F = create_F(count_r, count_t,f)  # задаем начальные условия для потока падающего излучения -- сложность count_t
        alpha = np.zeros((count_t, count_r))
        for m in range(count_t - 1):
            for l in range(count_r):
                count_n(F_matrix(F[count_t - 1 - m][l],R, R_t, T, T_2, T_st, T_ts), tau, count_t - 1 - m, l, n)
            count_alpha(alpha, R, R_t, N, count_t - 2 - m, count_r)
            count_F(count_r, F, alpha, count_t - 2 - m)
        T_t = create_T_t(count_t, count_r, F, f)  # пропускание
        output = np.trapz(T_t, t)
        T_integral.append(output / (tpulse * k[0]))
        Iin.append(f * coeff_F_into_I)


def plot_T_integral(T_integral, Iin):
    plt.plot(Iin, T_integral, 'ro')
    plt.xlabel("Iin", fontsize=14, fontweight="bold")
    plt.ylabel("T_integral", fontsize=14, fontweight="bold")
    plt.xscale('log')
    #plt.ylim(0, 0.7)
    plt.show()
    print('Hello')
    print(T_integral)
    print(Iin)


def diff_time(F, n, sigma, t):
    F_1 = []
    F_0 = []
    F_2 = []
    my_t = []
    delta_arr = []
    # Cделаем счетчик
    for m in range(len(t)):
        #F_1.append(sigma[1] * n[count_t - 1 - m][int(count_r/2) - 1][1])
        #F_0.append(sigma[0] * n[count_t - 1 - m][int(count_r/2) - 1][0])
        n_1 = 0
        n_0 = 0
        n_2 = 0
        for l in range(len(r)):
            n_1 += n[count_t - 1 - m][l][1]
            n_0 += n[count_t - 1 - m][l][0]
            n_2 += n[count_t - 1 - m][l][3]
        F_1.append(sigma[1] * n_1)
        F_0.append(sigma[0] * n_0)
        F_2.append(sigma[2] * n_2)
        my_t.append(t[m])
        #delta_2 = F_1[m] + F_2[m] - F_0[m]
        delta_2 = F_1[m] - F_2[m]
        if m != 0:
            delta_arr.append(abs(delta_2))
    find_time = t[np.argmin(delta_arr[1:]) + 1]
    plt.plot(t, F_1, 'r', t, F_0, 'b')
    plt.xlabel("t", fontsize=14, fontweight="bold")
    plt.ylabel("Diff", fontsize=14, fontweight="bold")
    plt.show()
    return find_time


def integral_time(F, n, sigma, t):
    F_2 = []
    F_1 = []
    F_0 = []
    my_t = []
    delta_arr = []
    I_1 = []
    I_0 = []
    I_2 = []
    # Сделаем счетчик
    for m in range(len(t)):
        n_1 = 0
        n_0 = 0
        n_2 = 0
        for l in range(len(r)):
            n_1 += n[count_t - 1 - m][l][1]
            n_0 += n[count_t - 1 - m][l][0]
            n_2 += n[count_t - 1 - m][l][3]
        F_1.append(F[count_t - 1 - m][count_r - 1] * sigma[1] * n_1)
        F_0.append(F[count_t - 1 - m][count_r - 1] * sigma[0] * n_0)
        F_2.append(F[count_t - 1 - m][count_r - 1] * sigma[2] * n_2)
        my_t.append(t[m])
        I_1.append(np.trapz(F_1, my_t))
        I_0.append(np.trapz(F_0, my_t))
        I_2.append(np.trapz(F_2, my_t))
        delta_2 = I_1[m] - I_2[m]
        if m != 0:
            delta_arr.append(abs(delta_2))
    find_time = t[np.argmin(delta_arr[1:]) + 1]
    #plt.plot(t, I_1, 'r', t, I_0, 'b')
    #plt.xlabel("t", fontsize=14, fontweight="bold")
    #plt.ylabel("Integral", fontsize=14, fontweight="bold")
    #plt.show()
    return find_time


def part_two(f0_start, count_r, count_t, R, R_t, T, T_2, T_st, T_ts, tau):
    F = create_F(count_r, count_t, f0_start)  # задаем начальные условия для потока падающего излучения -- сложность count_t
    alpha = np.zeros((count_t, count_r))
    for m in range(count_t - 1):
        for l in range(count_r):
            count_n(F_matrix(F[count_t - 1 - m][l],R, R_t, T, T_2, T_st, T_ts), tau, count_t - 1 - m, l, n)
        count_alpha(alpha, R, R_t,  N, count_t - 2 - m, count_r)
        count_F(count_r, F, alpha, count_t - 2 - m)
    itime = integral_time(F, n, sigma, t)
    #itime = 0
    dtime = diff_time(F, n, sigma, t)
    array = []
    array.append(itime)
    array.append(dtime)
    return array


#Зададим параметры системы

sigma = np.array([2.4, 30, 48]) # sigma0 и sigmaS --- 2.4 и 30 --- сечения поглощения
k = np.array([0.01, 0.001, 0, 1, 1]) #t0 сидит здесь --- 0,53 * 330 millisec и 1 picseс -- ??? большая разница -- скорости переходов =1/t

degree = np.array([-22, 12, -6]) # first - degree of sigma , second - degree of k, third - degree of I

F01 = (k[0] / sigma[0]) * 10 ** int(degree[1] - degree[0]) # относительно этой величины измеряем поток
Fcr = (k[0] / sigma[1]) * 10 ** int(degree[1] - degree[0]) # характерная величина системы, отвечающая за переключение механизмов поглощения
#Найдем myF0
I_J_start = 0.01 #microJ
I_J_finish = 1000
h = 1.0545726 * 10**(-34)
lamb = 532 * 10**(-9) # длина волны падающего излучения
c = 3 * 10**8
omega = 2 * np.pi * c / lamb
w0 = 81 * 10**(-6) # параметр Гауссового пучка
tp = 8 * 10**(-9) # длительность импульса
coeff_I_into_F = (10**(-6))/(h * omega * w0**2 * tp)
coeff_F_into_I = F01/coeff_I_into_F
F_start = I_J_start * coeff_I_into_F # ??
F_finish = I_J_finish * coeff_I_into_F

f0_start = F_start / F01# интенсивность на входе в долях f0
f0_finish = F_finish / F01
fcr = Fcr / F01

count = 30 # число отрезков
param = 300
count_r = count
count_t = int(count * param) # число шагов ---N


tpulse = 8000
dt = tpulse/count_t
tau = dt*k[0]
rL = 0.01 # толщина образца
dz = rL/count_r
N = -np.log(0.7)/(sigma[0] * rL) # число электронов в начальный момент - задает плотность образца  -
h = dz * N * sigma[0] # шаг по координате
print('tau:', tau) # Шаг должен быть меньше самого маленького времени
print('1*k[0]:', 1*k[0])

R = sigma[1]/sigma[0]
R_t = sigma[2]/sigma[0]
T = k[3]/k[0]
T_2 = k[4]/k[0]
T_st = k[1]/k[0]
T_ts = k[2]/k[0]
r = np.linspace(0, rL*(N * sigma[0]), count_r)
t = np.linspace(0, tpulse*k[0], count_t)
n = create_n(count_t, count_r, N) #  задаем начальные условия для населенности -- сложность count_r

T_integral = []
Iin = []
count_T_integral(T_integral, Iin, f0_start, f0_finish,count_r,count_t,R, R_t, T, T_2, T_st, T_ts, tau,n, coeff_F_into_I)
plot_T_integral(T_integral, Iin)

jj = math.ceil(pow((f0_finish / f0_start), 1 / 3))
dtime_ = []
itime_ = []
f_ = []
dtime_a = [] # не забыть домножить
f_a = []
i_ = []
i_a = []
for j in range(1, jj):
    f = f0_start * j ** 3  # f0 = F* k0/sigma0
    dtime_.append(part_two(f, count_r, count_t, R, R_t, T, T_2, T_st, T_ts, tau)[1]/k[0])
    itime_.append(part_two(f, count_r, count_t, R, R_t, T, T_2, T_st, T_ts, tau)[0]/k[0])
    podln = (f*(1+R)/(fcr*R))/(f/fcr - 1)
    if podln > 0:
        coeff = 1 / (k[0] * (1 + f / (fcr * R))) * math.log(podln)
        dtime_a.append(coeff)
        f_a.append(f)
        i_a.append(f*coeff_F_into_I)
    f_.append(f)
    i_.append(f * coeff_F_into_I)
# построим гафик my_itime от f
#plt.plot(i_, dtime_, 'ro', i_a, dtime_a,'r')
plt.plot(i_, dtime_, 'ro')
plt.xlabel("I, microJ", fontsize=14, fontweight="bold")
plt.ylabel("D_time, picsec", fontsize=14, fontweight="bold")
#plt.legend(['chisl ','analitic '], loc=2)
#plt.ylim(0,0.1)
plt.show()
# График для интегрального времени
plt.plot(i_, itime_, 'ro')
plt.xlabel("I, microJ", fontsize=14, fontweight="bold")
plt.ylabel("I_time, picsec", fontsize=14, fontweight="bold")
#plt.legend(['chisl ','analitic '], loc=2)
#plt.ylim(0,0.1)
plt.show()