import numpy as np
import matplotlib.pyplot as plt

N=120
L_uranium = 100
L_water = 20  #
S_rem_fast = np.zeros(N)
S_rem_fast[0:L_uranium] = 1.79E-02  # Активная зона (диоксид урана)
S_rem_fast[L_uranium:N] = 5.72E-02  # Отражатель (вода)

S_gen_fast = np.zeros(N)
S_gen_fast[0:L_uranium] = 1.54E-02  # Активная зона
S_gen_fast[L_uranium:N] = 0  # Отражатель

S_s = np.zeros(N)
S_s[0:L_uranium] = 1.43E-03  # Активная зона
S_s[L_uranium:N] = 5.67E-02  # Отражатель

D_fast = np.zeros(N)
D_fast[0:L_uranium] = 1.07E+00  # Активная зона
D_fast[L_uranium:N] = 1.33E+00  # Отражатель


S_rem_th = np.zeros(N)
S_rem_th[0:L_uranium] = 2.49E-01  # Активная зона
S_rem_th[L_uranium:N] = 1.90E-02  # Отражатель

S_gen_th = np.zeros(N)
S_gen_th[0:L_uranium] = 4.68E-01  # Активная зона
S_gen_th[L_uranium:N] = 0  # Отражатель

D_th = np.zeros(N)
D_th[0:L_uranium] = 5.22E-01  # Активная зона
D_th[L_uranium:N] = 2.85E-01  # Отражатель

# Потоки нейтронов
Flow_fast = np.ones(N)  # Поток быстрой группы
Flow_th = np.ones(N)  # Поток тепловой группы


k_eff = 1.0
eps = 1e-5


def solve_diffusion(flow, source, diff_len, S_rem, N):
    a = np.zeros(N)
    for i in range(1, N):
        a[i] = -(diff_len[i - 1] + diff_len[i]) / 2

    c = np.zeros(N)
    for i in range(0, N - 1):
        c[i] = -(diff_len[i + 1] + diff_len[i]) / 2

    q = source.copy()

    b = np.zeros(N)
    b[0] = S_rem[0] - c[0]
    b[N - 1] = S_rem[N - 1] - a[N - 1] + 0.5
    b[1:N - 1] = S_rem[1:N - 1] - c[1:N - 1] - a[1:N - 1]

    A = np.zeros(N)
    A[0] = q[0] / b[0]

    B = np.zeros(N)
    B[0] = c[0] / b[0]

    for i in range(N - 1):
        A[i + 1] = (q[i + 1] - a[i + 1] * A[i]) / (b[i + 1] - a[i + 1] * B[i])
        B[i + 1] = c[i + 1] / (b[i + 1] - a[i + 1] * B[i])

    flow[N - 1] = A[N - 1]
    for i in range(N - 2, -1, -1):
        flow[i] = A[i] - B[i] * flow[i + 1]


# Основной цикл для нахождения k_eff
def find_k_eff(max_iter=1000):
    global k_eff, Flow_fast, Flow_th

    for iteration in range(max_iter):
        last_norm = np.linalg.norm(Flow_th * S_gen_th + Flow_fast * S_gen_fast)
        k = k_eff

        # Источник нейтронов
        Q = (S_gen_fast * Flow_fast + S_gen_th * Flow_th) / k_eff

        # Решение уравнений диффузии для быстрой и тепловой групп
        solve_diffusion(Flow_fast, Q, D_fast, S_rem_fast, N)
        solve_diffusion(Flow_th, Flow_fast * S_s, D_th, S_rem_th, N)

        # Обновление k_eff
        k_eff = k_eff * np.linalg.norm(Flow_th * S_gen_th + Flow_fast * S_gen_fast) / last_norm


    print(f"Эффективный коэффициент размножения k_eff = {k_eff}")


    x = np.arange(N)
    plt.figure(figsize=(10, 6))
    plt.plot(x, Flow_th, color='green', linewidth=2, label="Тепловая группа")
    plt.plot(x, Flow_fast, color='red', linewidth=2, label="Быстрая группа")
    plt.xlabel('Расстояние (см)', fontsize=12)
    plt.ylabel('Поток нейтронов', fontsize=12)
    plt.title('Распределение потока нейтронов', fontsize=14)
    plt.legend(fontsize=12)
    plt.grid(True)
    plt.show()


# Основной блок
find_k_eff()

#Доп.задание:
'''
N_range = range(101, 151)  # Диапазон N (от 101 до 150)
k_eff_values = []

# Итерации по N
for N in N_range:
    k_eff = calculate_k_eff(N)
    k_eff_values.append(k_eff)

# Визуализация зависимости k_eff от N
plt.figure(figsize=(10, 6))
plt.grid()
plt.plot(N_range, k_eff_values, marker='o', linestyle='-', color='b')
plt.xlabel('N (толщина отражателя)', fontsize=12)
plt.show() '''


