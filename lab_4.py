import numpy as np
import matplotlib.pyplot as plt
import random


Sig_s = 0.3
Sig_a = 0.03
Nu_f = 2.4
Nu_f_Sig_f = 0.03
Sig_f = Nu_f_Sig_f / Nu_f
Sig_tot = Sig_s + Sig_a
Sig_cap = Sig_a - Sig_f

H = 100
h = 5
N_disc = int(H / h)

num_batches = int(input('Write a number of batches - queues: '))
neutrons = int(input('Write a number of neutrons for each batch: '))

No_active = int(0.5 * num_batches)
flux_arr = np.zeros((num_batches, N_disc))
reaction_count = np.zeros((2, N_disc))
sources = [np.random.uniform(0, N_disc, neutrons)]
keff = 1.0

def neutron_tracking(Y, keff):
    direction = 2 * random.random() - 1
    life_span = -np.log(random.random()) / Sig_tot
    Y_append = Y + direction * life_span

    num_fiss = 0
    if Y_append >= H:
        return num_fiss, None
    if Y_append <= 0:
        Y_append = -Y_append  # Отражение

    rand_inter = random.random()
    if rand_inter < Sig_s / Sig_tot:  # Рассеяние
        cell_id = int(Y_append / h)
        reaction_count[0, cell_id] += 1
        return neutron_tracking(Y_append, keff)

    if rand_inter < (Sig_s + Sig_cap) / Sig_tot:  # Поглощение
        cell_id = int(Y_append / h)
        reaction_count[1, cell_id] += 1
        return num_fiss, None

    # ДЕЛЕНИЕ
    fis_neutron = int(np.floor(Nu_f / keff))
    if random.random() < (Nu_f / keff) - fis_neutron:
        fis_neutron += 1
    num_fiss += fis_neutron
    return num_fiss, Y_append

for i in range(num_batches):
    generated = []
    for j in sources[-1]:
        created, position = neutron_tracking(j, keff)
        if position is not None:
            generated.extend([position] * created)

    for k in range(N_disc):
        flux_arr[i, k] = reaction_count[0, k] / (Sig_tot * h * neutrons)

    keff = keff * len(generated) / neutrons
    neutrons = len(generated)
    sources.append(generated)
    reaction_count[:,:]=0.

    print(f"Batch {i+1}: keff = {keff:}, neutrons: {neutrons}")


flux_avg = np.mean(flux_arr[No_active:], axis=0)

y = np.arange(N_disc) * h

plt.figure(figsize=(8, 5))
plt.plot(y, flux_avg, color='green', label='Метод Монте-Карло ', linewidth=2)
plt.xlabel('Координаты')
plt.ylabel(' поток')
plt.title('Gjnjrb')
plt.legend()
plt.grid()
plt.show()


#################################################################################
'''
""" Constants"""
L = 100
D = 1.0
Sigma_a = 0.03
nuSigma_f = 0.03

N = 100  # Number_of_lattice_points
dx = L / (N - 1)  # Step
x = np.linspace(0, L, N)  # Coordinates
A = np.zeros((N, N))
F = np.zeros((N, N))

for i in range(N):
    F[i, i] = nuSigma_f

for i in range(1, N - 1):
    A[i, i - 1] = -D / dx ** 2
    A[i, i] = 2 * D / dx ** 2 + Sigma_a
    A[i, i + 1] = -D / dx ** 2

# Left Boundary
A[0, 0] = D / dx ** 2 + Sigma_a
A[0, 1] = -D / dx ** 2

# Right boundary (zero flux condition)
A[-1, -1] = D / dx ** 2 + Sigma_a + 1/2
A[-1, -2] = -D / dx
def accelerated_source_iteration(A, nuSigma_f, delta_lambda=0.1, eps=1e-6, max_iter=1000):
    N = A.shape[0]
    phi = np.random.rand(N)
    keff = 1.0
    Lambda = 1

    for iteration in range(max_iter):
        lambda_ = keff + delta_lambda

        Q_f = nuSigma_f * phi / Lambda

        phi_0 = np.linalg.solve(A - 1 / lambda_ * F, Q_f)  # solve

        Lambda *= np.sum(nuSigma_f * phi_0) / np.sum(nuSigma_f * phi)
        keff_0 = (Lambda * lambda_) / (Lambda + lambda_)

        if np.abs(keff_0 - keff) < eps:  # check
            return keff_0, phi_0, iteration + 1

        keff = keff_0
        phi = phi_0


    return keff, phi, max_iter

keff, phi, iterations = accelerated_source_iteration(A, nuSigma_f)


print(f"Эффективный коэффициент размножения k_eff: {keff}")


"""delta_lambdas = np.linspace(0, 0.01, 10)
for delta_lambda in delta_lambdas:
    keff, phi, iterations = accelerated_source_iteration(A, nuSigma_f, delta_lambda)
    print(f"delta_lambda = {delta_lambda:.2f}, k_eff = {keff:.6f}, сходимость на {iterations}-й итерации")"""



flux_avg_norm = flux_avg / flux_avg[0]
phi_norm = phi / phi[0]


plt.figure(figsize=(8, 5))
plt.plot(y, flux_avg_norm, color='green', label='Метод Монте-Карло ', linewidth=2)
plt.plot(x, phi_norm, label='Асимптотический поток', color='red', linewidth=2)
plt.xlabel('Координаты')
plt.ylabel(' поток')
plt.title('Gjnjrb')
plt.legend()
plt.grid()
plt.show()  '''
"""import numpy as np
import matplotlib.pyplot as plt
import random

Sig_s = 0.3
Sig_a = 0.03
Nu_f = 2.4
Nu_f_Sig_f = 0.03
Sig_f = Nu_f_Sig_f / Nu_f
Sig_tot = Sig_s + Sig_a
Sig_cap = Sig_a - Sig_f

# Для дельта-трекинга используем максимальное сечение
Sig_max = 2 * Sig_tot  # Берем с запасом, как рекомендуется

H = 100
h = 5
N_disc = int(H / h)

num_batches = int(input('Write a number of batches - queues: '))
neutrons = int(input('Write a number of neutrons for each batch: '))

No_active = int(0.5 * num_batches)
flux_arr = np.zeros((num_batches, N_disc))
reaction_count = np.zeros((2, N_disc))
distance_traveled = np.zeros(N_disc)  # Для хранения пройденных расстояний в каждой ячейке
sources = [np.random.uniform(0, N_disc, neutrons)]
keff = 1.0


def neutron_tracking(Y, keff):
    direction = 2 * random.random() - 1
    life_span = -np.log(random.random()) / Sig_max  # Используем Sig_max для дельта-трекинга
    Y_append = Y + direction * life_span

    # Проверяем, было ли реальное взаимодействие
    if random.random() > Sig_tot / Sig_max:
        # Виртуальное взаимодействие - продолжаем движение
        if 0 <= Y_append < H:
            return neutron_tracking(Y_append, keff)
        else:
            return 0, None

    num_fiss = 0

    # Обработка границ
    if Y_append >= H:
        return num_fiss, None
    if Y_append <= 0:
        Y_append = -Y_append  # Отражение

    # Реальное взаимодействие
    cell_id = int(Y_append / h)
    start_cell = int(Y / h)

    # Вычисляем пройденное расстояние в каждой ячейке
    if start_cell == cell_id:
        distance_traveled[cell_id] += abs(Y_append - Y)
    else:
        if Y_append > Y:
            distance_traveled[start_cell] += (start_cell * h + h) - Y
            for i in range(start_cell + 1, cell_id):
                distance_traveled[i] += h
            distance_traveled[cell_id] += Y_append - cell_id * h
        else:
            distance_traveled[start_cell] += Y - start_cell * h
            for i in range(cell_id + 1, start_cell):
                distance_traveled[i] += h
            distance_traveled[cell_id] += (cell_id * h + h) - Y_append

    rand_inter = random.random()
    if rand_inter < Sig_s / Sig_tot:  # Рассеяние
        reaction_count[0, cell_id] += 1
        return neutron_tracking(Y_append, keff)

    if rand_inter < (Sig_s + Sig_cap) / Sig_tot:  # Поглощение
        reaction_count[1, cell_id] += 1
        return num_fiss, None

    # ДЕЛЕНИЕ
    fis_neutron = int(np.floor(Nu_f / keff))
    if random.random() < (Nu_f / keff) - fis_neutron:
        fis_neutron += 1
    num_fiss += fis_neutron
    return num_fiss, Y_append


for i in range(num_batches):
    generated = []
    distance_traveled[:] = 0  # Обнуляем пройденные расстояния для каждого батча

    for j in sources[-1]:
        created, position = neutron_tracking(j, keff)
        if position is not None:
            generated.extend([position] * created)

    # Расчет потока через пройденное расстояние (вместо счетчиков реакций)
    for k in range(N_disc):
        flux_arr[i, k] = distance_traveled[k] / (h * neutrons)

    keff = keff * len(generated) / neutrons
    neutrons = len(generated)
    sources.append(generated)
    reaction_count[:, :] = 0.

    print(f"Batch {i + 1}: keff = {keff:}, neutrons: {neutrons}")

flux_avg = np.mean(flux_arr[No_active:], axis=0)
y = np.arange(N_disc) * h

plt.figure(figsize=(8, 5))
plt.plot(y, flux_avg, color='green', label='Метод Монте-Карло с дельта-трекингом', linewidth=2)
plt.xlabel('Координаты')
plt.ylabel('Поток')
plt.title('Распределение потока нейтронов')
plt.legend()
plt.grid()
plt.show()"""
