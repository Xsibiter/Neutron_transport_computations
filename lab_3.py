import numpy as np
import matplotlib.pyplot as plt

# Первая часть кода (ваш оригинальный код)
sigma_a = 0.03
sigma_gen = sigma_a
H = 100
h = 1
D = 1
N = int(H / h)

Q = np.zeros(N)
Flow = np.ones(N)
J_input = np.zeros((N, 2))
J_output = np.zeros((N, 2))
Flow_last = np.zeros(N)

k_eff = 1.0
k = 1.0
count = 0
tol = 1e-10

D_tilde = D / h
A = (6 * D_tilde * (1 + 4 * D_tilde)) / (1 + 16 * D_tilde + 48 * (D_tilde ** 2))
B = (1 - 48 * (D_tilde ** 2)) / (1 + 16 * D_tilde + 48 * (D_tilde ** 2))
C = (-8 * D_tilde) / (1 + 16 * D_tilde + 48 * (D_tilde ** 2))

while not np.abs(k - k_eff) < tol or np.max(np.abs(Flow_last - Flow)) > tol:
    count += 1
    Last_Norm = np.linalg.norm(Flow)
    k = k_eff

    Q[:] = sigma_gen * Flow[:] / k_eff
    Flow_last[:] = Flow[:]

    for i in range(N):
        Flow[i] = (Q[i] + (1 - B - C) * (J_input[i][1] + J_input[i][0])) / (2 * A + sigma_a)

        J_output[i, 1] = B * J_input[i][1] + C * J_input[i][0] + A * Flow[i]
        J_output[i, 0] = C * J_input[i][1] + B * J_input[i][0] + A * Flow[i]

    for i in range(1, N - 1):
        J_input[i][0] = J_output[i - 1][1]  # Левая граница текущей ячейки = выход правой соседа
        J_input[i][1] = J_output[i + 1][0]  # Правая граница текущей ячейки = вход левой соседа

    # Граничные условия
    J_input[0][0] = J_output[0][0]  # Левая граница (условие отражения)
    J_input[0][1] = J_output[1][0]
    J_input[N - 1][1] = 0  # Правая граница (абсорбция)
    J_input[N - 1][0] = J_output[N - 2][1]

    k_eff = k * (np.linalg.norm(Flow) / Last_Norm)

print(f"Эффективный коэффициент размножения k_eff: {k_eff:.6f}")


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
A[-1, -1] = D / dx ** 2 + Sigma_a + 1 / 2
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


Flow_norm = Flow / np.linalg.norm(Flow)  # Нодальный
phi_norm = phi / np.linalg.norm(phi)  # Ускоренный


error = np.abs(Flow_norm - phi_norm)

# Построение графиков
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), gridspec_kw={'height_ratios': [1, 1]})

# График потоков
x_flow = np.linspace(h / 2, H - h / 2, N)
ax1.plot(x_flow, Flow_norm, color='red', linewidth=1, label="Нодальный метод")
ax1.plot(x, phi_norm, color='yellow', linewidth=1, label="Ускоренный метод")
ax1.set_xlabel('Координата')
ax1.set_ylabel('Нормированный поток нейтронов')
ax1.legend()
ax1.grid()

ax2.plot(x_flow, error, color='red', linewidth=1, label="Погрешность (|Нодальный - Ускоренный|)")
ax2.set_xlabel('Координата')
ax2.set_ylabel('Погрешность')
ax2.legend()
ax2.grid()

plt.tight_layout()
plt.show()