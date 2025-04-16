import numpy as np
import matplotlib.pyplot as plt

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

plt.plot(x, phi, label='Асимптотический поток', color='green')
plt.xlabel('Положение (см)')
plt.ylabel('Поток нейтронов')
plt.title('Распределение потока нейтронов')
plt.legend()
plt.grid()
plt.show()

'''delta_lambdas = np.linspace(0, 0.01, 10)
for delta_lambda in delta_lambdas:
    keff, phi, iterations = accelerated_source_iteration(A, nuSigma_f, delta_lambda)
    print(f"delta_lambda = {delta_lambda:.2f}, k_eff = {keff:.6f}, сходимость на {iterations}-й итерации")'''