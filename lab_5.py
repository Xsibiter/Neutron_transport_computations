import numpy as np
import matplotlib.pyplot as plt
import random

H = 100
h = 2.5
N_disc = int(H / h)

materials = [  # ПОГЛОТИЛЬЕ ДЕЛЕНИЕ РАССЕНИЕ
    {"Σ_s": 0.08, "Σ_a": 0.2, "ν_f": 0, "Σ_f": 0},  # 0 - 10
    {"Σ_s": 0.3, "Σ_a": 0.03, "ν_f": 2.5, "Σ_f": 0.015},  # 10-80
    {"Σ_s": 0.8, "Σ_a": 0.01, "ν_f": 0, "Σ_f": 0},  # 80-100
]

boundaries = [0,2, 80,100]
total_neutrons = 0
total_inters = 0
cell_inter = np.zeros(N_disc)

def get_material(y):
    for i in range(len(boundaries) - 1):
        if boundaries[i] <= y < boundaries[i + 1]:
            return materials[i]
    return materials[-1]
def get_next_boundary(Y, direction):
    if direction > 0:
        return min([b for b in boundaries if b > Y], default=H)
    else:
        return max([b for b in boundaries if b < Y], default=0)

num_batches = int(input('Write a number of batches - queues: '))
neutrons = int(input('Write a number of neutrons for each batch: '))

No_active = int(0.5 * num_batches)
flux_arr = np.zeros((num_batches, N_disc))
reaction_count = np.zeros((2, N_disc))
sources = [np.random.uniform(0, H, neutrons)]
keff = 1.0
b_10 = 0
b_80 = 0

def neutron_tracking(Y, keff):
    global b_10, b_80
    num_fiss = 0
    prev_Y = Y

    while True:
        material = get_material(Y)
        Σ_tot = material["Σ_s"] + material["Σ_a"]
        direction = 2 * random.random() - 1
        life_span = -np.log(random.random()) / Σ_tot
        Y_new = Y + direction * life_span

        Y_new = max(0, min(Y_new, H - 1e-6))

        next_boundary = get_next_boundary(Y, direction)

        if Y_new >= H:
            return num_fiss, None
        if Y_new <= 0:
            Y_new = -Y_new  # Отражение

        if (direction > 0 and Y_new >= next_boundary) or (direction < 0 and Y_new <= next_boundary):
            Y = next_boundary
            continue

        if (prev_Y < 10 and Y_new >= 10) or (prev_Y > 10 and Y_new <= 10):
                b_10 += 1


        if (prev_Y < 80 and Y_new >= 80) or (prev_Y > 80 and Y_new <= 80):
                b_80 += 1

        rand_inter = random.random()
        cell_id = min(int(Y_new / h), N_disc - 1)
        cell_inter[cell_id] += 1

        if rand_inter < material["Σ_s"] / Σ_tot:
            reaction_count[0, cell_id] += 1 / material["Σ_s"]
            direction = 2 * random.random() - 1
            Y_new = Y + direction * life_span
            Y = Y_new
            continue

        elif rand_inter < (material["Σ_s"] + material["Σ_a"] - material["Σ_f"]) / Σ_tot:
            reaction_count[1, cell_id] += 1
            return num_fiss, None
            print()
        fis_neutron = int(np.floor(material["ν_f"] / keff))
        if random.random() < (material["ν_f"] / keff) - fis_neutron:
            direction = 2 * random.random() - 1
            Y_new = Y + direction * life_span
            fis_neutron += 1

        num_fiss += fis_neutron
        return num_fiss, Y_new


for i in range(num_batches):
    generated = []
    for j in sources[-1]:
        created, position = neutron_tracking(j, keff)
        if position is not None:
            generated.extend([position] * created)

    for k in range(N_disc):
        flux_arr[i, k] = reaction_count[0, k] / (h * neutrons)

    total_neutrons += neutrons

    keff = keff * len(generated) / neutrons
    neutrons = len(generated)
    sources.append(generated)
    reaction_count[:, :] = 0


    print(f"Batch {i + 1}: keff = {keff:.5f}, neutrons: {neutrons}")

total_inters = np.sum(cell_inter)
print(f"Пересечения границы Y = 10: {b_10}")
print(f"Пересечения границы Y = 80: {b_80}")
print(f"Общее количество нейтронов за все циклы: {total_neutrons}")
print(f"Общее число пересечений: {total_inters}")
flux_avg = np.mean(flux_arr[No_active:], axis=0)

y = np.arange(N_disc) * h

fig, ax1 = plt.subplots()

ax1.plot(y, flux_avg, color='green', label="Нейтронный поток")
ax1.set_xlabel('Координата (Y)')
ax1.set_ylabel('Нейтронный поток', color='green')
ax1.tick_params(axis='y', labelcolor='green')

ax1.axvline(x=10, color='red', linestyle='dashed', label="Граница Y=10")
ax1.axvline(x=80, color='blue', linestyle='dashed', label="Граница Y=80")

ax2 = ax1.twinx()
ax2.set_ylabel('Число пересечений границ', color='black')
ax2.scatter([10], [b_10], color="red", s=100, label=f"Пересечения Y=10: {b_10}")
ax2.scatter([80], [b_80], color="blue", s=100, label=f"Пересечения Y=80: {b_80}")


plt.grid()
plt.show()

