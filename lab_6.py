import numpy as np
import matplotlib.pyplot as plt
import random

Sig_s = 0.3
Sig_a = 0.03
Nu_f = 2.4
Nu_f_Sig_f = 0.03026952
Sig_f = Nu_f_Sig_f / Nu_f
Sig_tot = Sig_s + Sig_a
Sig_cap = Sig_a - Sig_f
H = 100
h = 5
N_disc = int(H / h)
w_min = 0.1
d = 10
keff_list = []

num_batches = int(input('Write a number of batches - queues: '))
neutrons = int(input('Write a number of neutrons for each batch: '))

No_active = int(0.5 * num_batches)

flux_arr = np.zeros((num_batches, N_disc))
reaction_count = np.zeros((2, N_disc))
sources = [np.random.uniform(0, H, neutrons)]
keff = 1.0

w_vac = 0.0
w_lost = 0.0
fission_weight_total = 0.0

def neutron_tracking(Y, w, keff):
    global w_vac, w_lost, fission_weight_total

    direction = 2 * random.random() - 1
    life_span = -np.log(random.random()) / Sig_tot
    Y_new = Y + direction * life_span

    if Y_new >= H:
        w_vac += w
        return []

    if Y_new <= 0:
        Y_new = -Y_new

    cell = int(Y_new / h)

    w_new = w * (Sig_s / Sig_tot)
    w_lost += w - w_new
    reaction_count[0, cell] += w

    w_fiss = w * (Nu_f * Sig_f / Sig_tot)
    fission_weight_total += w_fiss

    n_new = int(np.floor(w_fiss / keff))
    if random.random() < (w_fiss / keff - n_new):
        n_new += 1

    nt_new = [Y_new] * n_new

    if w_new >= w_min:
        weighted = neutron_tracking(Y_new, w_new, keff)
        nt_new.extend(weighted)
    else:
        if random.random() <= 1 / d:
            w_new *= d
            weighted = neutron_tracking(Y_new, w_new, keff)
            nt_new.extend(weighted)

    return nt_new

for i in range(num_batches):
    generated = []

    for nt_pos in sources[-1]:
        created = neutron_tracking(nt_pos, w=1.0, keff=keff)
        generated.extend(created)

    for k in range(N_disc):
        flux_arr[i, k] = reaction_count[0, k] / (Sig_s * h * len(sources[-1]))

    keff = keff * len(generated) / len(sources[-1])
    keff_list.append(keff)
    sources.append(generated)
    reaction_count[:, :] = 0.

    print(f"Batch {i+1}: keff = {keff:.5f}, neutrons: {len(generated)}")

avg_keff = np.mean(keff_list)
avg_keff_active = np.mean(keff_list[No_active:])
print(f"\n Средниий кефф  {avg_keff:.5f}")

flux_avg = np.mean(flux_arr[No_active:], axis=0)
x = np.arange(N_disc) * h

print(f"\n Взаимодействия  {w_lost:.3f}")
print(f"Вакуум {w_vac:.3f}")
print(f"Деление  {fission_weight_total:.3f}")

plt.figure(figsize=(9, 5))
plt.plot(x, flux_avg, color='darkblue', label='Average Flux', linewidth=2)
plt.xlabel('Coordinate')
plt.ylabel('Neutron Flux')
plt.title('Non-Analog Monte')
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()