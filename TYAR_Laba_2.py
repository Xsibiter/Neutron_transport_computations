import os
import subprocess
import re
import pandas as pd

# 📌 Концентрации для различных замедлителей:
moderators = {
    "H2O": {
        "U5": [1.44997E-05, 7.29858E-06, 4.87665E-06, 3.66160E-06, 1.83390E-06, 1.22329E-06, 9.17723E-07, 7.34303E-07],
        "U8": [6.43257E-04, 3.23792E-04, 2.16346E-04, 1.62442E-04, 8.13583E-05, 5.42694E-05, 4.07135E-05, 3.25763E-05],
        "H": [6.59075E-02, 6.63508E-02, 6.64998E-02, 6.65746E-02, 6.66871E-02, 6.67247E-02, 6.67435E-02, 6.67548E-02],
        "O": [3.29538E-02, 3.31754E-02, 3.32499E-02, 3.32873E-02, 3.33436E-02, 3.33624E-02, 3.33718E-02, 3.33774E-02]
    }
}

# 📌 Функция изменения lab2.txt
def generate_input(moderator, idx):
    with open("lab2.txt", "r") as template:
        data = template.read()

    # Заменяем концентрации в файле
    for nuclide, values in moderator.items():
        data = data.replace(f"@ {nuclide}  @", f"@ {nuclide}  @   {values[idx]:.5E},")

    input_filename = f"lab2_{idx}.txt"
    with open(input_filename, "w") as new_input:
        new_input.write(data)

    return input_filename


# 📌 Функция запуска GETERA-93
def run_getera(input_filename):
    subprocess.run(["getera.exe", input_filename], check=True)


# 📌 Функция для парсинга макро- и микросечений
def parse_output():
    with open("lab2.out", "r") as output_file:
        lines = output_file.readlines()

    # Извлекаем макросечения (Σa1, Σa2, Σs 1-->2)
    sabs_values = [None, None]
    sigma_s1_2 = None

    for i, line in enumerate(lines):
        if "i / j -->" in line:
            row_1 = re.findall(r"[\d\.E+-]+", lines[i + 1])  # Строка 1 (Σa1)
            row_2 = re.findall(r"[\d\.E+-]+", lines[i + 2])  # Строка 2 (Σa2)
            sabs_values[0] = float(row_1[1])  # Σa1
            sabs_values[1] = float(row_2[1])  # Σa2
            sigma_s1_2 = float(row_1[2])  # Σs 1-->2

    # Извлекаем микросечения (capture для U-235 и U-238)
    sigma_values = []
    for i in range(154, 178):
        if "u235" in lines[i] or "u238" in lines[i]:
            values = re.findall(r"[\d\.E+-]+", lines[i])
            sigma_values.append(float(values[5]))  # capture-сечение

    return sigma_values, sabs_values, sigma_s1_2  # Возвращаем нужные параметры


# 📌 Запуск всех расчетов
results = []

for moderator_name, data in moderators.items():
    for idx in range(8):  # Прогоняем все случаи концентраций
        input_fQile = generate_input(data, idx)
        run_getera(input_file)

        sigma, sa_values, s12 = parse_output()

        row = [sigma[0], sigma[1], sigma[2], sigma[3], sa_values[0], sa_values[1], s12]
        results.append(row)

# 📌 Сохранение данных в таблицу
columns = ["σa1(5)", "σa2(5)", "σa1(8)", "σa2(8)", "Σa1", "Σa2", "Σs 1--->2"]
df = pd.DataFrame(results, columns=columns)
df.to_excel("results.xlsx", index=False)

print("✅ Готово! Данные сохранены в `results.xlsx`")