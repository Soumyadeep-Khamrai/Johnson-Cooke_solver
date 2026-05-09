#This is a material model solver for metals under high strain rate at elevated temp

'''
σ=(A+Bε^n )(1+Clnε*)(1-〖T*〗^m ) - Fornula
σ is the equivalent flow stress
ε is the equivalent plastic strain
A - yield stress of the material under reference conditions
B - strain hardening const
n - strain hardening coeff
ε* - dimensionless strain rate, ε * = ε./ ε.ref
T* - homologous temp, T* = (T - Tref)/(Tm – Tref), Tm - Melting temp
'''

from cmath import log
import numpy as np
import matplotlib.pyplot as plt
import csv

# Material properties input

young_modulus = float(input("Enter the Young’s Modulus value in MPa: "))
A = float(input("Enter the yield strength: "))
B = float(input("Enter the value of B: "))
n = float(input("Enter the value of n: "))
C = float(input("Enter the value of C: "))
m = float(input("Enter the value of m: "))
e_ref = float(input("Enter the reference strain rate: "))
t_ref = float(input("Enter the reference temperature (Kelvin): "))
t_melt = float(input("Enter the melting temperature (Kelvin): "))

# Test conditions

e = float(input("Enter the strain rate: "))
t = float(input("Enter the test temperature (Kelvin): "))

e_dimless = e / e_ref
t_dimless = (t - t_ref) / (t_melt - t_ref)

strain_list = []
stress_list = []

stress = 0
stress_increment = 48.9

while stress < A:
    strain = stress / young_modulus

    strain_list.append(np.real(strain))
    stress_list.append(np.real(stress))

    print(strain, stress)

    stress += stress_increment

# Plastic region (Johnson-Cook)

for strain_plastic in range(0, 100):

    stress = (
        (A + B * (strain_plastic / 1000) ** n)
        * (1 + C * log(e_dimless))
        * (1 - t_dimless ** m)
    )

    strain_elastic = stress / young_modulus
    strain = strain_elastic + (strain_plastic / 1000)

    strain_list.append(np.real(strain))
    stress_list.append(np.real(stress))

    print(np.real(strain), np.real(stress))

#CSV
csv_filename = f"JC_curve_SR_{e}_T_{t}K.csv"
with open(csv_filename, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Strain Rate", e])
    writer.writerow(["Temperature (K)", t])
    writer.writerow([])
    writer.writerow(["Strain", "Stress (MPa)"])
    for s1, s2 in zip(strain_list, stress_list):
        writer.writerow([s1, s2])
print(f"\nData written to: {csv_filename}")

#PLOTTED DATA
plt.figure(figsize=(8, 6))
plt.plot(
    strain_list,
    stress_list,
    linewidth=2,
    label=f"SR={e}, T={t}K"
)
plt.xlabel("Strain")
plt.ylabel("Stress (MPa)")
plt.title("Johnson-Cook Stress-Strain Curve")
plt.grid(True)
plt.legend()
plt.show()