from sympy import Symbol, solve
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt




# 기본 버전 (2023.03.15)
def calc_temp(insol, surface_temp, wall=True):
    x = Symbol('x')
    if wall == True:
        eq1 = (297-x)/0.5*1.6 + insol + 0.94*5.67*10**(-8)*(81-x**4)*0.5 + 0.94*5.67*10**(-8)*(surface_temp**4-x**4)*0.5
    else:
        eq1 = (297 - x) / 0.5 * 1.6 + insol + 0.94 * 5.67 * 10 ** (-8) * (81 - x ** 4)
    res = solve((eq1), dict=True)
    print(res)
    # print(res)
    # 결과 중에서 reasonable한 것을 선정
    reasonable_arr = []
    for item in res:
        try:
            if item[x] > 0:
                reasonable_arr.append(item[x])
        except:
            print('error occurred')
    if len(reasonable_arr) == 1:
        return reasonable_arr[0]
    elif len(reasonable_arr) == 0:
        # print('no result')
        return 0
    else:
        raise IOError

# TEG 부착 버전 (2023.03.23)
def calc_temp_teg(insol, surface_temp, wall=True):
    t_teg_out = Symbol('t_teg_out')
    t_teg_in = Symbol('t_teg_in')
    I_max = 11.3
    V_max = 27.3
    Q_max = 188.7
    dt_max = 78

    num_teg = 127             # N
    num_thermocouple = 216     # n
    seebeck_coeff = Q_max*(t_teg_out - dt_max)/num_thermocouple*t_teg_out**2*I_max       # alpha
    R_load = 0.02           # R_load
    R_teg = 1               # R (R of the TEM)
    K_teg = 1.61            # K (conductivity of the TEM)

    # delta_t = abs(t_teg_out - t_teg_in)
    if wall == True:
        # outside 가 hotside인 경우
        delta_t = t_teg_out - t_teg_in
        I_teg = seebeck_coeff * delta_t/(R_load + R_teg)
        eq1 = insol + 0.94*5.67*10**(-8)*(81-t_teg_out**4)*0.5 + 0.94*5.67*10**(-8)*(surface_temp**4-t_teg_out**4)*0.5 \
              - num_teg * num_thermocouple * (
                      seebeck_coeff * I_teg * t_teg_out
                      - 0.5 * I_teg **2 * R_teg
                      + K_teg * delta_t
              )
        eq2 = num_teg * num_thermocouple * (
                      seebeck_coeff * I_teg * t_teg_out
                      + 0.5 * I_teg **2 * R_teg
                      + K_teg * delta_t
              ) + (297 - t_teg_in) / 0.5 * 1.6
    else:
        eq1 = insol + 0.94 * 5.67 * 10 ** (-8) * (81 - t_teg_out ** 4)
    res = solve((eq1, eq2), dict=True)
    res = solve(eq1, dict=True)
    print(res)
    # print(res)
    # 결과 중에서 reasonable한 것을 선정
    # reasonable_arr = []
    # for item in res:
    #     try:
    #         if item[x] > 0:
    #             reasonable_arr.append(item[x])
    #     except:
    #         print('error occurred')
    # if len(reasonable_arr) == 1:
    #     return reasonable_arr[0]
    # elif len(reasonable_arr) == 0:
    #     # print('no result')
    #     return 0
    # else:
    #     raise IOError

calc_temp_teg(1360, 350, wall=True)
raise IOError
def calc_surface_temp(insol, alt, initial_temp):
    # emissivity = 0.7
    # albedo = 0.12
    # print(insol*(1-albedo)*math.cos(math.radians(alt)))
    # return (insol*(1-albedo)*math.cos(math.radians(alt))/(emissivity*5.6704*10**(-8)))**0.25

    # # Define physical constants
    # sigma = 5.67e-8  # Stefan-Boltzmann constant (W/m^2/K^4)
    #
    # # Define regolith properties
    # porosity = 0.5  # Porosity of the regolith
    # beta = 0.5  # Coefficient that takes into account the scattering of radiation within the regolith
    # rho = 1600.0  # Bulk density of the regolith (kg/m^3)
    # c = 680.0  # Specific heat capacity of the regolith (J/kg/K)
    # emissivity = 0.9  # Emissivity of the regolith
    #
    # # Define time steps and duration of lunar night
    # dt = 60.0  # Time step (s)
    # tmax = 14.0 * 24.0 * 3600.0  # Duration of lunar night (s)
    # nt = int(tmax / dt)  # Number of time steps
    #
    # # Define initial temperature and depth of regolith
    # T0 = 350.0  # Initial temperature of the regolith (K)
    # depth = 0.1  # Depth of the regolith that is being considered (m)
    #
    # # Calculate thermal conductivity and diffusivity of regolith
    # k = 2.0 * (1 - porosity) * (1 - 0.5 * beta) * sigma * T0 ** 3 / (3 * rho * c)
    # K = k / (rho * c)
    #
    # # Calculate time constant of regolith
    # tau = rho * c * depth ** 2 / (4 * k)
    #
    # # Calculate equilibrium temperature of regolith
    # q_rad = 0  # Assume no net radiative flux at start of lunar night
    # T_inf = (q_rad / (emissivity * sigma)) ** (1 / 4)
    #
    # # Initialize temperature array
    # T = np.zeros(nt)
    # T[0] = T0
    #
    # # Integrate temperature equation using forward Euler method
    # for i in range(1, nt):
    #     t = i * dt
    #     T[i] = T_inf + (T0 - T_inf) * np.exp(-t / tau)  # Temperature of regolith at time t
    #     q_abs = emissivity * sigma * T[i] ** 4  # Flux of radiation absorbed by regolith at time t
    #     q_em = emissivity * sigma * T[i - 1] ** 4  # Flux of radiation emitted by regolith at time t-dt
    #     q_rad = q_abs - q_em  # Net radiative flux at surface of regolith at time t
    #     T0 = T[i]  # Update initial temperature for next time step

    delta_temp = ((1-0.12)*insol*math.cos(math.radians(alt))-0.95*5.67*(10**(-8))*(initial_temp**4) + 0.016)/(50/3.6)*0.680
    print(initial_temp, delta_temp, (1-0.12)*insol*math.cos(math.radians(alt)), 0.9*5.67*10**(-8)*initial_temp**4)
    return initial_temp + delta_temp

# print(calc_surface_temp(0, 90))

data = pd.read_csv('./face_insol_moon.csv', encoding='utf-8 sig')
data = data[0:1400]
temp = 350
temp_arr = []
zen_arr = []
calc_temp_arr = []
for iter, row in data.iterrows():
    if row['alt_2'] >= 90:
        alt = 90
    else:
        alt = row['alt_2']
    temp = calc_surface_temp(row['normal_2'], alt, temp)
    temp_arr.append(temp)
    zen_arr.append(alt)
    calc_temp_arr.append(calc_temp(row['res_2'], temp))
    print(temp)

plt.plot(range(len(temp_arr)), temp_arr, label='temperature_moon')
plt.plot(range(len(temp_arr)), data['res_2']/3, label='insolation/3')
# plt.plot(range(len(temp_arr)), zen_arr)
plt.plot(range(len(temp_arr)), calc_temp_arr, label='temperature_wall')
plt.legend()
plt.show()