from sympy import Symbol, solve
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def calc_surface_temp(insol, alt, initial_temp):
    delta_temp = ((1-0.12)*insol*math.cos(math.radians(alt))-0.95*5.67*(10**(-8))*(initial_temp**4) + 0.016)/(50/3.6)*0.680
    # print(initial_temp, delta_temp, (1-0.12)*insol*math.cos(math.radians(alt)), 0.9*5.67*10**(-8)*initial_temp**4)
    return initial_temp + delta_temp

class object:
    def __init__(self, initial_temp, rho, Cp, k, Area, thickness, edge_cond):
        self.temp = initial_temp
        self.rho = rho
        self.Cp = Cp
        self.k = k
        self.Area = Area
        self.thickness = thickness
        self.volume = Area * thickness
        self.edge_cond = edge_cond

        if edge_cond == True:
            self.volume = self.volume * 0.5

        self.q_in = []
        self.q_out = []

    def update_temp(self):
        # in/out 어레이의 값을 가지고 비체적을 통해 온도변화를 반영한 뒤, 어레이 값을 초기화
        self.delta_q = sum(self.q_in) - sum(self.q_out)
        self.temp = self.temp + self.delta_q /(self.rho * self.Cp * self.volume)
        self.q_in = []
        self.q_out = []

# class t_con:
#     def __init__(self, obj_a, obj_b):
#         self.object_a = obj_a
#         self.object_b = obj_b
# 
# 
# class t_con_conv(t_con):
#     def __init__(self, obj_a, obj_b):
#         super(t_con, self).__init__(obj_a, obj_b)
# 
# class t_con_cond(t_con):
#     def __init__(self, obj_a, obj_b):
#         super(t_con, self).__init__(obj_a, obj_b)
#     def process(self):
#         if self.obj_a.temp > self.obj_b.temp:
#             self.obj_a.q_out.append(self.obj_a.k * (self.obj_a - self.obj_b) / ((self.obj_a.thickness + self.obj_b.thickness)*0.5))
#             self.obj_b.q_in.append(self.obj_b.k * (self.obj_a - self.obj_b) / ((self.obj_a.thickness + self.obj_b.thickness) * 0.5))
#         elif self.obj_b.temp > self.obj_a.temp:
#             self.obj_b.q_out.append(self.obj_b.k * (self.obj_b - self.obj_a) / ((self.obj_a.thickness + self.obj_b.thickness) * 0.5))
#             self.obj_a.q_in.append(self.obj_a.k * (self.obj_b - self.obj_a) / ((self.obj_a.thickness + self.obj_b.thickness) * 0.5))
# class t_con_rad(t_con):
#     def __init__(self, obj_a, obj_b):
#         super(t_con, self).__init__(obj_a, obj_b)

def process_cond(obj_a, obj_b):
    if obj_a.temp > obj_b.temp:
        obj_a.q_out.append(
            obj_a.k * (obj_a.temp - obj_b.temp) / ((obj_a.thickness + obj_b.thickness) * 0.5) * obj_a.Area)
        obj_b.q_in.append(
            obj_b.k * (obj_a.temp - obj_b.temp) / ((obj_a.thickness + obj_b.thickness) * 0.5) * obj_a.Area)
    elif obj_b.temp > obj_a.temp:
        obj_b.q_out.append(
            obj_b.k * (obj_b.temp - obj_a.temp) / ((obj_a.thickness + obj_b.thickness) * 0.5) * obj_b.Area)
        obj_a.q_in.append(
            obj_a.k * (obj_b.temp - obj_a.temp) / ((obj_a.thickness + obj_b.thickness) * 0.5) * obj_b.Area)
    else:
        pass

T_teg_arr = [[250, 250, 250]]
T_conc_arr = [[250, 250, 250]]

# TEG Characteristics
R_load = 0.02
A_TEM=0.5*0.5
f_TEM = 0.5                     # The fraction of A covered by the TE elements of the TEC
L_TEM=0.003                     # m %TE element length or height
count_thermocouple = 127        # Number of couples in the TEC module
kappa_teg = 32      # W/mK
rho_teg = 3950      # kg/m^3
cp_teg = 800        # J/kgK
I_max = 11.3    # A
V_max = 24.6    # V
dT_max = 200     # When I_max, Q = 0, dT_max %C
Q_max = 172.0   # When I_max, Q_max, dT = 0 %W
delta_x_teg = L_TEM/(len(T_teg_arr[-1])+1)

# Concrete Wall Characteristics
L_CONC = 0.5
kappa_conc = 1.6
rho_conc = 2500
cp_conc = 880
delta_x_conc = L_CONC/(len(T_conc_arr[-1])+1)

# Regolith Characteristics
initial_temp = 300

total_length = 1680       # Hours

M = 360*total_length         # Number of time step
t = 3600*total_length            # s %Simulatino duration %seconds
timestep = t/(M)       # time step duration(s)



#data = pd.read_csv('./face_insol_moon.csv', encoding='utf-8 sig')
data = pd.read_excel('./year_result_wall(east).xlsx', engine='openpyxl')
data = data[0+140:total_length+1+140]
data = data[['normal_2', 'res_2', 'alt_2']]
data.index = data.index*3600
data.index = pd.to_datetime(data.index, unit='s')
data_resampled = data.resample(rule=str(timestep) + 'S').last().interpolate()

# 중요!!! kappa하고 K값이 제대로 구분되어 고려하였는지 확인해 볼 것! (한솔선배 논문 9-10페이지)


insol_arr = data_resampled['normal_2'].values.tolist()
tilted_insol_arr = data_resampled['res_2'].values.tolist()
alt_arr = data_resampled['alt_2'].values.tolist()

df_dict = {
    'tilted_insol':tilted_insol_arr[0:-1],
    'altitude':alt_arr[0:-1],
    'normal_insol':insol_arr[0:-1],
    'surface_temp':[]
}

for i in range(len(T_teg_arr[0])):
    df_dict['TEG_' + str(i+1)] = []
for i in range(len(T_conc_arr[0])):
    df_dict['CONC_' + str(i+1)] = []

for k in range(M):
    T_teg = T_teg_arr[-1]
    T_conc = T_conc_arr[-1]

    # TEG 계산
    dt_teg = T_teg[0] - T_teg[-1]
    T_h = T_teg[0]
    # T_c = min(T_teg[0], T_teg[-1])
    rho = 2*A_TEM*f_TEM*((T_h-dT_max)**2)*Q_max/(2*(T_h**2)*L_TEM*(count_thermocouple**2)*(I_max**2))
    alpha = 2*Q_max*(T_h-dT_max)/(count_thermocouple*(T_h**2)*I_max)
    R = 2 * (count_thermocouple**2) * L_TEM * rho / (A_TEM * f_TEM) / count_thermocouple

    I_teg = alpha * (dt_teg) / (R_load + R)
    V_teg = count_thermocouple * alpha * dt_teg / ((R_load/R + 1)*R_load/R)
    P_dc_teg = I_teg * V_teg
    kappa_teg = ((T_h-dT_max)**2)*Q_max/((T_h**2)*dT_max)

    T_next_teg = T_teg
    T_next_conc = T_conc

    # 복사계산을 위한 전처리
    # loop수에 맞추어 insolation, altitude값을 로드
    insol = insol_arr[k]
    tilted_insol = tilted_insol_arr[k]
    alt = alt_arr[k]
    surface_temp = calc_surface_temp(insol, alt, initial_temp)


    # 발전량에 따라 변하는 열전소자 양단 끝쪽 처리
    if P_dc_teg >= 0:
        T_next_teg[0] = T_teg[0] + (timestep/delta_x_teg**2 * kappa_teg / rho_teg / cp_teg * (T_teg[1] - T_teg[0])
                             + (0.94*5.67*10**(-8)*(81-T_teg[0]**4)*0.5)*(timestep/3600) + (0.94*5.67*10**(-8)*(surface_temp**4-T_teg[0]**4)*0.5)*(timestep/3600) + (tilted_insol)*(timestep/3600)) # insolation cond
        T_next_teg[-1] = T_teg[-1] + (timestep / delta_x_teg ** 2 * kappa_teg / rho_teg / cp_teg * (T_teg[-2] - T_teg[-1])
                               + timestep / delta_x_conc ** 2 * kappa_conc / rho_conc / cp_conc * (T_conc[1] - T_teg[-1])
                               - P_dc_teg)
    else:
        T_next_teg[0] = T_teg[0] + (timestep/delta_x_teg**2 * kappa_teg / rho_teg / cp_teg * (T_teg[1] - T_teg[0])
                             + 0.94*5.67*10**(-8)*(81-T_teg[0]**4)*0.5 + 0.94*5.67*10**(-8)*(surface_temp**4-T_teg[0]**4)*0.5 + tilted_insol
                             + P_dc_teg)
        T_next_teg[-1] = T_teg[-1] + (timestep / delta_x_teg ** 2 * kappa_teg / rho_teg / cp_teg * (T_teg[-2] - T_teg[-1])
                               + timestep / delta_x_conc ** 2 * kappa_conc / rho_conc / cp_conc * (T_conc[1] - T_teg[-1]))

    # 열전소자 나머지 셀 처리
    for i in range(len(T_teg)-2):
        T_teg[i+1] = T_teg[i+1] + (timestep / delta_x_teg ** 2 * kappa_teg / rho_teg / cp_teg * (T_teg[i] - 2*T_teg[i+1] + T_teg[i+2]))

    # 열전소자-콘크리트 벽 연결부 처리 (T_conc[0]으로 놓고 하지만, 이게 T_teg[-1]이랑 같은 거 아닌지...?)
    T_next_conc[0] = T_conc[0] + (timestep / delta_x_conc ** 2 * kappa_conc / rho_conc / cp_conc * (T_conc[1] - T_conc[0])
                               + timestep / delta_x_teg ** 2 * kappa_teg / rho_teg / cp_teg * (T_teg[-2] - T_conc[0]))

    # 콘크리트 벽 내부면 처리 (온도를 24도로 놓고 고정)
    T_next_conc[-1] = 24+273

    # 콘크리트 나머지 셀 처리
    for i in range(len(T_conc)-2):
        T_conc[i + 1] = T_conc[i + 1] + (timestep / delta_x_conc ** 2 * kappa_conc / rho_conc / cp_conc * (T_conc[i] - 2 * T_conc[i + 1] + T_conc[i + 2]))

    # print(T_teg)
    # print(T_next_teg)
    for i in range(len(T_teg)):
        df_dict['TEG_' + str(i + 1)].append(T_teg[i])
    for i in range(len(T_next_teg)):
        df_dict['CONC_' + str(i + 1)].append(T_conc[i])
    df_dict['surface_temp'].append(surface_temp)
    T_teg_arr.append(T_next_teg)
    T_conc_arr.append(T_next_conc)




# df_dict = {
#     'tilted_insol':tilted_insol_arr,
#     'altitude':alt_arr,
#     'normal_insol':insol_arr
# }
#
# for i in range(len(T_teg_arr[0])):
#     df_dict['TEG_' + str(i+1)] = []
# for i in range(len(T_conc_arr[0])):
#     df_dict['CONC_' + str(i+1)] = []
#
# for i in range(len(T_teg_arr)):
#     for j in range(len(T_teg_arr[0])):
#         df_dict['TEG_' + str(j+1)].append(T_teg_arr[i][j])
#     for j in range(len(T_conc_arr[0])):
#         df_dict['CONC_' + str(j+1)].append(T_conc_arr[i][j])

for key in df_dict.keys():
    print(len(df_dict[key]))

df = pd.DataFrame(df_dict)
print(df)
# 초 단위 데이터 저장
# df.to_csv('./result.csv', encoding='ANSI')

# 시간 단위 데이터 저장
df.index = pd.to_datetime(df.index, unit='s')
df = df.resample(rule='H').mean().interpolate()
df.to_csv('./result.csv', encoding='ANSI')


# 일사량 비교를 위한 Scaling
df['tilted_insol'] = df['tilted_insol']/3

plt.plot(df.index, df['tilted_insol'], label='Insolation(1/3)')
plt.plot(df.index, df['surface_temp'], label='Surface_temp')
plt.plot(df.index, df['TEG_1'], label='TEG_1')
plt.plot(df.index, df['CONC_1'], label='CONC_1')
plt.plot(df.index, df['CONC_2'], label='CONC_2')
plt.plot(df.index, df['CONC_3'], label='CONC_3')
plt.legend()
plt.show()