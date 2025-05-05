
import math
import pandas as pd

# 读取参数 CSV 文件
U = pd.read_csv('params.csv')

# 提取 U_up 和 U_low 值
U_up = U['U_up'].iloc[217]
U_low = U['U_low'].iloc[217]

# 输出提取的参数
print(f"U_up: {U_up}, U_low: {U_low}")

# 参数设定
d = 0.00235  # 直径 (单位: 米)
l = 0.005  # 间距 (单位: 米)
A_side = 0.000036895  # 侧面积 (单位: 平方米)
R_h = 0.001175  # 半径 (单位: 米)
alpha = 0.0404  # 热导率常数
q = 1444  # 热源 (单位: W/m^2)
rho_C_f = 4200000  # 物质密度乘以比热容
rho_f = 1000  # 流体密度 (单位: kg/m^3)
C_p_f = 4200  # 流体比热 (单位: J/(kg·K))
k_f = 0.598  # 流体导热系数 (单位: W/(m·K))
epsilon = 0.588  # 孔隙度
rho_s = 552  # 固体密度 (单位: kg/m^3)
C_s = 1644  # 固体比热 (单位: J/(kg·K))
k_s = 0.2  # 固体导热系数 (单位: W/(m·K))
L_h = 0.005  # 长度 (单位: 米)


# 计算换热系数 h
def calculate_heat_transfer_coefficient(v):
    return 105.387 * math.sqrt(v)


# 计算面积（导热作用）
def calculate_area(d, L_h):
    return math.pi * d * L_h


# 计算体积（流体作用）
def calculate_volume(d, R_h, L_h, epsilon):
    return (d * l - math.pi * R_h ** 2) * L_h / epsilon


# 计算温度
def calculate_temperature(U_up, U_low, alpha):
    return (U_up + U_low) / alpha


# 计算总导热
def calculate_total_conduction(k, A_p, T, l):
    return k * A_p * T / l


# 计算系数 B
def calculate_B(A_p, h, k):
    return -A_p * h / k


# 计算 E
def calculate_E(B, L_h):
    return math.exp(B * L_h)


# 计算流体热
def calculate_fluid_heat(rho_C_f, v, T, l, V):
    return (rho_C_f) * v * T / l * V


# 计算体积热源
def calculate_volume_heat_source(q, L_h, k, rho_C_f, l, v):
    return q * L_h * k / (rho_C_f) * l * v


# 计算导热
def calculate_conduction(T, A_side, A_p, h, L_h, E, q, l):
    return T * A_side / (A_p * h * L_h ** 2) * E + (A_p * q) / (2 * l)


# 主计算逻辑
def calculate_v(U_up, U_low):
    # 假设初始值
    v = 0.00002  # 初始假设值
    h = calculate_heat_transfer_coefficient(v)

    # 计算面积
    A_p = calculate_area(d, L_h)

    # 计算温度
    T = calculate_temperature(U_up, U_low, alpha)

    # 计算总导热
    P = calculate_total_conduction(k_s, A_p, T, l)

    # 计算系数 B
    B = calculate_B(A_p, h, k_s)

    # 计算 E
    E = calculate_E(B, L_h)

    # 计算体积
    V = calculate_volume(d, R_h, L_h, epsilon)

    # 计算流体热
    P_f = calculate_fluid_heat(rho_C_f, v, T, l, V)

    # 计算体积热源
    q_prime = calculate_volume_heat_source(q, L_h, k_s, rho_C_f, l, v)

    # 计算导热
    P_s = calculate_conduction(T, A_side, A_p, h, L_h, E, q, l)

    # 总方程
    P = P_f + P_s

    return v, P


# 计算 v 和总功率
v_calculated, total_power = calculate_v(U_up, U_low)

# 输出结果
print(f"计算得到的未知量 v: {v_calculated:.6f} m/s")
print(f"计算得到的总功率 P: {total_power:.6f} W")
