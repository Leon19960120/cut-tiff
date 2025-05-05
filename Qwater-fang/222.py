import math
import pandas as pd

# 读取参数 CSV 文件
U = pd.read_csv('params.csv')

# 提取 U_up 和 U_low 值
U_up = U['U_up'].iloc[0]
U_low = U['U_low'].iloc[0]

# 输出提取的参数
print(f"U_up: {U_up}, U_low: {U_low}")

# 参数设定
d = 0.00235  # 直径 (单位: 米)
l = 0.005  # 间距 (单位: 米)
A_side = 0.000036895  # 侧面积 (单位: 平方米)
R_h = 0.001175  # 半径 (单位: 米)
alpha = 0.0404  # 热导率常数
q = 1444  # 热源 (单位: W/m^2)
rhoC_f = 4200000  # 水热容
rho_f = 1000  # 流体密度 (单位: kg/m^3)
C_f = 4200  # 流体比热 (单位: J/(kg·K))
k_f = 0.598  # 流体导热系数 (单位: W/(m·K))
epsilon = 0.588  # 孔隙度
rho_s = 552  # 固体密度 (单位: kg/m^3)
C_s = 1644  # 固体比热 (单位: J/(kg·K))
k_s = 0.2  # 固体导热系数 (单位: W/(m·K))
L_h = 0.005  # 长度 (单位: 米)

# 计算有效导热系数 k
def calculate_k(k_s, epsilon, k_f):
    return k_s ** (1 - epsilon) * k_f ** epsilon

# 计算有效热容 RhoC_eff
def calculate_rhoC_eff(epsilon, rho_s, C_s, rho_f, C_f):
    return (1 - epsilon) * rho_s * C_s + epsilon * rho_f * C_f

# 计算换热系数 h
def calculate_h(v):
    return 105.387 * math.sqrt(v)

# 计算面积（导热作用A_p）
def calculate_A_p(d, L_h):
    return math.pi * d * L_h

# 计算体积（流体作用）
def calculate_V(d, l, R_h, L_h, epsilon):
    return (d * l - math.pi * R_h ** 2) * L_h / epsilon

# 计算温度
def calculate_T(U_up, U_low, alpha):
    return (U_up + U_low) / alpha

# 计算总导热
def calculate_P(k, A_p, T, l):
    return k * A_p * T / l

# 计算系数 B
def calculate_B(A_p, h, k):
    return -A_p * h / k

# 计算 E
def calculate_E(B, L_h):
    return math.exp(B / L_h)

# 计算流体热
def calculate_P_f(rhoC_f, v, T, l, V):
    return rhoC_f * v * T / l * V

# 计算体积热源
def calculate_q_triple(q, L_h, k, rhoC_eff, l, v):
    return q / L_h * k / rhoC_eff / l * v

# 计算导热
def calculate_P_s(T, A_side, A_p, h, L_h, E, q_triple, l):
    return T * A_side * A_p * h / (L_h ** 2) * E + A_p * q_triple * l / 2

# 主计算逻辑
def calculate_v(U_up, U_low):
    # 初始假设值
    v = 0.00002  # 初始假设值

    # 计算中间参数
    k = calculate_k(k_s, epsilon, k_f)
    rhoC_eff = calculate_rhoC_eff(epsilon, rho_s, C_s, rho_f, C_f)
    h = calculate_h(v)
    A_p = calculate_A_p(d, L_h)
    T = calculate_T(U_up, U_low, alpha)
    P = calculate_P(k, A_p, T, l)
    B = calculate_B(A_p, h, k)
    E = calculate_E(B, L_h)
    V = calculate_V(d, l, R_h, L_h, epsilon)
    P_f = calculate_P_f(rhoC_f, v, T, l, V)
    q_triple = calculate_q_triple(q, L_h, k, rhoC_eff, l, v)
    P_s = calculate_P_s(T, A_side, A_p, h, L_h, E, q_triple, l)

    # 总方程
    total_power = P_f + P_s

    return v, total_power

# 计算 v 和总功率
v_calculated, total_power = calculate_v(U_up, U_low)

# 输出结果
print(f"计算得到的未知量 v: {v_calculated:.6f} m/s")
print(f"计算得到的总功率 P: {total_power:.6f} W")