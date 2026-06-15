# Runge-Kutta-Fehlberg 4(5) 自适应步长模板。  # 数值分析：用嵌套四阶与五阶解之差估计局部误差。
import numpy as np  # NumPy 用于动态结果数组的最终转换。
import matplotlib.pyplot as plt  # Matplotlib 用于画被接受步长的变化。

f = lambda t, y: -y + t**2 + 1  # 【必须替换】ODE 右端。
t0, T, y0 = 0.0, 1.0, 1.0  # 【必须替换】区间与初值。
h, h_min, h_max, tol = 0.1, 1e-6, 0.2, 1e-6  # 【必须替换】初始/最小/最大步长与每步局部误差限。
t_values, y_values, accepted_h = [t0], [y0], []  # Python 列表可动态追加，适合节点数量未知的自适应算法。
while t_values[-1] < T:  # Python 语法：[-1] 取最后元素；尚未到终点时继续。
    hn = min(h, h_max, T-t_values[-1])  # 保证不超过最大步长，并让最后一步恰好到 T。
    tn, yn = t_values[-1], y_values[-1]  # 取当前时间和当前数值解。
    K1 = f(tn, yn)  # Fehlberg 第 1 阶段斜率。
    K2 = f(tn+hn/4, yn+hn*K1/4)  # Fehlberg 第 2 阶段斜率。
    K3 = f(tn+3*hn/8, yn+hn*(3*K1+9*K2)/32)  # Fehlberg 第 3 阶段斜率。
    K4 = f(tn+12*hn/13, yn+hn*(1932*K1-7200*K2+7296*K3)/2197)  # Fehlberg 第 4 阶段斜率。
    K5 = f(tn+hn, yn+hn*(439*K1/216-8*K2+3680*K3/513-845*K4/4104))  # Fehlberg 第 5 阶段斜率。
    K6 = f(tn+hn/2, yn+hn*(-8*K1/27+2*K2-3544*K3/2565+1859*K4/4104-11*K5/40))  # Fehlberg 第 6 阶段斜率。
    y4 = yn+hn*(25*K1/216+1408*K3/2565+2197*K4/4104-K5/5)  # 嵌套四阶近似。
    y5 = yn+hn*(16*K1/135+6656*K3/12825+28561*K4/56430-9*K5/50+2*K6/55)  # 嵌套五阶近似。
    error = abs(y5-y4)  # 两个嵌套近似的差作为可计算的局部误差估计。
    if error <= tol:  # 只有误差合格时才接受本步。
        t_values.append(tn+hn)  # list.append 在列表末尾加入新时间点。
        y_values.append(y5)  # 接受精度更高的五阶近似。
        accepted_h.append(hn)  # 保存本次被接受的实际步长。
    factor = 4.0 if error == 0 else 0.9*(tol/error)**(1/5)  # 条件表达式处理零误差；五阶局部误差使用 1/5 次幂调步长。
    factor = min(4.0, max(0.1, factor))  # 限制单次步长变化倍数，防止步长剧烈振荡。
    h = min(h_max, max(h_min, hn*factor))  # 把下一候选步长限制到 [h_min,h_max]。
    assert error <= tol or h > h_min, "达到 h_min 仍不能满足误差限。"  # 无法继续时明确停止。
t = np.asarray(t_values, dtype=float)  # 把动态 Python 列表转换为 NumPy 浮点数组。
y = np.asarray(y_values, dtype=float)  # 转换数值解列表。
hs = np.asarray(accepted_h, dtype=float)  # 转换被接受步长列表。
for tn, yn in zip(t, y):  # 配对输出自适应节点与数值解。
    print(f"t={tn:.8f}, RKF45_y={yn:.10f}")  # 自适应节点通常不是整齐小数，因此多打印几位。
plt.step(t[:-1], hs, where="post", linewidth=1.3, label="accepted h")  # step 画分段常数步长变化。
plt.xlabel("t")  # 设置横轴名称。
plt.ylabel("h")  # 设置纵轴名称。
plt.grid(True)  # 打开网格。
plt.legend()  # 显示图例。
plt.tight_layout()  # 自动调整布局。
plt.show()  # 显示图像。
