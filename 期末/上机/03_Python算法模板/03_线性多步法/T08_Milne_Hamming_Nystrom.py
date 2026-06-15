# Milne、Hamming、Nyström 多步法模板。  # 数值分析：这些方法需要额外起步值，并可能存在寄生根振荡。
import math  # 标准库 math 用于默认例题的指数函数。
import numpy as np  # NumPy 用于固定网格和数值解数组。
import matplotlib.pyplot as plt  # Matplotlib 用于绘图。

f = lambda t, y: -y + math.exp(-t)  # 【必须替换】ODE 右端。
t0, T, y0, h = 0.0, 1.0, 0.0, 0.1  # 【必须替换】区间、初值、步长。
method = "Milne"  # 【必须替换】可选 Milne/Hamming/Nystrom。
N = round((T-t0)/h)  # 计算固定步长总步数。
assert abs(t0+N*h-T) < 1e-12, "步长必须整除区间长度。"  # 检查网格终点。
t = t0 + np.arange(N+1)*h  # 创建等距时间网格。
y = np.zeros(N+1, dtype=float)  # 预分配数值解。
y[0] = y0  # 写入初值。
start_count = 1 if method == "Nystrom" else 3  # Nyström 需两个点；Milne/Hamming 需四个点。
for n in range(start_count):  # 用 RK4 计算所需起步值。
    k1 = f(t[n], y[n])  # RK4 第一阶段。
    k2 = f(t[n]+h/2, y[n]+h*k1/2)  # RK4 第二阶段。
    k3 = f(t[n]+h/2, y[n]+h*k2/2)  # RK4 第三阶段。
    k4 = f(t[n]+h, y[n]+h*k3)  # RK4 第四阶段。
    y[n+1] = y[n] + h*(k1+2*k2+2*k3+k4)/6  # 保存起步值。
if method == "Nystrom":  # Nyström/Leapfrog：y_{n+1}=y_{n-1}+2h*f_n。
    for n in range(1, N):  # 从已有 y_0,y_1 开始推进。
        y[n+1] = y[n-1] + 2*h*f(t[n], y[n])  # 二步显式中心公式；寄生根 -1 可能引起振荡。
elif method in {"Milne", "Hamming"}:  # Python 集合语法：判断方法是否属于两个可选名称。
    for n in range(3, N):  # 已有 y_0,...,y_3 后开始推进。
        fp = f(t[n], y[n])  # 当前函数值 f_n。
        fm1 = f(t[n-1], y[n-1])  # 历史函数值 f_{n-1}。
        fm2 = f(t[n-2], y[n-2])  # 历史函数值 f_{n-2}。
        y_pred = y[n-3] + 4*h*(2*fp-fm1+2*fm2)/3  # Milne 四步显式预测公式。
        f_pred = f(t[n+1], y_pred)  # 在预测点计算 f_{n+1}。
        if method == "Milne":  # Milne 校正公式。
            y[n+1] = y[n-1] + h*(f_pred+4*fp+fm1)/3  # Simpson 型隐式校正，使用预测函数值一次。
        else:  # Hamming 校正分支。
            y[n+1] = (9*y[n]-y[n-2]+3*h*(f_pred+2*fp-fm1))/8  # Hamming 四阶校正公式。
else:  # 捕获未知方法名称。
    raise ValueError(f"未知 method：{method}")  # 主动抛出异常。
for tn, yn in zip(t, y):  # 配对输出结果。
    print(f"t={tn:.3f}, {method}={yn:.10f}")  # 格式化打印。
plt.plot(t, y, "o-", linewidth=1.3, label=method)  # 画数值解。
plt.xlabel("t")  # 设置横轴名称。
plt.ylabel("y")  # 设置纵轴名称。
plt.grid(True)  # 打开网格。
plt.legend()  # 显示图例。
plt.tight_layout()  # 自动调整布局。
plt.show()  # 显示图像。
