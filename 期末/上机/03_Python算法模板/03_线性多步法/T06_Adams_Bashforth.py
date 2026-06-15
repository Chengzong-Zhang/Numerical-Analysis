# Adams-Bashforth 显式外插公式：统一覆盖 1 至 6 步。  # 数值分析：s 步 AB 通常为 s 阶，需要 s 个起步值。
import math  # 标准库 math 用于默认例题中的指数函数。
import numpy as np  # NumPy 用于系数数组、网格和点积。
import matplotlib.pyplot as plt  # Matplotlib 用于绘图。

f = lambda t, y: -y + math.exp(-t)  # 【必须替换】ODE 右端。
t0, T, y0, h, steps = 0.0, 1.0, 0.0, 0.1, 4  # 【必须替换】区间、初值、步长、AB 步数 1...6。
AB = [np.array([1]),np.array([3,-1])/2,np.array([23,-16,5])/12,np.array([55,-59,37,-9])/24,np.array([1901,-2774,2616,-1274,251])/720,np.array([4277,-7923,9982,-7298,2877,-475])/1440]  # 系数按 [f_n,f_{n-1},...] 排列。
assert steps in range(1, 7), "steps 必须是 1 到 6。"  # 检查所选步数是否合法。
coef = AB[steps-1]  # Python 下标从 0 开始，因此第 steps 个系数位于 steps-1。
N = round((T-t0)/h)  # 计算固定步长总步数。
assert N >= steps-1 and abs(t0+N*h-T) < 1e-12, "区间或步长不满足多步法要求。"  # 检查终点与起步点数量。
t = t0 + np.arange(N+1)*h  # 创建等距时间网格。
y = np.zeros(N+1, dtype=float)  # 预分配标量数值解。
y[0] = y0  # 写入初值。
for n in range(steps-1):  # 使用 RK4 计算 y_1,...,y_{steps-1}。
    k1 = f(t[n], y[n])  # RK4 第一阶段斜率。
    k2 = f(t[n]+h/2, y[n]+h*k1/2)  # RK4 第二阶段斜率。
    k3 = f(t[n]+h/2, y[n]+h*k2/2)  # RK4 第三阶段斜率。
    k4 = f(t[n]+h, y[n]+h*k3)  # RK4 第四阶段斜率。
    y[n+1] = y[n] + h*(k1+2*k2+2*k3+k4)/6  # 高阶起步，避免启动误差降低多步法整体阶数。
for n in range(steps-1, N):  # 从已经拥有 steps 个历史点的位置开始 AB 推进。
    f_history = np.array([f(t[n-j], y[n-j]) for j in range(steps)], dtype=float)  # 列表推导式收集 [f_n,f_{n-1},...]。
    y[n+1] = y[n] + h*(coef @ f_history)  # NumPy @ 点积实现 Adams-Bashforth 外插公式。
for tn, yn in zip(t, y):  # 配对输出节点与数值解。
    print(f"t={tn:.3f}, AB{steps}_y={yn:.10f}")  # 格式化输出所选 AB 方法结果。
plt.plot(t, y, "o-", linewidth=1.3, label=f"AB{steps}")  # 画数值解。
plt.xlabel("t")  # 设置横轴名称。
plt.ylabel("y")  # 设置纵轴名称。
plt.grid(True)  # 打开网格。
plt.legend()  # 显示图例。
plt.tight_layout()  # 自动调整布局。
plt.show()  # 显示图像。
