# 方程组与高阶 ODE：先化为一阶系统，再统一使用经典 RK4。  # 数值分析：m 阶标量方程可转换为 m 维一阶系统。
import math  # 标准库 math 用于默认方程中的指数与正弦函数。
import numpy as np  # NumPy 用于状态向量、网格和二维结果数组。
import matplotlib.pyplot as plt  # Matplotlib 用于一次画出所有分量。

F = lambda t, w: np.array([w[1], math.exp(2*t)*math.sin(t)-2*w[0]+2*w[1]], dtype=float)  # 【必须替换】二阶 y''=g(t,y,y') 写成 [y',g]；一般方程组直接返回各分量导数。
t0, T, w0, h = 0.0, 1.0, np.array([-0.4,-0.6], dtype=float), 0.1  # 【必须替换】w0=[y(t0),y'(t0),...]、区间和步长。
N = round((T-t0)/h)  # 计算固定步长总步数。
assert abs(t0+N*h-T) < 1e-12, "步长必须整除区间长度。"  # 检查最终节点。
t = t0 + np.arange(N+1)*h  # 创建等距时间网格。
m = w0.size  # 系统一阶方程个数等于状态向量维数。
w = np.zeros((N+1, m), dtype=float)  # 每一行保存一个时间点的全部状态分量。
w[0] = w0  # 写入初始状态。
for n in range(N):  # 对向量系统逐步推进。
    K1 = F(t[n], w[n])  # 第一阶段斜率向量。
    K2 = F(t[n]+h/2, w[n]+h*K1/2)  # 第二阶段在预测中点计算。
    K3 = F(t[n]+h/2, w[n]+h*K2/2)  # 第三阶段再次修正中点。
    K4 = F(t[n]+h, w[n]+h*K3)  # 第四阶段在右端点计算。
    w[n+1] = w[n] + h*(K1+2*K2+2*K3+K4)/6  # 向量 RK4 加权公式与标量公式完全相同。
for index, tn in enumerate(t):  # enumerate 同时给出数组下标和时间节点。
    print(f"t={tn:.3f}, w={w[index]}")  # NumPy 数组会直接打印所有状态分量。
for component in range(m):  # 遍历并绘制每一个状态分量。
    plt.plot(t, w[:,component], linewidth=1.3, label=f"w{component+1}")  # NumPy 切片 w[:,component] 取指定分量的全部时间值。
plt.xlabel("t")  # 设置横轴名称。
plt.ylabel("components")  # 设置纵轴名称。
plt.grid(True)  # 打开网格。
plt.legend()  # 显示图例。
plt.tight_layout()  # 自动调整布局。
plt.show()  # 显示图像。
