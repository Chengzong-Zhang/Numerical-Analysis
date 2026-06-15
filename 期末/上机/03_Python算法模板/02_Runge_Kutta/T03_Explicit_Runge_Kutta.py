# 一般显式 Runge-Kutta：用 Butcher 表覆盖三阶 Heun、三阶 RK 和经典 RK4。  # 数值分析：显式 RK 的 A 矩阵严格下三角。
import numpy as np  # NumPy 用于向量状态、矩阵乘法和时间网格。
import matplotlib.pyplot as plt  # Matplotlib 用于画第一分量。

f = lambda t, y: -y + t**2 + 1  # 【必须替换】标量或向量 ODE 右端；向量题应返回 NumPy 一维数组。
t0, T, y0, h = 0.0, 1.0, 1.0, 0.1  # 【必须替换】区间、初值、步长；向量初值写 np.array([a,b])。
method = "RK4"  # 【必须替换】可选 Heun3/RK3/RK4。
if method == "Heun3":  # 教材三级三阶 Heun 方法。
    A = np.array([[0,0,0],[1/3,0,0],[0,2/3,0]], dtype=float)  # A 的第 r 行表示第 r 阶段对前面阶段的组合。
    b = np.array([1/4,0,3/4], dtype=float)  # b 是最后一步组合各阶段的权重。
    c = np.array([0,1/3,2/3], dtype=float)  # c 是各阶段所在的相对时间位置。
elif method == "RK3":  # 教材公式 (5.3.13) 的三级三阶 RK。
    A = np.array([[0,0,0],[1/2,0,0],[0,3/4,0]], dtype=float)  # 严格下三角结构保证阶段可依次算出。
    b = np.array([2/9,3/9,4/9], dtype=float)  # 最终权重和为 1，保证相容性。
    c = np.array([0,1/2,3/4], dtype=float)  # 阶段节点。
elif method == "RK4":  # 经典四级四阶 Runge-Kutta。
    A = np.array([[0,0,0,0],[1/2,0,0,0],[0,1/2,0,0],[0,0,1,0]], dtype=float)  # 经典 RK4 阶段矩阵。
    b = np.array([1/6,1/3,1/3,1/6], dtype=float)  # 经典 RK4 最终权重。
    c = np.array([0,1/2,1/2,1], dtype=float)  # 经典 RK4 阶段时间位置。
else:  # 捕获未知方法名。
    raise ValueError(f"未知 method：{method}")  # 抛出异常并显示错误方法名。
N = round((T-t0)/h)  # 计算总步数。
assert abs(t0+N*h-T) < 1e-12, "步长必须整除区间长度。"  # 检查固定网格终点。
t = t0 + np.arange(N+1)*h  # 创建等距时间网格。
y0 = np.atleast_1d(np.asarray(y0, dtype=float))  # NumPy 语法：把标量或列表统一成一维浮点数组。
y = np.zeros((N+1, y0.size), dtype=float)  # 每一行保存一个时间点的状态向量。
y[0] = y0  # 写入初始状态。
s = b.size  # RK 级数等于权重向量长度。
for n in range(N):  # 逐步推进。
    K = np.zeros((s, y0.size), dtype=float)  # K[r] 保存第 r 个阶段的斜率向量。
    for r in range(s):  # 依次计算显式阶段。
        stage_y = y[n] + h*(A[r] @ K)  # Python/NumPy 语法：@ 表示矩阵乘法；组合已算阶段。
        K[r] = np.atleast_1d(f(t[n]+c[r]*h, stage_y))  # 在阶段点计算斜率，并统一为一维数组。
    y[n+1] = y[n] + h*(b @ K)  # 按权重 b 组合全部阶段，完成一步。
for tn, yn in zip(t, y[:,0]):  # 只打印第一分量；方程组可自行增加其他分量。
    print(f"t={tn:.3f}, y1={yn:.10f}")  # 格式化输出第一分量。
plt.plot(t, y[:,0], "o-", linewidth=1.3, label=method)  # 画第一分量数值解。
plt.xlabel("t")  # 设置横轴名称。
plt.ylabel("y[0]")  # 设置纵轴名称。
plt.grid(True)  # 打开网格。
plt.legend()  # 显示图例。
plt.tight_layout()  # 自动调整布局。
plt.show()  # 显示图像。
