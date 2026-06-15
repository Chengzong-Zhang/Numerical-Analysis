# BDF/Gear 1 至 6 阶隐式模板。  # 数值分析：BDF 特别适合刚性问题，但每一步需解隐式方程。
import math  # 标准库 math 用于默认例题中的三角函数。
import numpy as np  # NumPy 用于系数、网格和历史值数组。
import matplotlib.pyplot as plt  # Matplotlib 用于画数值解与精确解。

f = lambda t, y: -10*(y-math.cos(t))-math.sin(t)  # 【必须替换】ODE 右端；默认例题精确解为 cos(t)。
fy = lambda t, y: -10.0  # 【必须替换】偏导 df/dy；Newton 解隐式方程时必须提供。
t0, T, y0, h, k = 0.0, 1.0, 1.0, 0.05, 2  # 【必须替换】区间、初值、步长、BDF 阶数 1...6。
tol, max_iter = 1e-12, 200  # Newton 迭代停止阈值与最大迭代次数。
alpha = [np.array([1,-1]),np.array([1,-4/3,1/3]),np.array([1,-18/11,9/11,-2/11]),np.array([1,-48/25,36/25,-16/25,3/25]),np.array([1,-300/137,300/137,-200/137,75/137,-12/137]),np.array([1,-360/147,450/147,-400/147,225/147,-72/147,10/147])]  # 归一化 BDF1...BDF6 左端系数。
beta = np.array([1,2/3,6/11,12/25,60/137,60/147], dtype=float)  # 右端 h*beta_k*f(t_{n+1},y_{n+1}) 系数。
assert k in range(1, 7), "k 必须是 1 到 6。"  # 检查 BDF 阶数。
a, b = alpha[k-1], beta[k-1]  # 取所选 BDF 的左端系数和右端系数。
N = round((T-t0)/h)  # 计算固定步长总步数。
assert N >= k-1 and abs(t0+N*h-T) < 1e-12, "区间或步长不满足 BDF 要求。"  # 检查起步点数量与终点。
t = t0 + np.arange(N+1)*h  # 创建等距时间网格。
y = np.zeros(N+1, dtype=float)  # 预分配数值解。
y[0] = y0  # 写入初值。
for n in range(k-1):  # 用 RK4 提供 k-1 个起步值。
    q1 = f(t[n], y[n])  # RK4 第一阶段。
    q2 = f(t[n]+h/2, y[n]+h*q1/2)  # RK4 第二阶段。
    q3 = f(t[n]+h/2, y[n]+h*q2/2)  # RK4 第三阶段。
    q4 = f(t[n]+h, y[n]+h*q3)  # RK4 第四阶段。
    y[n+1] = y[n] + h*(q1+2*q2+2*q3+q4)/6  # 保存起步值。
for n in range(k-1, N):  # 从已有 k 个历史值后开始 BDF 推进。
    known = -sum(a[j+1]*y[n-j] for j in range(k))  # Python 生成器求和；把历史项移到方程右端。
    z = y[n]  # 用最近数值解作为 y_{n+1} 的 Newton 初猜。
    for iteration in range(max_iter):  # 解 G(z)=z-known-h*b*f(t_{n+1},z)=0。
        G = z-known-h*b*f(t[n+1], z)  # 隐式 BDF 方程残差。
        dG = 1-h*b*fy(t[n+1], z)  # 残差导数；向量题应换成 I-h*b*Jacobian。
        z_new = z-G/dG  # 标量 Newton 更新；向量题使用 np.linalg.solve(dG,G)。
        if abs(z_new-z) < tol: break  # Newton 修正足够小时停止。
        z = z_new  # 更新猜测并继续 Newton。
    assert iteration < max_iter-1, "BDF Newton 未收敛，请检查 fy、减小 h 或调整初猜。"  # 显式报告隐式求解失败。
    y[n+1] = z_new  # 保存本步 BDF 解。
for tn, yn in zip(t, y):  # 配对输出结果。
    print(f"t={tn:.3f}, BDF{k}_y={yn:.10f}")  # 格式化打印。
plt.plot(t, y, "o-", linewidth=1.3, label=f"BDF{k}")  # 画 BDF 数值解。
plt.plot(t, np.cos(t), "k--", linewidth=1.3, label="exact")  # NumPy cos 向量化计算默认精确解。
plt.xlabel("t")  # 设置横轴名称。
plt.ylabel("y")  # 设置纵轴名称。
plt.grid(True)  # 打开网格。
plt.legend()  # 显示图例。
plt.tight_layout()  # 自动调整布局。
plt.show()  # 显示图像。
