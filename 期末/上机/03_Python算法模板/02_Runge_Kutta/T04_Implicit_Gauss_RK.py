# 隐式 Gauss-Runge-Kutta：覆盖一级二阶、二级四阶和三级六阶方法。  # 数值分析：每一步必须联立求解所有阶段。
import math  # 标准库 math 用于平方根。
import numpy as np  # NumPy 用于 Butcher 系数、阶段向量和时间网格。
import matplotlib.pyplot as plt  # Matplotlib 用于绘图。

f = lambda t, y: -10*(y-math.cos(t))-math.sin(t)  # 【必须替换】ODE 右端；默认例题带快速衰减项，精确解为 cos(t)。
t0, T, y0, h = 0.0, 1.0, 1.0, 0.1  # 【必须替换】区间、初值、步长。
method, tol, max_iter = "Gauss4", 1e-12, 100  # 【必须替换】Gauss2/Gauss4/Gauss6；隐式迭代容差与上限。
if method == "Gauss2":  # 一级二阶隐式中点法。
    A = np.array([[1/2]], dtype=float)  # 单阶段仍包含自身斜率，因此是隐式方法。
    b = np.array([1.0], dtype=float)  # 最终阶段权重。
    c = np.array([1/2], dtype=float)  # Gauss 单节点位于区间中点。
elif method == "Gauss4":  # 二级四阶 Gauss-RK。
    q = math.sqrt(3)  # 两节点 Gauss 公式中的根号系数。
    A = np.array([[1/4,1/4-q/6],[1/4+q/6,1/4]], dtype=float)  # A 含对角和上三角项，阶段彼此耦合。
    b = np.array([1/2,1/2], dtype=float)  # 两个阶段等权组合。
    c = np.array([1/2-q/6,1/2+q/6], dtype=float)  # 两个 Gauss 节点。
elif method == "Gauss6":  # 三级六阶 Gauss-RK。
    q = math.sqrt(15)  # 三节点 Gauss 公式中的根号系数。
    A = np.array([[5/36,2/9-q/15,5/36-q/30],[5/36+q/24,2/9,5/36-q/24],[5/36+q/30,2/9+q/15,5/36]], dtype=float)  # 三级耦合阶段矩阵。
    b = np.array([5/18,4/9,5/18], dtype=float)  # 六阶方法的最终权重。
    c = np.array([1/2-q/10,1/2,1/2+q/10], dtype=float)  # 三个 Gauss 节点。
else:  # 捕获未知方法名。
    raise ValueError(f"未知 method：{method}")  # 主动抛出异常。
N = round((T-t0)/h)  # 计算固定步长总步数。
assert abs(t0+N*h-T) < 1e-12, "步长必须整除区间长度。"  # 检查网格终点。
t = t0 + np.arange(N+1)*h  # 创建等距时间网格。
y = np.zeros(N+1, dtype=float)  # 预分配标量数值解。
y[0] = y0  # 写入初值。
s = b.size  # 阶段数等于权重个数。
for n in range(N):  # 逐步推进。
    K = np.full(s, f(t[n], y[n]), dtype=float)  # 用当前斜率填充全部阶段，作为不动点初猜。
    for iteration in range(max_iter):  # 对耦合阶段方程做不动点迭代。
        K_new = np.zeros(s, dtype=float)  # 预分配本轮新阶段。
        for r in range(s):  # 逐个更新阶段。
            stage_y = y[n] + h*(A[r] @ K)  # 当前全部阶段共同决定第 r 个阶段状态。
            K_new[r] = f(t[n]+c[r]*h, stage_y)  # 计算新的第 r 阶段斜率。
        if np.linalg.norm(K_new-K, ord=np.inf) < tol: break  # 无穷范数控制所有阶段的最大变化。
        K = K_new  # 接受本轮阶段并继续迭代。
    assert iteration < max_iter-1, "隐式阶段未收敛：请减小 h 或改用 Newton。"  # 防止静默使用错误阶段值。
    y[n+1] = y[n] + h*(b @ K_new)  # 用最终阶段加权组合完成一步。
for tn, yn in zip(t, y):  # 配对输出各节点结果。
    print(f"t={tn:.3f}, {method}={yn:.10f}")  # 格式化输出。
plt.plot(t, y, "o-", linewidth=1.3, label=method)  # 画隐式 Gauss-RK 数值解。
plt.xlabel("t")  # 设置横轴名称。
plt.ylabel("y")  # 设置纵轴名称。
plt.grid(True)  # 打开网格。
plt.legend()  # 显示图例。
plt.tight_layout()  # 自动调整布局。
plt.show()  # 显示图像。
