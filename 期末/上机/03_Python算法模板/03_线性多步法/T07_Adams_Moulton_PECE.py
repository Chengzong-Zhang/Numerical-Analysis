# Adams-Moulton 隐式内插与 Adams 预测-校正 PECE。  # 数值分析：AB 负责预测，AM 负责校正。
import math  # 标准库 math 用于默认例题的指数函数。
import numpy as np  # NumPy 用于系数数组、历史值数组与点积。
import matplotlib.pyplot as plt  # Matplotlib 用于绘图。

f = lambda t, y: -y + math.exp(-t)  # 【必须替换】ODE 右端。
t0, T, y0, h, order = 0.0, 1.0, 0.0, 0.1, 4  # 【必须替换】区间、初值、步长、AM 阶数 1...6。
mode, tol, max_iter = "PECE", 1e-12, 100  # 【必须替换】可选 PECE/implicit；隐式迭代容差与上限。
AB = [np.array([1]),np.array([3,-1])/2,np.array([23,-16,5])/12,np.array([55,-59,37,-9])/24,np.array([1901,-2774,2616,-1274,251])/720,np.array([4277,-7923,9982,-7298,2877,-475])/1440]  # AB 预测权重。
AM = [np.array([1]),np.array([1,1])/2,np.array([5,8,-1])/12,np.array([9,19,-5,1])/24,np.array([251,646,-264,106,-19])/720,np.array([475,1427,-798,482,-173,27])/1440]  # AM 权重按 [f_{n+1},f_n,...] 排列。
assert order in range(1, 7), "order 必须是 1 到 6。"  # 检查阶数是否合法。
ab, am = AB[order-1], AM[order-1]  # 取相同阶数的预测与校正系数。
N = round((T-t0)/h)  # 计算固定步长总步数。
assert N >= order-1 and abs(t0+N*h-T) < 1e-12, "区间或步长不满足多步法要求。"  # 检查起步点数量和终点。
t = t0 + np.arange(N+1)*h  # 创建固定网格。
y = np.zeros(N+1, dtype=float)  # 预分配数值解。
y[0] = y0  # 写入初值。
for n in range(order-1):  # 用 RK4 提供多步法所需起步值。
    k1 = f(t[n], y[n])  # RK4 第一阶段。
    k2 = f(t[n]+h/2, y[n]+h*k1/2)  # RK4 第二阶段。
    k3 = f(t[n]+h/2, y[n]+h*k2/2)  # RK4 第三阶段。
    k4 = f(t[n]+h, y[n]+h*k3)  # RK4 第四阶段。
    y[n+1] = y[n] + h*(k1+2*k2+2*k3+k4)/6  # 保存高精度起步值。
for n in range(order-1, N):  # 主多步循环。
    f_history = np.array([f(t[n-j], y[n-j]) for j in range(order)], dtype=float)  # 收集 [f_n,f_{n-1},...]。
    z_pred = y[n] + h*(ab @ f_history)  # P：用 Adams-Bashforth 显式预测 y_{n+1}。
    old_history = f_history[:max(order-1, 0)]  # AM 除 f_{n+1} 外使用 [f_n,...]；切片不会越界。
    if mode == "PECE":  # PECE 只使用预测值校正一次，计算量较小。
        z = y[n] + h*(am[0]*f(t[n+1], z_pred) + am[1:] @ old_history)  # E+C：把预测函数值代入 AM。
    elif mode == "implicit":  # 完全隐式 AM 反复求解封闭公式。
        z = z_pred  # 用预测值作为隐式迭代初猜。
        for _ in range(max_iter):  # 对 AM 隐式方程做不动点迭代。
            z_new = y[n] + h*(am[0]*f(t[n+1], z) + am[1:] @ old_history)  # 更新隐式未知量。
            if abs(z_new-z) < tol: break  # 相邻迭代差小于容差时停止。
            z = z_new  # 更新猜测并继续迭代。
        z = z_new  # 保存最终隐式迭代值。
    else:  # 捕获未知模式。
        raise ValueError(f"未知 mode：{mode}")  # 主动抛出异常。
    y[n+1] = z  # E：保存校正值，供后续步骤使用。
for tn, yn in zip(t, y):  # 配对输出结果。
    print(f"t={tn:.3f}, AM{order}_{mode}={yn:.10f}")  # 格式化打印。
plt.plot(t, y, "o-", linewidth=1.3, label=f"AM{order} {mode}")  # 画数值解。
plt.xlabel("t")  # 设置横轴名称。
plt.ylabel("y")  # 设置纵轴名称。
plt.grid(True)  # 打开网格。
plt.legend()  # 显示图例。
plt.tight_layout()  # 自动调整布局。
plt.show()  # 显示图像。
