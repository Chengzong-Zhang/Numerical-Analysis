# Euler 与教材二阶单步法统一模板。  # 覆盖向前/向后 Euler、隐式梯形、改进 Euler、中点、二阶 Heun。
import numpy as np  # NumPy 用于创建固定网格和预分配数值解数组。
import matplotlib.pyplot as plt  # Matplotlib 用于画数值解曲线。

f = lambda t, y: -y + t**2 + 1  # 【必须替换】ODE 右端 f(t,y)；Python 的 ** 表示乘方。
t0, T, y0, h = 0.0, 1.0, 1.0, 0.1  # 【必须替换】区间、初值与固定步长。
method = "improved_euler"  # 【必须替换】可选 forward_euler/backward_euler/trapezoid/improved_euler/midpoint/heun2。
tol, max_iter = 1e-12, 100  # 隐式法不动点迭代的停止阈值与最大迭代次数。
N = round((T-t0)/h)  # 固定步长方法总步数。
assert abs(t0+N*h-T) < 1e-12, "步长必须整除区间长度。"  # 检查最后节点是否恰好等于 T。
t = t0 + np.arange(N+1)*h  # 创建等距节点 t_0,...,t_N。
y = np.zeros(N+1, dtype=float)  # 预分配数值解，避免循环中反复扩容。
y[0] = y0  # 写入初值 y(t0)=y0。
for n in range(N):  # 逐步推进到下一节点。
    tn, yn = t[n], y[n]  # 保存当前时间和当前近似值，便于与公式对应。
    if method == "forward_euler":  # 向前 Euler：显式、一阶，实轴稳定区间为 (-2,0)。
        y[n+1] = yn + h*f(tn, yn)  # 数值公式 y_{n+1}=y_n+h*f_n。
    elif method == "backward_euler":  # 向后 Euler：隐式、一阶、A 稳定。
        z = yn + h*f(tn, yn)  # 用向前 Euler 预测值作为隐式未知量的初猜。
        for _ in range(max_iter):  # 下划线 _ 表示循环计数本身不需要使用。
            z_new = yn + h*f(t[n+1], z)  # 不动点更新 z^{k+1}=y_n+h*f(t_{n+1},z^k)。
            if abs(z_new-z) < tol: break  # Python 单行 if；相邻迭代差足够小时停止。
            z = z_new  # 接受新猜测并继续迭代。
        y[n+1] = z_new  # 保存向后 Euler 隐式解。
    elif method == "trapezoid":  # 隐式梯形：二阶、A 稳定，来自积分形式的梯形求积。
        z = yn + h*f(tn, yn)  # 用 Euler 预测值启动迭代。
        for _ in range(max_iter):  # 对隐式梯形方程做不动点迭代。
            z_new = yn + h*(f(tn, yn)+f(t[n+1], z))/2  # 梯形公式 z=y_n+h/2*(f_n+f_{n+1})。
            if abs(z_new-z) < tol: break  # 达到迭代容差后停止。
            z = z_new  # 更新隐式未知量。
        y[n+1] = z_new  # 保存梯形法结果。
    elif method == "improved_euler":  # 改进 Euler，又称显式梯形或 Euler 预测校正，整体二阶。
        k1 = f(tn, yn)  # 第一阶段：当前点斜率。
        k2 = f(tn+h, yn+h*k1)  # 第二阶段：在 Euler 预测点计算右端斜率。
        y[n+1] = yn + h*(k1+k2)/2  # 对两端斜率取平均，完成校正。
    elif method == "midpoint":  # 变形 Euler/中点 RK2，整体二阶。
        k1 = f(tn, yn)  # 当前点斜率。
        k2 = f(tn+h/2, yn+h*k1/2)  # 使用 k1 预测中点，再求中点斜率。
        y[n+1] = yn + h*k2  # 使用中点斜率推进整步。
    elif method == "heun2":  # 教材二阶 Heun：a2=2/3，最终权重为 1/4 和 3/4。
        k1 = f(tn, yn)  # 第一阶段斜率。
        k2 = f(tn+2*h/3, yn+2*h*k1/3)  # 在 2/3 步位置计算第二阶段斜率。
        y[n+1] = yn + h*(k1+3*k2)/4  # 加权组合满足二阶 Runge-Kutta 阶条件。
    else:  # 捕获拼写错误或未支持的方法名称。
        raise ValueError(f"未知 method：{method}")  # Python 语法：raise 主动抛出异常并终止。
for tn, yn in zip(t, y):  # 配对输出每个节点的数值结果。
    print(f"t={tn:.3f}, {method}={yn:.10f}")  # 使用 f-string 打印格式化结果。
plt.plot(t, y, "o-", linewidth=1.3, label=method)  # 画所选方法的数值解。
plt.xlabel("t")  # 设置横轴名称。
plt.ylabel("y")  # 设置纵轴名称。
plt.grid(True)  # 打开网格。
plt.legend()  # 显示图例。
plt.tight_layout()  # 自动调整布局。
plt.show()  # 显示图像。
