# 误差、收敛阶、Richardson 外推与绝对稳定性模板。  # 数值分析：用步长缩小后的误差比估计整体收敛阶。
import math  # 标准库 math 用于指数、对数运算。
import numpy as np  # NumPy 用于步长、误差和稳定函数数组。
import matplotlib.pyplot as plt  # Matplotlib 用于画稳定函数。

exact = lambda t: math.exp(-t)  # 【必须替换】精确解；未知时可用更小步长高阶解作参考。
approx = lambda h: (1-h)**round(1/h)  # 【必须替换】返回终点数值解；默认是 Euler 解 y'=-y 到 t=1。
hs = np.array([0.2,0.1,0.05,0.025], dtype=float)  # 【必须替换】逐次减半最方便估计阶数。
errors = np.zeros_like(hs)  # NumPy 语法：zeros_like 创建与 hs 同形状的零数组。
for i, hi in enumerate(hs):  # enumerate 同时给出下标和当前步长。
    errors[i] = abs(approx(hi)-exact(1.0))  # 终点绝对误差；若题目要求最大误差，应对所有节点取 max。
orders = np.full_like(hs, np.nan)  # 第一组误差没有前一组可比较，因此先填充 NaN。
for i in range(1, hs.size):  # 从第二个步长开始估计收敛阶。
    orders[i] = math.log(errors[i-1]/errors[i])/math.log(hs[i-1]/hs[i])  # 一般步长比的阶数估计公式。
assumed_order = 1  # 【必须替换】Richardson 外推依据的方法理论阶数。
y_coarse, y_fine = approx(hs[-2]), approx(hs[-1])  # Python 负下标 -2/-1 取最后两组粗细网格结果。
ratio = hs[-2]/hs[-1]  # 计算粗细步长之比。
y_richardson = y_fine+(y_fine-y_coarse)/(ratio**assumed_order-1)  # 消去主误差项，提高精度。
for hi, error, order in zip(hs, errors, orders):  # 配对输出误差和估计阶数。
    print(f"h={hi:.6g}, error={error:.6e}, estimated_order={order}")  # 格式化打印；第一项 order 为 nan。
print(f"Richardson 外推终值 = {y_richardson:.15g}")  # 打印外推结果。
z = np.linspace(-5, 1, 1200)  # 【必须替换】z=lambda*h；沿实轴检查稳定区间。
R_euler = 1+z  # 向前 Euler 稳定函数 R(z)=1+z。
R_rk2 = 1+z+z**2/2  # 改进 Euler/RK2 稳定函数。
R_rk4 = 1+z+z**2/2+z**3/6+z**4/24  # 经典 RK4 稳定函数。
plt.plot(z, np.abs(R_euler), label="Euler")  # 画 Euler 的 |R(z)|。
plt.plot(z, np.abs(R_rk2), label="RK2")  # 画 RK2 的 |R(z)|。
plt.plot(z, np.abs(R_rk4), label="RK4")  # 画 RK4 的 |R(z)|。
plt.axhline(1, color="black", linestyle="--", label="|R|=1")  # |R(z)|<1 的区域对应绝对稳定。
plt.ylim(0, 3)  # 限制纵轴范围，突出稳定边界。
plt.xlabel("z=lambda*h")  # 设置横轴名称。
plt.ylabel("|R(z)|")  # 设置纵轴名称。
plt.grid(True)  # 打开网格。
plt.legend()  # 显示图例。
plt.tight_layout()  # 自动调整布局。
plt.show()  # 显示图像。
