# Taylor 级数法：y_{n+1}=y_n+sum_{j=1}^p h^j/j!*y^(j)(t_n)。  # 数值分析：截断到 p 阶导数时通常得到 p 阶方法。
import math  # Python 语法：导入标准库 math，用 factorial 计算阶乘。
import numpy as np  # Python 语法：给 NumPy 起别名 np，用于数组与等距网格。
import matplotlib.pyplot as plt  # Python 语法：导入绘图库 pyplot，并使用别名 plt。

derivatives = [lambda t, y: -y + t**2 + 1, lambda t, y: y - t**2 + 2*t - 1]  # 【必须替换】依次填写 y'、y''、...；lambda 创建匿名函数。
t0, T, y0, h = 0.0, 1.0, 1.0, 0.1  # 【必须替换】左端点、右端点、初值、固定步长；逗号可一次赋多个变量。
p = len(derivatives)  # Python 语法：len 返回列表长度；数值分析：导数函数个数就是 Taylor 方法阶数 p。
N = round((T-t0)/h)  # Python 语法：round 取最近整数；固定步长方法需要整数步数。
assert abs(t0+N*h-T) < 1e-12, "步长必须整除区间长度。"  # Python 语法：assert 条件失败时停止，防止实际终点错误。
t = t0 + np.arange(N+1)*h  # NumPy 语法：arange(N+1) 生成 0,...,N，再向量化得到全部节点。
y = np.zeros(N+1, dtype=float)  # NumPy 语法：zeros 预分配浮点数组；N 步需要 N+1 个节点值。
y[0] = y0  # Python 下标从 0 开始，因此 y[0] 对应数学记号 y_0。
for n in range(N):  # Python 语法：range(N) 生成 0,...,N-1；数值分析：逐步从 t_n 推进到 t_{n+1}。
    increment = 0.0  # 本步 Taylor 增量先置零。
    for j, derivative in enumerate(derivatives, start=1):  # enumerate 同时给出阶数 j 和第 j 阶导数函数。
        increment += h**j/math.factorial(j)*derivative(t[n], y[n])  # 数值分析：累加 h^j/j!*y^(j)(t_n)。
    y[n+1] = y[n] + increment  # Taylor 截断公式完成一步推进，局部截断误差通常为 O(h^(p+1))。
for tn, yn in zip(t, y):  # Python 语法：zip 把时间与数值解逐项配对，方便输出表格。
    print(f"t={tn:.3f}, Taylor_y={yn:.10f}")  # Python 语法：f-string 控制小数位并插入变量。
plt.plot(t, y, "o-", linewidth=1.3, label=f"Taylor order {p}")  # 绘图语法：o- 表示节点圆点加连线。
plt.xlabel("t")  # 绘图语法：设置横轴名称。
plt.ylabel("y")  # 绘图语法：设置纵轴名称。
plt.grid(True)  # 绘图语法：打开网格，便于观察数值变化。
plt.legend()  # 绘图语法：显示图例。
plt.tight_layout()  # 绘图语法：自动调整边距，防止文字被裁切。
plt.show()  # 绘图语法：显示图像；批量验证时使用 Agg 后端，不会弹窗。
