"""教材第五章 5.3.4：变形 Euler（二阶中点法）与三阶 RK 组成的自适应嵌套 RK2(3) 模板。"""  # Python 语法：三引号字符串放在文件开头作为模块说明；数值分析：嵌套方法共享阶段斜率并用二、三阶结果之差估计局部误差。

import numpy as np  # Python 语法：import 导入第三方库并取别名 np；这里用 NumPy 计算绝对值、保存浮点数组并支持以后改成方程组。
import matplotlib.pyplot as plt  # Python 语法：从 Matplotlib 导入 pyplot 并取别名 plt；这里用于观察自适应算法实际接受的步长。


def f(t, y):  # 【必须替换】Python 语法：def 定义函数，t 是自变量，y 是当前数值解；数值分析：返回初值问题 y'=f(t,y) 的右端值。
    return -y + t**2 + 1.0  # 【必须替换】Python 语法：return 返回表达式，** 表示乘方；请把本行替换成题目给出的 f(t,y)。


t0 = 0.0  # 【必须替换】Python 语法：= 给变量赋值；数值分析：t0 是求解区间左端点，也是初值所在点。
T = 1.0  # 【必须替换】数值分析：T 是求解区间右端点，算法将通过截短最后一步恰好到达 T。
y0 = 1.0  # 【必须替换】数值分析：y0 是初值 y(t0)；若求方程组，可改成 np.array([...], dtype=float)。
h_initial = 0.1  # 【必须替换】数值分析：初始步长只影响算法起步，之后步长会由误差估计自动调整。
h_min = 1.0e-6  # 【必须替换】数值分析：允许的最小步长；达到它仍不满足容差时停止，避免无限拒绝。
h_max = 0.25  # 【必须替换】数值分析：允许的最大步长；限制解变化平缓区间中步长增长过快。
atol = 1.0e-7  # 【必须替换】数值分析：绝对误差容限，主要控制数值解接近零时的误差。
rtol = 1.0e-5  # 【必须替换】数值分析：相对误差容限，使允许误差随数值解量级变化。
safety = 0.9  # 【通常不替换】数值分析：安全系数小于 1，使理论建议步长略微缩小，减少下一步被拒绝的概率。
factor_min = 0.2  # 【通常不替换】数值分析：限制一次调整时步长最多缩小到原来的 0.2 倍，避免剧烈振荡。
factor_max = 5.0  # 【通常不替换】数值分析：限制一次调整时步长最多增大到原来的 5 倍，避免误差偶然很小时激增。

t_values = [float(t0)]  # Python 语法：方括号创建列表，float 强制转成浮点数；自适应节点数未知，所以使用可动态增长的列表。
y_values = [np.asarray(y0, dtype=float)]  # Python 语法：np.asarray 把标量或列表统一转为数组；这样同一模板可处理标量和方程组。
accepted_h = []  # Python 语法：[] 创建空列表；这里保存每个被接受的实际步长，便于检查自适应效果。
error_history = []  # Python 语法：空列表稍后用 append 追加元素；这里保存每个被接受步骤的无量纲误差指标。
h = float(h_initial)  # Python 语法：创建当前候选步长变量；数值分析：每次尝试后都会依据后验误差更新它。
rejected_steps = 0  # Python 语法：整数变量用于计数；数值分析：拒绝次数过多通常说明容差过严、问题刚性或初始步长过大。

while t_values[-1] < T:  # Python 语法：while 在条件为真时循环，[-1] 取列表最后元素；数值分析：逐步推进直到右端点。
    tn = t_values[-1]  # Python 语法：读取列表末元素；数值分析：tn 是当前已经接受的时间节点。
    yn = y_values[-1]  # Python 语法：读取列表末元素；数值分析：yn 是当前节点处的数值近似。
    hn = min(h, h_max, T - tn)  # Python 语法：min 取最小值；数值分析：同时满足当前建议、最大步长和最后一步不越过终点。

    k1 = f(tn, yn)  # 数值分析：第一阶段斜率 k1=f(tn,yn)，被二阶与三阶公式共同使用。
    k2 = f(tn + hn / 2.0, yn + hn * k1 / 2.0)  # Python 语法：函数参数可写算术表达式；数值分析：在半步预测点计算第二阶段斜率。
    k3 = f(tn + 3.0 * hn / 4.0, yn + 3.0 * hn * k2 / 4.0)  # 数值分析：在四分之三步处计算第三阶段斜率，供三阶公式使用。
    y2 = yn + hn * k2  # 数值分析：变形 Euler/显式中点法的二阶近似，其权重为 (0,1,0)。
    y3 = yn + hn * (2.0 * k1 / 9.0 + k2 / 3.0 + 4.0 * k3 / 9.0)  # 数值分析：教材给出的三阶 RK 近似，与 y2 共享 k1、k2、k3。

    scale = atol + rtol * np.maximum(np.abs(yn), np.abs(y3))  # Python 语法：np.maximum 逐元素取较大值；数值分析：绝对与相对容差共同形成误差尺度。
    error_ratio = float(np.max(np.abs(y3 - y2) / scale))  # Python 语法：np.max 取各分量最大值并用 float 转标量；数值分析：嵌套解之差是二阶局部误差的后验估计。

    if error_ratio <= 1.0:  # Python 语法：if 仅在条件成立时执行缩进块；数值分析：无量纲误差不超过 1 表示本步满足用户容差。
        t_values.append(tn + hn)  # Python 语法：list.append 在末尾追加元素；数值分析：接受本步并保存新节点。
        y_values.append(y3)  # 数值分析：接受精度更高的三阶近似，而二阶近似只用于误差估计。
        accepted_h.append(hn)  # 数值分析：保存实际采用的步长，用于观察解变化快慢与步长之间的关系。
        error_history.append(error_ratio)  # 数值分析：保存已接受步的归一化误差，理论上每项都应不超过 1。
    else:  # Python 语法：else 处理 if 条件不成立的情况；数值分析：误差超限时不推进时间，并重新尝试当前步骤。
        rejected_steps += 1  # Python 语法：+= 1 表示计数器自增；数值分析：记录因误差超限而拒绝的候选步骤数。

    factor = factor_max if error_ratio == 0.0 else safety * error_ratio ** (-1.0 / 3.0)  # Python 语法：条件表达式 a if 条件 else b；数值分析：二阶方法局部误差约为 C*h^3，因此步长按误差的负三分之一次幂调整。
    factor = min(factor_max, max(factor_min, factor))  # Python 语法：嵌套 min/max 把调整倍数截断；数值分析：避免步长在相邻尝试间变化过猛。
    h = min(h_max, max(h_min, hn * factor))  # 数值分析：生成下一次候选步长，并强制限制在 [h_min,h_max] 内。

    if error_ratio > 1.0 and h <= h_min:  # Python 语法：and 要求两个条件同时成立；数值分析：误差超限且无法继续缩小步长时判定失败。
        raise RuntimeError("达到 h_min 后局部误差仍超出容限；请减小 h_min、放宽容差或改用适合刚性问题的方法。")  # Python 语法：raise 主动抛出异常；数值分析：避免返回不满足误差要求的结果。

t = np.asarray(t_values, dtype=float)  # Python 语法：把动态列表转为 NumPy 浮点数组；数值分析：t 保存所有被接受的非等距节点。
y = np.asarray(y_values, dtype=float)  # Python 语法：把解列表转为数组；标量问题得到一维数组，方程组问题得到二维数组。
hs = np.asarray(accepted_h, dtype=float)  # Python 语法：把步长列表转为数组；hs[i] 对应从 t[i] 到 t[i+1] 的步长。
errors = np.asarray(error_history, dtype=float)  # Python 语法：把误差列表转为数组；其元素应在区间 [0,1] 内。

print("方法：教材嵌套 RK2(3) = 变形 Euler（二阶） + 三阶 Runge-Kutta")  # Python 语法：print 向终端输出字符串；数值分析：说明当前模板采用的低阶与高阶组合。
print(f"接受步数={len(hs)}, 拒绝步数={rejected_steps}, 终点近似={y[-1]}")  # Python 语法：f 字符串把表达式嵌入文本；数值分析：汇总自适应计算量与终点结果。
for ti, yi in zip(t, y):  # Python 语法：zip 将两个序列按位置配对，for 逐对遍历；数值分析：输出每个自适应节点上的近似解。
    print(f"t={ti:.8f}, y={yi}")  # Python 语法：:.8f 指定浮点数显示 8 位小数；数值分析：自适应节点一般不是等距小数。

plt.step(t[:-1], hs, where="post", linewidth=1.3, label="accepted h")  # Python 语法：切片 [:-1] 排除末元素；数值分析：步长在每个接受区间上视作常数。
plt.xlabel("t")  # Python 语法：调用 xlabel 设置横轴标题；这里横轴表示时间节点。
plt.ylabel("accepted step size h")  # Python 语法：调用 ylabel 设置纵轴标题；这里纵轴表示被接受步长。
plt.grid(True)  # Python 语法：布尔值 True 开启网格；便于观察步长在哪些区间自动缩小或增大。
plt.legend()  # Python 语法：显示带 label 的图例；这里标明曲线表示被接受步长。
plt.tight_layout()  # Python 语法：自动调整图形边距；避免坐标轴文字被裁切。
plt.show()  # Python 语法：显示图窗；批量验证脚本会使用无窗口后端，因此不会阻塞。
