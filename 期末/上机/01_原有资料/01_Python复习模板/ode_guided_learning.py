r"""
数值 ODE 上机考试：边学、边读、边改、边运行的教学脚本

第 0 部分：我该怎么运行这个文件
================================
在 PowerShell 中运行：
    cd C:\coding\matlab\期末\上机
    python ode_guided_learning.py

正常现象：
1. 终端依次出现 Euler、RK2、RK4、收敛阶、RKF45、二阶系统、Adams 的短输出。
2. 当前文件夹生成 ode_guided_scalar.png 和 ode_guided_system.png。
3. 收敛阶大致为 Euler=1、改进 Euler=2、RK4=4。

第一遍先不用懂：
- _prepare_y、_call_f、_fixed_step_solver 是为了同时兼容标量和向量的工具。
- RKF45 的分数系数直接从模板复制，先只看“算误差 -> 接受/拒绝 -> 调步长”。
- Adams 系数先不推导，只记住“RK4 起步 + PECE”。

本文件基于已有 ode_exam_template.py、第五章教材及老师 problem1/2/3.m 改造。
仅使用 numpy、matplotlib，不使用 scipy。
"""

import numpy as np
import matplotlib.pyplot as plt


# ============================================================================
# 第 1 部分：ODE 上机题的统一套路
# ============================================================================
#
# 原问题永远先认成：
#     y' = f(t, y),    y(t0) = y0,    从 t0 算到 T
#
# 考场最重要的输入只有：
#     f       右端函数，最常改
#     t0      初始时刻
#     y0      初值
#     T       终点
#     h/tol   固定步长 h，或自适应误差容限 tol
#
# 标量 y：例如 y'=t+y，y0=1，代码中 y 是一个数。
# 向量 y：例如方程组，y0=np.array([1, 0])，代码中 y 是数组。
#
# 二阶方程不能直接喂给一阶方法。先令 w1=y, w2=y'，化成一阶向量系统。
#
# 固定步长：一直迈同样大的 h，代码简单。
# 自适应步长：误差大就缩步，误差小就放步，实际步长每次可能不同。


def _prepare_y(y0):
    """内部工具：统一转成一维数组计算，并记住原问题是不是标量。"""
    scalar = np.asarray(y0).ndim == 0
    y = np.atleast_1d(np.asarray(y0, dtype=float)).copy()
    return y, scalar


def _call_f(f, t, y, scalar):
    """内部工具：标量问题把数组 y=[value] 临时还原成数字再调用 f。"""
    argument = float(y[0]) if scalar else y
    value = np.atleast_1d(np.asarray(f(t, argument), dtype=float))
    if value.shape != y.shape:
        raise ValueError(f"f 返回形状 {value.shape}，但 y 的形状是 {y.shape}")
    return value


def _finish(ts, ys, scalar):
    """内部工具：标量问题返回一维 ys；向量问题返回二维 ys。"""
    ts = np.asarray(ts, dtype=float)
    ys = np.asarray(ys, dtype=float)
    return ts, ys[:, 0] if scalar else ys


def _fixed_step_solver(f, t0, y0, T, h, one_step):
    """固定步长公共循环；one_step 决定这一小步用 Euler、RK2 还是 RK4。"""
    if h <= 0 or T < t0:
        raise ValueError("要求 h>0 且 T>=t0")
    y, scalar = _prepare_y(y0)
    t = float(t0)
    ts, ys = [t], [y.copy()]

    while t < T - 1e-14:
        hn = min(h, T - t)               # 最后一步自动截到 T
        y = one_step(f, t, y, hn, scalar)
        t += hn
        ts.append(t)
        ys.append(y.copy())
    return _finish(ts, ys, scalar)


# ============================================================================
# 第 2 部分：最简单的 Euler 法
# ============================================================================
#
# 数学只看三行：
#     当前斜率：      k = f(t_n, y_n)
#     向前走一步：    y_{n+1} = y_n + h*k
#     新时刻：        t_{n+1} = t_n + h
#
# 意义：把当前位置的切线斜率当成整小步都不变。简单，但只有一阶。


def euler_simple(f, t0, y0, T, h):
    """只支持标量的极短 Euler，先读懂这个。"""
    ts = np.arange(t0, T + h / 2, h)      # 建立等距时间点
    ys = np.zeros(len(ts))                # 准备存数值解
    ys[0] = y0                            # 放入初值
    for n in range(len(ts) - 1):          # 从第 0 步一直算到最后一步前
        ys[n + 1] = ys[n] + h * f(ts[n], ys[n])  # Euler 更新公式
    return ts, ys


def forward_euler(f, t0, y0, T, h):
    """通用向前 Euler：标量和向量都支持。考试中需要 Euler 时调用它。"""
    def one_step(f, t, y, h, scalar):
        slope = _call_f(f, t, y, scalar)  # 当前点斜率 f(t_n,y_n)
        return y + h * slope              # y_{n+1}=y_n+h*f_n

    return _fixed_step_solver(f, t0, y0, T, h, one_step)


# 小练习：把下面 demo 中的 f 改成 lambda t, y: t + y，再运行一次。


# ============================================================================
# 第 3 部分：从 Euler 到改进 Euler / RK2
# ============================================================================
#
# 改进 Euler 的人话：
# 1. 先用当前斜率 k1 预测终点位置；
# 2. 在预测终点算斜率 k2；
# 3. 用起点和终点斜率的平均值校正。
#
# k1：在当前点看到的斜率。
# k2：在试探点/预测点看到的斜率。


def improved_euler(f, t0, y0, T, h):
    """改进 Euler / 显式梯形：先预测终点，再平均校正，二阶。"""
    def one_step(f, t, y, h, scalar):
        k1 = _call_f(f, t, y, scalar)                  # 起点斜率
        y_predict = y + h * k1                         # Euler 预测终点
        k2 = _call_f(f, t + h, y_predict, scalar)      # 预测终点斜率
        return y + h * (k1 + k2) / 2                   # 两个斜率取平均

    return _fixed_step_solver(f, t0, y0, T, h, one_step)


def midpoint_rk2(f, t0, y0, T, h):
    """中点 RK2 / 变形 Euler：用中点斜率走完整一步，二阶。"""
    def one_step(f, t, y, h, scalar):
        k1 = _call_f(f, t, y, scalar)                  # 起点斜率
        k2 = _call_f(f, t + h / 2, y + h * k1 / 2, scalar)  # 中点斜率
        return y + h * k2                              # 用中点斜率更新

    return _fixed_step_solver(f, t0, y0, T, h, one_step)


def heun_rk2(f, t0, y0, T, h):
    """教材/老师 problem1 使用的 Heun RK2：a2=2/3，二阶。"""
    def one_step(f, t, y, h, scalar):
        k1 = _call_f(f, t, y, scalar)                       # 起点斜率
        k2 = _call_f(f, t + 2 * h / 3, y + 2 * h * k1 / 3, scalar)
        return y + h * (k1 + 3 * k2) / 4                    # 权重 1/4, 3/4

    return _fixed_step_solver(f, t0, y0, T, h, one_step)


# 考场用途：
# improved_euler：题目写“改进 Euler / 显式梯形 / 预测校正”时用。
# midpoint_rk2：题目写“中点方法 / 变形 Euler”时用。
# heun_rk2：老师题目明确写 Heun 时用，注意本教材权重是 1/4 和 3/4。
# 小练习：运行第 5 部分，观察三种 RK2 的阶数都应接近 2。


# ============================================================================
# 第 4 部分：经典 RK4
# ============================================================================
#
# RK4 是“四个斜率加权平均”：
#     k1 = f(t_n,       y_n)
#     k2 = f(t_n+h/2,   y_n+h*k1/2)
#     k3 = f(t_n+h/2,   y_n+h*k2/2)
#     k4 = f(t_n+h,     y_n+h*k3)
#     y_{n+1} = y_n + h*(k1+2*k2+2*k3+k4)/6


def rk4(f, t0, y0, T, h):
    """经典四阶 Runge-Kutta；通常是考场最稳的默认方法。"""
    def one_step(f, t, y, h, scalar):
        k1 = _call_f(f, t, y, scalar)                         # 斜率 1：起点
        k2 = _call_f(f, t + h / 2, y + h * k1 / 2, scalar)   # 斜率 2：半步
        k3 = _call_f(f, t + h / 2, y + h * k2 / 2, scalar)   # 斜率 3：再算半步
        k4 = _call_f(f, t + h, y + h * k3, scalar)           # 斜率 4：终点
        return y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6       # 加权平均更新

    return _fixed_step_solver(f, t0, y0, T, h, one_step)


# ===== 考场最常改的 5 个位置 =====
# f  = lambda t, y: ...     # 1. 抄题目右端
# t0 = ...                  # 2. 初始时刻
# y0 = ...                  # 3. 初值；方程组写 np.array([...])
# T  = ...                  # 4. 终点
# h  = ...                  # 5. 固定步长
# ts, ys = rk4(f, t0, y0, T, h)


# ============================================================================
# 第 5 部分：误差与收敛阶
# ============================================================================
#
# 最大误差 max_error：所有网格点误差中最大的一个。
# 用 h,h/2,h/4：若 E(h)≈C*h^p，步长减半后误差约缩小 2^p 倍。
#     p = log(E(h)/E(h/2)) / log(2)


def max_error(ts, ys, exact, component=None):
    """已知精确解时算最大误差；向量问题可用 component=0 只检查 y。"""
    numerical = np.asarray(ys, dtype=float)
    exact_values = np.asarray([exact(t) for t in ts], dtype=float)
    if component is not None:
        numerical = numerical[:, component]
        if exact_values.ndim > 1:
            exact_values = exact_values[:, component]
    difference = numerical - exact_values
    if difference.ndim == 1:
        return float(np.max(np.abs(difference)))
    return float(np.max(np.linalg.norm(difference, ord=np.inf, axis=1)))


def estimate_order(method, f, t0, y0, T, h, exact, component=None):
    """用 h,h/2,h/4 计算误差和两次观测阶。"""
    hs = np.array([h, h / 2, h / 4], dtype=float)
    errors = []
    for hi in hs:
        ts, ys = method(f, t0, y0, T, hi)
        errors.append(max_error(ts, ys, exact, component))
    errors = np.asarray(errors)
    orders = np.log(errors[:-1] / errors[1:]) / np.log(2)  # 收敛阶公式
    return hs, errors, orders


# 若阶数不接近理论值，先检查：
# 1. f、精确解或初值抄错；2. RK 公式系数错；3. 误差取错分量；
# 4. h 还不够小；5. h 太小导致舍入误差主导；6. 显式法因刚性而不稳定。


# ============================================================================
# 第 6 部分：自适应 RKF45，重点
# ============================================================================
#
# 固定步长：一直迈同样大的步。
# 自适应：误差大就缩步，误差小就放步。
#
# K1,...,K6 是同一小步内六个不同试探位置的斜率，四阶和五阶公式共享它们。
# 同时算 y4 和 y5，是为了不用精确解也能用二者差值估计本步误差：
#     err = norm(y5-y4)
#
# 最关键控制流：
#     err <= tol：接受步，才允许更新 t,y，并保存结果。
#     err >  tol：拒绝步，不许更新 t,y，只缩小 h 后重算。


def rkf45_adaptive(f, t0, y0, T, h0, tol, h_min=1e-8, h_max=None):
    """完整 RKF45。返回 ts, ys, hs；hs 只记录每次接受步实际使用的步长。"""
    if h0 <= 0 or tol <= 0 or h_min <= 0 or T < t0:
        raise ValueError("要求 h0,tol,h_min>0 且 T>=t0")
    if h_max is None:
        h_max = max(h0, T - t0)

    y, scalar = _prepare_y(y0)
    t, h = float(t0), min(h0, h_max)
    ts, ys, hs = [t], [y.copy()], []

    while t < T - 1e-14:
        h = min(h, h_max, T - t)  # 最后一步不能越过 T

        # 六个 K 都是斜率 f(...)；系数直接复制，不要在考场临时推导。
        K1 = _call_f(f, t, y, scalar)
        K2 = _call_f(f, t + h / 4, y + h * K1 / 4, scalar)
        K3 = _call_f(f, t + 3 * h / 8,
                     y + h * (3 * K1 / 32 + 9 * K2 / 32), scalar)
        K4 = _call_f(f, t + 12 * h / 13,
                     y + h * (1932 * K1 / 2197 - 7200 * K2 / 2197
                              + 7296 * K3 / 2197), scalar)
        K5 = _call_f(f, t + h,
                     y + h * (439 * K1 / 216 - 8 * K2 + 3680 * K3 / 513
                              - 845 * K4 / 4104), scalar)
        K6 = _call_f(f, t + h / 2,
                     y + h * (-8 * K1 / 27 + 2 * K2 - 3544 * K3 / 2565
                              + 1859 * K4 / 4104 - 11 * K5 / 40), scalar)

        # 用共享的 K 同时得到四阶值和五阶值。
        y4 = y + h * (25 * K1 / 216 + 1408 * K3 / 2565
                      + 2197 * K4 / 4104 - K5 / 5)
        y5 = y + h * (16 * K1 / 135 + 6656 * K3 / 12825
                      + 28561 * K4 / 56430 - 9 * K5 / 50 + 2 * K6 / 55)
        err = np.linalg.norm(y5 - y4, ord=np.inf)  # 向量也能用的无穷范数误差

        if err <= tol:
            # 接受步：此时才更新 t,y，并把本次实际步长 h 保存进 hs。
            t += h
            y = y5
            ts.append(t)
            ys.append(y.copy())
            hs.append(h)
        # 若 err>tol：这里故意没有 else 更新；拒绝步只允许在下方改变 h。

        # 按误差调整下一次尝试步长；0.84 是保险系数，0.1~4 防止变化太猛。
        factor = 4.0 if err == 0 else 0.84 * (tol / err) ** 0.25
        factor = min(4.0, max(0.1, factor))
        new_h = min(h_max, h * factor)
        if err > tol and new_h < h_min:
            raise RuntimeError("无法满足 tol：所需步长小于 h_min")
        h = max(h_min, new_h)

    ts, ys = _finish(ts, ys, scalar)
    return ts, ys, np.asarray(hs)


# 自适应算法考场最容易错的 8 点：
# 1. 拒绝步仍更新 t；2. 拒绝步仍更新 y；3. 把拒绝步 h 存进 hs；
# 4. 忘记最后一步截到 T；5. y4/y5 系数抄错；6. 向量误差没取范数；
# 7. err=0 时除零；8. 没有 h_min，失败时死循环。


# ============================================================================
# 第 7 部分：二阶方程化为一阶系统
# ============================================================================
#
# 原方程：y'' = g(t,y,y')
# 令：    w1=y, w2=y'
# 则：    w1'=w2, w2'=g(t,w1,w2)
#
# 例：y''+y=0，所以 y''=-y：


def harmonic_system(t, w):
    """y''+y=0 的一阶系统。"""
    return np.array([w[1], -w[0]])  # [w1',w2']=[w2,-w1]


# 用法：
# ts, ys = rk4(harmonic_system, 0, np.array([1, 0]), 2*np.pi, 0.1)
# ys[:,0] 是每个时刻的 w1，也就是 y。
# ys[:,1] 是每个时刻的 w2，也就是 y'。
#
# ===== 考场换题模板 =====
# def system_f(t, w):
#     y, dy = w
#     ddy = ...                 # 把题目 y''=g(t,y,y') 抄到这里
#     return np.array([dy, ddy])
# y0_system = np.array([题目给的y(t0), 题目给的y'(t0)])


# ============================================================================
# 第 8 部分：Adams PECE，只做考场够用版
# ============================================================================
#
# 多步法计算 y_{n+1} 时要用多个旧点；题目只给 y0，所以先用同为四阶的 RK4
# 算出 y1,y2,y3。随后：
# AB4 预测：
# yP = yn + h/24*(55fn-59f(n-1)+37f(n-2)-9f(n-3))
# AM4 校正：
# y(n+1) = yn + h/24*(9fP+19fn-5f(n-1)+f(n-2))


def adams_pece4(f, t0, y0, T, h):
    """四阶 Adams PECE；只用于固定等步长且 (T-t0)/h 为整数。"""
    n_float = (T - t0) / h
    n_steps = int(round(n_float))
    if h <= 0 or n_steps < 3 or not np.isclose(n_float, n_steps):
        raise ValueError("Adams 要求固定等步长、至少三步、终点正好落在网格上")

    y0_array, scalar = _prepare_y(y0)
    ts = t0 + h * np.arange(n_steps + 1)
    ys = np.zeros((n_steps + 1, y0_array.size))
    ys[0] = y0_array

    # RK4 起步：先得到 y1,y2,y3。
    for i in range(3):
        start_y = float(ys[i, 0]) if scalar else ys[i]
        _, small_ys = rk4(f, ts[i], start_y, ts[i + 1], h)
        ys[i + 1] = np.atleast_1d(small_ys[-1])

    fs = np.zeros_like(ys)
    for i in range(4):
        fs[i] = _call_f(f, ts[i], ys[i], scalar)

    for n in range(3, n_steps):
        # P = Predict：AB4 预报 yP。
        y_predict = ys[n] + h * (
            55 * fs[n] - 59 * fs[n - 1] + 37 * fs[n - 2] - 9 * fs[n - 3]
        ) / 24
        # E = Evaluate：计算预报点的 fP。
        f_predict = _call_f(f, ts[n + 1], y_predict, scalar)
        # C = Correct：AM4 用 fP 校正。
        ys[n + 1] = ys[n] + h * (
            9 * f_predict + 19 * fs[n] - 5 * fs[n - 1] + fs[n - 2]
        ) / 24
        # E = Evaluate：保存校正后真正的 f(n+1)，供下一步使用。
        fs[n + 1] = _call_f(f, ts[n + 1], ys[n + 1], scalar)

    return _finish(ts, ys, scalar)


# ============================================================================
# 第 9 部分：老师上机题型对应模板
# ============================================================================
#
# 题型 1（老师 problem1）五种方法比较 + 画图：
#     调 forward_euler / improved_euler / heun_rk2 / midpoint_rk2 / rk4，
#     再用 plot_results；第二问调 adams_pece4 与 rk4 比较。
#
# 题型 2（老师 problem2）二阶化系统 + RK4 + 收敛阶：
#     写 system_f 和向量 y0；调 rk4；误差只看 y 时用 component=0；
#     调 estimate_order(rk4,...,component=0)。
#
# 题型 3（老师 problem3）刚性方程组 + Euler/梯形/RK4：
#     最低限度写 f(t,y)=A@y+g(t)，Euler 和 RK4 可直接调用；
#     线性隐式梯形每步解：
#     (I-h*A/2)y_{n+1}=(I+h*A/2)y_n+h*(g(t_n)+g(t_{n+1}))/2
#     Python 用 np.linalg.solve(left, right)。显式法步长过大会发散，不等于代码错。
#
# 题型 4 自适应 RKF45：
#     ts, ys, hs = rkf45_adaptive(f,t0,y0,T,h0,tol)
#     ts/ys 只含接受点；hs 只含接受步，且 len(hs)==len(ts)-1。


def plot_results(results, exact, filename, title, component=None):
    """画多种方法与真解；component 用于向量解。"""
    plt.figure(figsize=(8, 5))
    for name, (ts, ys) in results.items():
        values = ys if np.asarray(ys).ndim == 1 else ys[:, component or 0]
        plt.plot(ts, values, marker="o", markersize=3, label=name)
    dense_t = np.linspace(min(v[0][0] for v in results.values()),
                          max(v[0][-1] for v in results.values()), 500)
    exact_values = np.asarray([exact(t) for t in dense_t])
    if exact_values.ndim > 1:
        exact_values = exact_values[:, component or 0]
    plt.plot(dense_t, exact_values, "k-", linewidth=2, label="exact")
    plt.xlabel("t")
    plt.ylabel("y")
    plt.title(title)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename, dpi=160)
    plt.close()


# ============================================================================
# 第 10 部分：主函数 main
# ============================================================================
# 每个 demo 都是独立函数。复习时可在 main 中注释掉已掌握部分，只运行当前部分。


def demo_euler_and_rk2():
    print("\n[2-3] Euler 与 RK2：练习方程 y'=t+y, y(0)=1")
    f = lambda t, y: t + y
    for name, method in [("Euler", euler_simple), ("Improved", improved_euler),
                         ("Midpoint", midpoint_rk2), ("Heun", heun_rk2)]:
        _, ys = method(f, 0.0, 1.0, 1.0, 0.1)
        print(f"{name:<10} y(1)={ys[-1]:.8f}")


def demo_rk4_and_order():
    print("\n[4-5] RK4 最小例子 + 看输出判断阶数")
    f = lambda t, y: -y + t**2 + 1
    exact = lambda t: t**2 - 2 * t + 3 - 2 * np.exp(-t)
    methods = [("Euler", forward_euler), ("Improved", improved_euler), ("RK4", rk4)]
    results = {}
    for name, method in methods:
        hs, errors, orders = estimate_order(method, f, 0.0, 1.0, 1.0, 0.2, exact)
        print(f"{name:<10} errors={np.array2string(errors, precision=2)} "
              f"orders={np.array2string(orders, precision=2)}")
    # 老师 problem1 的画图题：五种固定步长方法一起比较。
    for name, method in [("Euler", forward_euler), ("Improved", improved_euler),
                         ("Midpoint", midpoint_rk2), ("Heun", heun_rk2), ("RK4", rk4)]:
        results[name] = method(f, 0.0, 1.0, 1.0, 0.1)
    plot_results(results, exact, "ode_guided_scalar.png", "Fixed-step methods")


def demo_rkf45():
    print("\n[6] RKF45：tol 越小，通常接受步数越多")
    f = lambda t, y: -y + t**2 + 1
    for tol in [1e-4, 1e-6, 1e-8]:
        ts, ys, hs = rkf45_adaptive(f, 0.0, 1.0, 1.0, 0.1, tol, h_max=0.5)
        print(f"tol={tol:.0e}: 接受 {len(hs):2d} 步, hs={np.array2string(hs, precision=4)}")
        assert len(hs) == len(ts) - 1


def demo_second_order_system():
    print("\n[7] 二阶方程 y''+y=0 化系统后用 RK4")
    exact = lambda t: np.array([np.cos(t), -np.sin(t)])
    ts, ys = rk4(harmonic_system, 0.0, np.array([1.0, 0.0]), 2 * np.pi, 0.1)
    print(f"ys[:,0] 是 y；ys[:,1] 是 y'；最大误差={max_error(ts, ys, exact):.3e}")
    plot_results({"RK4": (ts, ys)}, exact, "ode_guided_system.png",
                 "Second-order ODE as a system", component=0)


def demo_adams():
    print("\n[8] Adams PECE 与 RK4 比较")
    f = lambda t, y: -y + t**2 + 1
    ts_rk, ys_rk = rk4(f, 0.0, 1.0, 1.0, 0.1)
    ts_ad, ys_ad = adams_pece4(f, 0.0, 1.0, 1.0, 0.1)
    print(f"终点差 |Adams-RK4|={abs(ys_ad[-1]-ys_rk[-1]):.3e}; "
          f"网格一致={np.allclose(ts_ad, ts_rk)}")


def main():
    print("=== 数值 ODE 教学型模板：按第 2 -> 8 部分运行 ===")
    print("[0-1] 命令：python ode_guided_learning.py；统一输入是 f,t0,y0,T,h/tol")
    demo_euler_and_rk2()       # 已掌握后可注释
    demo_rk4_and_order()       # 已掌握后可注释
    demo_rkf45()               # 自适应重点
    demo_second_order_system() # 二阶化系统
    demo_adams()               # 多步法最低限度
    print("\n[9] 老师题型映射见源码第 9 部分和 ode_guided_notes.md")
    print("\n完成：已生成 ode_guided_scalar.png、ode_guided_system.png")


if __name__ == "__main__":
    main()
