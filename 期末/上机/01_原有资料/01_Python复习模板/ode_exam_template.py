"""
数值分析第五章 ODE 上机考试模板

依赖：numpy、matplotlib（不使用 scipy）
直接运行：
    python ode_exam_template.py

约定：
- 固定步长方法返回 ts, ys。
- 自适应 RKF45 返回 ts, ys, hs，其中 hs[i] 是第 i 个接受步的实际步长。
- 标量初值返回一维 ys；向量初值返回形状 (步数+1, 维数) 的 ys。
"""

import numpy as np
import matplotlib.pyplot as plt


def _prepare_y(y0):
    """内部工具：统一按一维数组计算，并记录原问题是否为标量。"""
    scalar = np.asarray(y0).ndim == 0
    return np.atleast_1d(np.asarray(y0, dtype=float)).copy(), scalar


def _call_f(f, t, y, scalar):
    """内部工具：让 f(t,y) 同时支持标量和向量。"""
    arg = float(y[0]) if scalar else y
    value = np.atleast_1d(np.asarray(f(t, arg), dtype=float))
    if value.shape != y.shape:
        raise ValueError(f"f(t,y) 返回形状 {value.shape}，但 y 的形状是 {y.shape}")
    return value


def _finish(ts, ys, scalar):
    ts = np.asarray(ts, dtype=float)
    ys = np.asarray(ys, dtype=float)
    return ts, ys[:, 0] if scalar else ys


def _fixed_step_solver(f, t0, y0, T, h, step):
    """固定步长方法的公共循环；最后一步自动缩短到 T。"""
    if h <= 0 or T < t0:
        raise ValueError("要求 h > 0 且 T >= t0")
    y, scalar = _prepare_y(y0)
    t = float(t0)
    ts, ys = [t], [y.copy()]
    eps = 10 * np.finfo(float).eps * max(1.0, abs(T))

    while t < T - eps:
        hn = min(h, T - t)
        y = step(f, t, y, hn, scalar)
        t += hn
        ts.append(t)
        ys.append(y.copy())
    return _finish(ts, ys, scalar)


def forward_euler(f, t0, y0, T, h):
    """向前 Euler：一阶。"""
    def step(f, t, y, h, scalar):
        return y + h * _call_f(f, t, y, scalar)

    return _fixed_step_solver(f, t0, y0, T, h, step)


def improved_euler(f, t0, y0, T, h):
    """改进 Euler / 显式梯形 / Euler 预测校正：二阶。"""
    def step(f, t, y, h, scalar):
        k1 = _call_f(f, t, y, scalar)
        k2 = _call_f(f, t + h, y + h * k1, scalar)
        return y + h * (k1 + k2) / 2

    return _fixed_step_solver(f, t0, y0, T, h, step)


def midpoint_rk2(f, t0, y0, T, h):
    """中点 RK2 / 变形 Euler：二阶。"""
    def step(f, t, y, h, scalar):
        k1 = _call_f(f, t, y, scalar)
        k2 = _call_f(f, t + h / 2, y + h * k1 / 2, scalar)
        return y + h * k2

    return _fixed_step_solver(f, t0, y0, T, h, step)


def heun_rk2(f, t0, y0, T, h):
    """教材中的二阶 Heun：a2=2/3，权重为 1/4 和 3/4。"""
    def step(f, t, y, h, scalar):
        k1 = _call_f(f, t, y, scalar)
        k2 = _call_f(f, t + 2 * h / 3, y + 2 * h * k1 / 3, scalar)
        return y + h * (k1 + 3 * k2) / 4

    return _fixed_step_solver(f, t0, y0, T, h, step)


def rk4(f, t0, y0, T, h):
    """经典四阶 Runge-Kutta。"""
    def step(f, t, y, h, scalar):
        k1 = _call_f(f, t, y, scalar)
        k2 = _call_f(f, t + h / 2, y + h * k1 / 2, scalar)
        k3 = _call_f(f, t + h / 2, y + h * k2 / 2, scalar)
        k4 = _call_f(f, t + h, y + h * k3, scalar)
        return y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6

    return _fixed_step_solver(f, t0, y0, T, h, step)


def rkf45_adaptive(f, t0, y0, T, h0, tol, h_min=1e-6, h_max=None):
    """
    自适应 Runge-Kutta-Fehlberg 4(5)。

    返回 ts, ys, hs；hs 只记录接受步的实际步长。
    tol 是每一步的绝对局部误差限，向量误差使用无穷范数。
    """
    if h0 <= 0 or tol <= 0 or h_min <= 0 or T < t0:
        raise ValueError("要求 h0,tol,h_min > 0 且 T >= t0")
    if h_max is None:
        h_max = max(h0, T - t0)
    if h_max < h_min:
        raise ValueError("要求 h_max >= h_min")

    y, scalar = _prepare_y(y0)
    t, h = float(t0), min(h0, h_max)
    ts, ys, hs = [t], [y.copy()], []
    eps = 10 * np.finfo(float).eps * max(1.0, abs(T))

    while t < T - eps:
        h = min(h, h_max, T - t)  # 保证最后一步不超过 T

        K1 = _call_f(f, t, y, scalar)
        K2 = _call_f(f, t + h / 4, y + h * K1 / 4, scalar)
        K3 = _call_f(
            f, t + 3 * h / 8,
            y + h * (3 * K1 / 32 + 9 * K2 / 32), scalar
        )
        K4 = _call_f(
            f, t + 12 * h / 13,
            y + h * (1932 * K1 / 2197 - 7200 * K2 / 2197 + 7296 * K3 / 2197),
            scalar
        )
        K5 = _call_f(
            f, t + h,
            y + h * (439 * K1 / 216 - 8 * K2 + 3680 * K3 / 513 - 845 * K4 / 4104),
            scalar
        )
        K6 = _call_f(
            f, t + h / 2,
            y + h * (-8 * K1 / 27 + 2 * K2 - 3544 * K3 / 2565
                     + 1859 * K4 / 4104 - 11 * K5 / 40),
            scalar
        )

        y4 = y + h * (25 * K1 / 216 + 1408 * K3 / 2565
                      + 2197 * K4 / 4104 - K5 / 5)
        y5 = y + h * (16 * K1 / 135 + 6656 * K3 / 12825
                      + 28561 * K4 / 56430 - 9 * K5 / 50 + 2 * K6 / 55)
        err = np.linalg.norm(y5 - y4, ord=np.inf)

        if err <= tol:
            # 接受步：用五阶值更新 t,y，并保存实际使用的步长
            t += h
            y = y5
            ts.append(t)
            ys.append(y.copy())
            hs.append(h)
        # 拒绝步：不更新时间和解，只在下面缩小步长

        if err == 0:
            factor = 4.0
        else:
            factor = 0.84 * (tol / err) ** 0.25
            factor = min(4.0, max(0.1, factor))
        new_h = min(h_max, h * factor)

        if err > tol and new_h < h_min:
            raise RuntimeError(
                f"在 t={t:.6g} 无法满足 tol={tol:g}；所需步长小于 h_min={h_min:g}"
            )
        h = max(h_min, new_h)

    ts, ys = _finish(ts, ys, scalar)
    return ts, ys, np.asarray(hs, dtype=float)


def adams_pece4(f, t0, y0, T, h):
    """
    四阶 Adams PECE：RK4 提供 y1,y2,y3，AB4 预报，AM4 校正一次。
    仅适用于固定步长，且 (T-t0)/h 应为整数。
    """
    n_float = (T - t0) / h
    n_steps = int(round(n_float))
    if h <= 0 or n_steps < 3 or not np.isclose(n_float, n_steps):
        raise ValueError("Adams PECE 要求 h>0、至少三步，且 (T-t0)/h 为整数")

    y0_arr, scalar = _prepare_y(y0)
    ts = t0 + h * np.arange(n_steps + 1)
    ys = np.zeros((n_steps + 1, y0_arr.size), dtype=float)
    ys[0] = y0_arr

    # RK4 起步，得到 y1,y2,y3
    for i in range(3):
        _, start_ys = rk4(f, ts[i], float(ys[i, 0]) if scalar else ys[i], ts[i + 1], h)
        ys[i + 1] = np.atleast_1d(start_ys[-1])

    fs = np.zeros_like(ys)
    for i in range(4):
        fs[i] = _call_f(f, ts[i], ys[i], scalar)

    for n in range(3, n_steps):
        y_pred = ys[n] + h * (
            55 * fs[n] - 59 * fs[n - 1] + 37 * fs[n - 2] - 9 * fs[n - 3]
        ) / 24
        f_pred = _call_f(f, ts[n + 1], y_pred, scalar)
        ys[n + 1] = ys[n] + h * (
            9 * f_pred + 19 * fs[n] - 5 * fs[n - 1] + fs[n - 2]
        ) / 24
        fs[n + 1] = _call_f(f, ts[n + 1], ys[n + 1], scalar)

    return _finish(ts, ys, scalar)


def max_error(ts, ys, exact, component=None):
    """已知精确解时计算最大误差；component 可指定向量解的某个分量。"""
    exact_values = np.asarray([exact(t) for t in ts], dtype=float)
    numerical = np.asarray(ys, dtype=float)
    if component is not None:
        numerical = numerical[:, component]
        exact_values = exact_values[:, component] if exact_values.ndim > 1 else exact_values
    diff = numerical - exact_values
    if diff.ndim == 1:
        return float(np.max(np.abs(diff)))
    return float(np.max(np.linalg.norm(diff, ord=np.inf, axis=1)))


def estimate_order(method, f, t0, y0, T, h, exact, component=None):
    """用 h,h/2,h/4 的最大误差估计整体收敛阶。"""
    hs = np.asarray([h, h / 2, h / 4], dtype=float)
    errors = []
    for hi in hs:
        ts, ys = method(f, t0, y0, T, hi)
        errors.append(max_error(ts, ys, exact, component))
    errors = np.asarray(errors)
    orders = np.log(errors[:-1] / errors[1:]) / np.log(2)
    return hs, errors, orders


def plot_solutions(results, exact=None, title="ODE 数值解比较", component=0, filename=None):
    """results 是 {方法名: (ts,ys)}；向量解默认画第 component 个分量。"""
    plt.figure(figsize=(8, 5))
    all_ts = []
    for name, (ts, ys) in results.items():
        values = ys if np.asarray(ys).ndim == 1 else ys[:, component]
        plt.plot(ts, values, marker="o", markersize=3, label=name)
        all_ts.extend(ts)
    if exact is not None:
        grid = np.linspace(min(all_ts), max(all_ts), 600)
        exact_values = np.asarray([exact(t) for t in grid])
        if exact_values.ndim > 1:
            exact_values = exact_values[:, component]
        plt.plot(grid, exact_values, "k-", linewidth=2, label="exact")
    plt.xlabel("t")
    plt.ylabel(f"y[{component}]" if any(np.asarray(v[1]).ndim > 1 for v in results.values()) else "y")
    plt.title(title)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    if filename:
        plt.savefig(filename, dpi=160)
    plt.close()


def solve_and_compare(f, t0, y0, T, h, exact=None, component=None):
    """运行五种固定步长方法，打印终值和最大误差。"""
    methods = {
        "Euler": forward_euler,
        "Improved Euler": improved_euler,
        "Midpoint RK2": midpoint_rk2,
        "Heun RK2": heun_rk2,
        "RK4": rk4,
    }
    results = {}
    print(f"{'method':<18} {'y(T)':>18} {'max error':>14}")
    print("-" * 54)
    for name, method in methods.items():
        ts, ys = method(f, t0, y0, T, h)
        results[name] = (ts, ys)
        final_value = ys[-1] if np.asarray(ys).ndim == 1 else ys[-1, component or 0]
        err_text = "---" if exact is None else f"{max_error(ts, ys, exact, component):.6e}"
        print(f"{name:<18} {final_value:>18.10g} {err_text:>14}")
    return results


def demo_scalar():
    """例1：标量方程 y'=-y+t^2+1，精确解 y=t^2-2t+3-2e^{-t}。"""
    print("\n=== 例1：标量方程 ===")
    f = lambda t, y: -y + t**2 + 1
    exact = lambda t: t**2 - 2 * t + 3 - 2 * np.exp(-t)
    results = solve_and_compare(f, 0.0, 1.0, 1.0, 0.1, exact)
    plot_solutions(results, exact, "Scalar ODE comparison", filename="ode_demo_scalar.png")

    hs, errors, orders = estimate_order(rk4, f, 0.0, 1.0, 1.0, 0.2, exact)
    print("\nRK4 收敛阶：")
    for i, (hi, err) in enumerate(zip(hs, errors)):
        p = "---" if i == 0 else f"{orders[i - 1]:.4f}"
        print(f"h={hi:<8g} max_error={err:.6e}  p={p}")

    ts, ys, used_hs = rkf45_adaptive(f, 0.0, 1.0, 1.0, 0.1, 1e-7)
    print(f"\nRKF45：接受 {len(used_hs)} 步，最大误差 {max_error(ts, ys, exact):.6e}")
    print("实际步长：", np.array2string(used_hs, precision=5))


def demo_system():
    """例2：y''+y=0 化为 w1'=w2, w2'=-w1。"""
    print("\n=== 例2：二阶方程化一阶系统 ===")
    f = lambda t, w: np.array([w[1], -w[0]])
    exact = lambda t: np.array([np.cos(t), -np.sin(t)])
    ts, ys = rk4(f, 0.0, np.array([1.0, 0.0]), 2 * np.pi, 0.1)
    err = max_error(ts, ys, exact)
    print(f"RK4 向量解最大无穷范数误差：{err:.6e}")
    plot_solutions({"RK4": (ts, ys)}, exact, "y'' + y = 0", component=0,
                   filename="ode_demo_system.png")


if __name__ == "__main__":
    demo_scalar()
    demo_system()
    print("\n已生成 ode_demo_scalar.png 和 ode_demo_system.png")
