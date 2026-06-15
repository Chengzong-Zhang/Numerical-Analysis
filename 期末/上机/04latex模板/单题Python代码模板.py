"""单道数值 ODE 上机题模板：修改“必须替换”区域后直接运行。"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

OUTPUT_DIR = Path(__file__).resolve().parent

# ======================== 必须替换 ========================
f = lambda x, y: -y + x**2 + 1  # 方程右端 f(x, y)。
exact = lambda x: x**2 - 2 * x + 3 - 2 * np.exp(-x)  # 无精确解时改为 None。
x0, x_end, y0, h = 0.0, 1.0, 1.0, 0.1  # 区间、初值、步长。
# =========================================================


def rk4(f, x0, x_end, y0, h):
    """经典四阶 Runge--Kutta；标量与 NumPy 向量初值均可使用。"""
    n_steps = round((x_end - x0) / h)
    if not np.isclose(x0 + n_steps * h, x_end):
        raise ValueError("步长 h 必须整除积分区间长度。")

    x = x0 + np.arange(n_steps + 1) * h
    y0_array = np.asarray(y0, dtype=float)
    y = np.zeros((n_steps + 1,) + y0_array.shape, dtype=float)
    y[0] = y0_array

    for n in range(n_steps):
        k1 = f(x[n], y[n])
        k2 = f(x[n] + h / 2, y[n] + h * k1 / 2)
        k3 = f(x[n] + h / 2, y[n] + h * k2 / 2)
        k4 = f(x[n] + h, y[n] + h * k3)
        y[n + 1] = y[n] + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return x, y


def write_latex_table(x, y, y_exact=None):
    """生成可由 LaTeX 主文件直接 input 的三线表。"""
    lines = [
        r"\begin{table}[H]",
        r"  \centering",
        r"  \caption{各节点数值结果}",
    ]
    if y_exact is None:
        lines += [
            r"  \begin{tabular}{cc}",
            r"    \toprule",
            r"    $x_n$ & 数值解 $y_n$ \\",
            r"    \midrule",
        ]
        lines += [f"    {xn:.6f} & {yn:.10f} \\\\" for xn, yn in zip(x, y)]
    else:
        lines += [
            r"  \begin{tabular}{cccc}",
            r"    \toprule",
            r"    $x_n$ & 数值解 $y_n$ & 精确解 $y(x_n)$ & 绝对误差 \\",
            r"    \midrule",
        ]
        lines += [
            f"    {xn:.6f} & {yn:.10f} & {ye:.10f} & {abs(yn - ye):.3e} \\\\"
            for xn, yn, ye in zip(x, y, y_exact)
        ]
    lines += [r"    \bottomrule", r"  \end{tabular}", r"\end{table}"]
    (OUTPUT_DIR / "result_table.tex").write_text("\n".join(lines), encoding="utf-8")


def save_figures(x, y, y_exact=None):
    """生成供 LaTeX 自动插入的 PDF 图片。"""
    plt.plot(x, y, "ro--", label="numerical")
    if y_exact is not None:
        plt.plot(x, y_exact, "k-", label="exact")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "solution.pdf")
    plt.close()

    if y_exact is not None:
        plt.semilogy(x, np.maximum(np.abs(y - y_exact), np.finfo(float).eps), "bo-")
        plt.xlabel("x")
        plt.ylabel("absolute error")
        plt.grid(True, which="both")
        plt.tight_layout()
        plt.savefig(OUTPUT_DIR / "error.pdf")
        plt.close()


def main():
    x, y = rk4(f, x0, x_end, y0, h)
    if y.ndim != 1:
        raise ValueError("当前结果表按标量方程编写；方程组请选取要报告的分量。")

    y_exact = None if exact is None else exact(x)
    write_latex_table(x, y, y_exact)
    save_figures(x, y, y_exact)

    print("x            numerical")
    for xn, yn in zip(x, y):
        print(f"{xn:8.4f}     {yn:.10f}")
    if y_exact is not None:
        errors = np.abs(y - y_exact)
        print(f"maximum absolute error = {errors.max():.6e}")
        print(f"endpoint absolute error = {errors[-1]:.6e}")


if __name__ == "__main__":
    main()
