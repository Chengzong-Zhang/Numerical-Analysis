# -*- coding: utf-8 -*-  # Python 语法：声明源文件采用 UTF-8 编码，使中文注释能够被正确读取。

# =============================================================================
# 211870125.pdf 对应问题、数学推导与经典四阶 Runge-Kutta（RK4）算法
# =============================================================================
#
# 【PDF 中的问题】
# 使用经典四阶 Runge-Kutta 方法数值求解二阶两点边值问题
#
#     y''(x) + (2/x)y'(x) - (6/x^2)y(x) = 5 - 6x + 7x^2,    1 < x < 2,
#     y(1) = 1/2,
#     y(2) = 4 + 4 ln(2).
#
# 题目给出的真解为
#
#     y(x) = x^2 - x^3 + (1/2)x^4 + x^2 ln(x).
#
# 分别取步长
#
#     h = 1/200, 1/100, 1/50, 1/25, 1/10, 1/5, 1/4,
#
# 计算 x=2 处的绝对误差，并根据数值结果估计误差收敛阶、分析结果。
#
# 【为什么可以把本题化成一阶初值问题】
# 令 z(x) = y'(x) + (3/x)y(x)，则原二阶方程可写成
#
#     z'(x) = (1/x)z(x) + 5 - 6x + 7x^2.
#
# 与两端边界条件相容的 z(x) 为
#
#     z(x) = 6x - 6x^2 + (7/2)x^3 + 5x ln(x).
#
# 因而本题可以转化为下面的一阶初值问题：
#
#     y'(x) = -(3/x)y(x) + 6x - 6x^2 + (7/2)x^3 + 5x ln(x),
#     y(1) = 1/2.
#
# 随后使用 RK4 从 x=1 推进到 x=2，并用原右端边界条件 y(2)=4+4ln(2) 做核对。
#
# 【重要适用限制】
# 上述降阶是利用本题特殊结构和边界条件得到的，不能直接套到任意二阶边值问题。
# 如果原题改变，通常要重新推导一阶初值问题，或改用打靶法、有限差分法等边值问题算法。
#
# 【换题时必须替换的位置】
# 1. exact_solution：替换成新题真解；若题目没有真解，就不能直接计算真误差与真收敛阶。
# 2. ode_rhs：替换成重新推导后的一阶方程右端 f(x,y)，不能只改原二阶方程的系数。
# 3. X_LEFT、X_RIGHT、Y_LEFT、Y_RIGHT_EXPECTED：替换区间、左端初值和右端边界值。
# 4. STEP_DENOMINATORS：替换题目指定的步长；本程序使用 h=1/分母。
# =============================================================================

from pathlib import Path  # Python 语法：从标准库 pathlib 导入 Path，用于稳健地表示文件路径。

import matplotlib.pyplot as plt  # Matplotlib 语法：导入绘图库 pyplot，并使用常用别名 plt。
import numpy as np  # NumPy 语法：导入数值计算库，并使用常用别名 np。


X_LEFT = 1.0  # 【换题时必须替换】数值区间左端点；浮点数 1.0 对应数学中的 x=1。
X_RIGHT = 2.0  # 【换题时必须替换】数值区间右端点；RK4 将从 X_LEFT 推进到这里。
Y_LEFT = 0.5  # 【换题时必须替换】左端初值 y(1)=1/2；0.5 是 1/2 的浮点表示。
Y_RIGHT_EXPECTED = 4.0 + 4.0 * np.log(2.0)  # 【换题时必须替换】原边值问题给出的右端边界值；np.log 是自然对数 ln。
STEP_DENOMINATORS = np.array([200, 100, 50, 25, 10, 5, 4], dtype=int)  # 【换题时必须替换】题目步长 h=1/分母；整数数组可防止误写小数。
SHOW_FIGURES = True  # 【可按需替换】True 表示运行结束时显示图像，False 表示只计算并打印结果。
SAVE_FIGURES = False  # 【可按需替换】True 表示把两幅图保存在本 Python 文件所在目录。


def exact_solution(x):  # 【换题时必须替换】定义真解函数；Python 的 def 关键字用于定义可重复调用的函数。
    """返回题目给出的真解 y(x)，用于计算真实误差。"""  # Python 语法：三引号字符串是函数文档字符串，可由 help() 查看。
    return x**2 - x**3 + 0.5 * x**4 + x**2 * np.log(x)  # 数值分析：真解作为误差基准；Python 的 ** 表示乘方。


def ode_rhs(x, y):  # 【换题时必须替换】定义降阶后的一阶方程 y'=f(x,y) 的右端函数。
    """返回降阶后一阶初值问题的右端 f(x,y)。"""  # 数值分析：单步法只需要在若干阶段点计算右端函数值。
    return -(3.0 / x) * y + 6.0 * x - 6.0 * x**2 + 3.5 * x**3 + 5.0 * x * np.log(x)  # 数值分析：这是由本题特殊结构推导出的 f(x,y)。


def classical_rk4(rhs, x_left, x_right, y_left, h):  # Python 语法：函数参数把方程、区间、初值和步长显式传入，避免依赖全局 h。
    """使用固定步长的经典四级四阶 Runge-Kutta 方法求解标量一阶初值问题。"""  # 数值分析：RK4 的局部截断误差为 O(h^5)，整体误差通常为 O(h^4)。
    interval_length = x_right - x_left  # 数值分析：区间长度除以步长应当得到整数步数，才能恰好走到右端点。
    step_count_float = interval_length / h  # Python 语法：/ 执行浮点除法；这里先得到理论步数的浮点表示。
    step_count = int(round(step_count_float))  # 数值分析：round 消除二进制浮点表示造成的微小偏差，再转成整数步数。
    if not np.isclose(step_count * h, interval_length, rtol=0.0, atol=1e-12):  # Python 语法：if 检查条件；np.isclose 用容差判断两个浮点数是否近似相等。
        raise ValueError("步长 h 必须整除区间长度，才能恰好到达右端点。")  # Python 语法：raise 主动抛出异常，防止使用错误网格继续计算。
    x_grid = np.linspace(x_left, x_right, step_count + 1, dtype=float)  # 数值分析：生成含两个端点的等距网格 x_0,...,x_N。
    y_values = np.zeros(step_count + 1, dtype=float)  # NumPy 语法：预分配全零数组，避免循环中反复扩展数组。
    y_values[0] = y_left  # 数值分析：把初值 y(x_0)=y_left 写入数值解数组的第一个位置。
    for n in range(step_count):  # Python 语法：range(step_count) 依次产生 0,...,N-1，用于完成 N 次推进。
        x_n = x_grid[n]  # Python 语法：方括号执行数组索引；x_n 是当前网格点。
        y_n = y_values[n]  # 数值分析：y_n 是真解 y(x_n) 的当前数值近似。
        k1 = rhs(x_n, y_n)  # RK4 第一阶段：在当前点计算斜率 k1=f(x_n,y_n)。
        k2 = rhs(x_n + 0.5 * h, y_n + 0.5 * h * k1)  # RK4 第二阶段：用 k1 预测半步状态，再计算半步斜率。
        k3 = rhs(x_n + 0.5 * h, y_n + 0.5 * h * k2)  # RK4 第三阶段：用 k2 再次预测半步状态，改善中点斜率。
        k4 = rhs(x_n + h, y_n + h * k3)  # RK4 第四阶段：用 k3 预测整步终点状态，再计算终点斜率。
        weighted_slope = (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0  # 数值分析：按 1:2:2:1 加权平均四个阶段斜率。
        y_values[n + 1] = y_n + h * weighted_slope  # RK4 更新公式：用加权斜率把数值解从 x_n 推进到 x_{n+1}。
    return x_grid, y_values  # Python 语法：return 返回两个数组；调用者可同时得到网格和整条数值解。


def calculate_errors(step_sizes):  # Python 语法：定义批量计算不同步长右端误差的函数。
    """对每个步长运行 RK4，并返回 x=右端点处的绝对误差。"""  # 数值分析：固定终点误差可用于研究方法的整体收敛阶。
    errors = np.zeros_like(step_sizes, dtype=float)  # NumPy 语法：创建与 step_sizes 形状相同的浮点零数组。
    endpoint_values = np.zeros_like(step_sizes, dtype=float)  # 数值分析：保存每个步长对应的右端数值解，便于输出核对。
    exact_endpoint = float(exact_solution(X_RIGHT))  # Python 语法：float 把 NumPy 标量转成普通浮点数；此值是真解 y(2)。
    for index, h in enumerate(step_sizes):  # Python 语法：enumerate 同时给出数组下标 index 和当前元素 h。
        _, y_values = classical_rk4(ode_rhs, X_LEFT, X_RIGHT, Y_LEFT, float(h))  # Python 语法：下划线 _ 表示此处不使用返回的网格数组。
        endpoint_values[index] = y_values[-1]  # Python 语法：索引 -1 表示数组最后一个元素，即右端点数值解。
        errors[index] = abs(endpoint_values[index] - exact_endpoint)  # 数值分析：绝对误差等于数值解与真解之差的绝对值。
    return endpoint_values, errors  # Python 语法：返回各步长的右端近似值和对应绝对误差。


def calculate_pairwise_orders(step_sizes, errors):  # Python 语法：定义利用相邻两组网格估计收敛阶的函数。
    """用相邻步长和误差计算观测收敛阶 p=ln(E2/E1)/ln(h2/h1)。"""  # 数值分析：若 E(h)≈C*h^p，该比值应趋近理论阶 p。
    orders = np.full(step_sizes.size, np.nan, dtype=float)  # NumPy 语法：用 nan 初始化；最细步长没有更细的相邻数据可比较。
    for index in range(1, step_sizes.size):  # Python 语法：从下标 1 开始，使 index-1 始终是合法下标。
        error_ratio = errors[index] / errors[index - 1]  # 数值分析：相邻误差之比近似为相邻步长之比的 p 次方。
        step_ratio = step_sizes[index] / step_sizes[index - 1]  # 数值分析：本题步长按从小到大排列，比例通常大于 1。
        orders[index] = np.log(error_ratio) / np.log(step_ratio)  # 数值分析：两边取自然对数后解出观测收敛阶 p。
    return orders  # Python 语法：返回长度与步长数组相同的阶数数组，第一个元素为 nan。


def print_results(step_sizes, endpoint_values, errors, pairwise_orders, fitted_order, fitted_intercept):  # Python 语法：长参数列表把输出所需数据集中传入。
    """以表格形式打印数值解、误差、相邻收敛阶和拟合收敛阶。"""  # 数值分析：同时看误差大小和阶数，才能判断是否进入渐近误差区间。
    print("\n经典四阶 Runge-Kutta 方法计算结果")  # Python 语法：print 输出字符串；开头的 \n 表示先换行。
    print("-" * 79)  # Python 语法：字符串乘法重复字符，用于生成表格分隔线。
    print(f"{'h':>12} {'RK4 y(2)':>20} {'绝对误差':>20} {'相邻观测阶':>18}")  # Python 语法：f-string 中的 > 表示右对齐，数字表示字段宽度。
    print("-" * 79)  # Python 语法：再次打印分隔线，使结果更易阅读。
    for h, endpoint, error, order in zip(step_sizes, endpoint_values, errors, pairwise_orders):  # Python 语法：zip 把四个数组的同位置元素组合起来遍历。
        order_text = "---" if np.isnan(order) else f"{order:.8f}"  # Python 语法：条件表达式根据是否为 nan 决定显示内容。
        print(f"{h:12.8f} {endpoint:20.12f} {error:20.12e} {order_text:>18}")  # Python 语法：.12e 使用科学计数法输出误差。
    print("-" * 79)  # Python 语法：打印结果表格底部分隔线。
    print(f"真解给出的 y(2)             = {exact_solution(X_RIGHT):.12f}")  # 数值分析：真解端点值是误差计算基准。
    print(f"原右端边界条件给出的 y(2) = {Y_RIGHT_EXPECTED:.12f}")  # 数值分析：该值应与真解端点值一致，否则题目数据或真解可能有误。
    print(f"log(error) 对 log(h) 的线性拟合：log(error) = {fitted_order:.12f} log(h) + {fitted_intercept:.12f}")  # 数值分析：拟合直线斜率就是总体观测收敛阶。
    print(f"拟合得到的误差收敛阶 p = {fitted_order:.12f}")  # 数值分析：经典 RK4 理论整体收敛阶为 4，拟合值应接近 4。


def plot_results(step_sizes, errors, fitted_order, fitted_intercept):  # Python 语法：定义绘制误差图和对数拟合图的函数。
    """绘制误差随步长变化的图，以及用于估计收敛阶的双对数拟合图。"""  # 数值分析：双对数坐标中 E≈C*h^p 表现为斜率 p 的直线。
    figure_error, axis_error = plt.subplots(figsize=(8.0, 5.0))  # Matplotlib 语法：创建一幅图 figure 和对应坐标轴 axis。
    axis_error.loglog(step_sizes, errors, "o-", label="RK4 endpoint error")  # 数值分析：loglog 将横纵轴都设为对数尺度，便于观察幂律关系。
    axis_error.set_xlabel("step size h")  # Matplotlib 语法：设置误差图横轴标签。
    axis_error.set_ylabel("absolute error at x=2")  # Matplotlib 语法：设置误差图纵轴标签。
    axis_error.set_title("RK4 endpoint error versus step size")  # Matplotlib 语法：设置误差图标题。
    axis_error.grid(True, which="both", linestyle="--", alpha=0.5)  # Matplotlib 语法：显示主次网格；alpha 控制透明度。
    axis_error.legend()  # Matplotlib 语法：根据 label 参数显示图例。
    figure_error.tight_layout()  # Matplotlib 语法：自动调整边距，减少标签被裁切的风险。
    log_h = np.log(step_sizes)  # 数值分析：对步长取自然对数，把幂律误差模型转化为线性模型。
    log_error = np.log(errors)  # 数值分析：对误差取自然对数；要求误差严格大于 0。
    fitted_log_error = fitted_order * log_h + fitted_intercept  # 数值分析：根据拟合斜率和截距计算拟合直线上的纵坐标。
    figure_fit, axis_fit = plt.subplots(figsize=(8.0, 5.0))  # Matplotlib 语法：创建第二幅图及其坐标轴。
    axis_fit.plot(log_h, log_error, "o-", label="computed log(error)")  # Matplotlib 语法：绘制实际对数误差数据并连接相邻点。
    axis_fit.plot(log_h, fitted_log_error, "--", label=f"linear fit, slope={fitted_order:.6f}")  # 数值分析：拟合线斜率即总体观测阶。
    axis_fit.set_xlabel("log(h)")  # Matplotlib 语法：设置拟合图横轴标签。
    axis_fit.set_ylabel("log(error)")  # Matplotlib 语法：设置拟合图纵轴标签。
    axis_fit.set_title("Linear fit for convergence-order estimation")  # Matplotlib 语法：设置拟合图标题。
    axis_fit.grid(True, linestyle="--", alpha=0.5)  # Matplotlib 语法：显示虚线网格以便观察线性程度。
    axis_fit.legend()  # Matplotlib 语法：显示实际数据与拟合直线的图例。
    figure_fit.tight_layout()  # Matplotlib 语法：自动优化第二幅图的边距。
    if SAVE_FIGURES:  # Python 语法：只有 SAVE_FIGURES 为 True 时才执行缩进块中的保存操作。
        output_directory = Path(__file__).resolve().parent  # Python 语法：__file__ 是当前脚本路径；resolve().parent 得到脚本所在目录。
        figure_error.savefig(output_directory / "211870125_RK4_endpoint_error.png", dpi=180)  # Matplotlib 语法：以 180 dpi 保存误差图。
        figure_fit.savefig(output_directory / "211870125_RK4_convergence_fit.png", dpi=180)  # Matplotlib 语法：以 180 dpi 保存对数拟合图。
    if SHOW_FIGURES:  # Python 语法：只有 SHOW_FIGURES 为 True 时才显示图形窗口。
        plt.show()  # Matplotlib 语法：显示此前创建的所有图；在交互式后端中会等待用户关闭窗口。
    else:  # Python 语法：当 SHOW_FIGURES 为 False 时执行此分支。
        plt.close("all")  # Matplotlib 语法：关闭全部图形，释放内存并避免后台运行时残留窗口。


def main():  # Python 语法：把程序主流程封装为 main 函数，便于测试和复用。
    """执行题目要求的全部计算、收敛阶估计、结果输出与绘图。"""  # 数值分析：主流程对应“多步长计算→误差→阶数→分析”。
    step_sizes = 1.0 / STEP_DENOMINATORS.astype(float)  # NumPy 语法：astype(float) 转为浮点数组，再逐元素计算 h=1/分母。
    if not np.isclose(exact_solution(X_LEFT), Y_LEFT, rtol=0.0, atol=1e-12):  # 数值分析：检查给定真解是否满足左端初值，防止数据抄写错误。
        raise ValueError("题目给出的真解不满足左端初值，请检查 exact_solution 或 Y_LEFT。")  # Python 语法：发现输入不一致时立即抛出异常。
    if not np.isclose(exact_solution(X_RIGHT), Y_RIGHT_EXPECTED, rtol=0.0, atol=1e-12):  # 数值分析：检查真解是否满足原问题右端边界条件。
        raise ValueError("题目给出的真解不满足右端边界条件，请检查 exact_solution 或 Y_RIGHT_EXPECTED。")  # Python 语法：防止在错误基准上分析误差。
    endpoint_values, errors = calculate_errors(step_sizes)  # Python 语法：多重赋值接收函数返回的两个数组。
    if np.any(errors <= 0.0):  # NumPy 语法：np.any 判断数组中是否至少有一个元素满足条件。
        raise ValueError("存在非正误差，无法对误差取对数并估计收敛阶。")  # 数值分析：log(error) 只对严格正误差有定义。
    pairwise_orders = calculate_pairwise_orders(step_sizes, errors)  # 数值分析：由相邻步长误差比计算局部观测收敛阶。
    fitted_order, fitted_intercept = np.polyfit(np.log(step_sizes), np.log(errors), deg=1)  # NumPy 语法：polyfit 做一次多项式最小二乘拟合；deg=1 表示直线。
    print_results(step_sizes, endpoint_values, errors, pairwise_orders, fitted_order, fitted_intercept)  # Python 语法：调用函数打印完整数值结果。
    plot_results(step_sizes, errors, fitted_order, fitted_intercept)  # Python 语法：调用函数绘制误差图和收敛阶拟合图。


if __name__ == "__main__":  # Python 语法：仅当该文件被直接运行时条件成立；被其他文件导入时不会自动执行。
    main()  # Python 语法：调用主函数，开始执行本题数值计算。
