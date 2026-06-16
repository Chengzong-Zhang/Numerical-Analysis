import numpy as np

# ==================== 1. 定义基础函数 (全局作用域) ====================
def f(t, y):
    """ODE 右端项：y' = -exp(10t)*(y-sin(t)) + cos(t)"""
    return -np.exp(10 * t) * (y - np.sin(t)) + np.cos(t)

def fy(t, y):
    """偏导数 ∂f/∂y：用于牛顿迭代和隐式方程求解"""
    return -np.exp(10 * t)

def y_exact(t):
    """精确解函数（解析推导结果）"""
    # 根据 ODE 解析解推导，特解为 sin(t)，初值 y(0)=0 满足 sin(0)=0
    return np.sin(t)

# ==================== 2. 数值求解器 Solve_Step(h, method) ====================
def Solve_Step(h, method):
    # 全局参数 (假设与主程序一致)
    t0 = 0.0
    T = 1.0
    y0_val = 0.0

    # 计算总步数与时间网格
    N = int(round((T - t0) / h))
    t = t0 + np.arange(N + 1) * h

    # 预分配数值解数组
    y = np.zeros(N + 1, dtype=float)
    y[0] = y0_val

    # 迭代参数 (常量定义)
    tol, max_iter = 1e-6, 2000

    if method == "trapezoid":
        for n in range(N):
            tn, yn = t[n], y[n]
            z = yn + h * f(tn, yn)
            # Newton 迭代求解隐式方程：z - yn - (h/2)[f(tn,yn) + f(tn+h,z)] = 0
            for _ in range(max_iter):
                z_new = yn + h*(f(tn, yn)+f(t[n+1], z))/2  # 梯形公式 z=y_n+h/2*(f_n+f_{n+1})。
                if abs(z_new-z) < tol: break  # 达到迭代容差后停止。
                z = z_new  # 更新隐式未知量。
            y[n + 1] = z


    elif method == "BDF2":
        k = 2
        # BDF2 系数：标准形式 y_{n+1} - 4/3*y_n + 1/3*y_{n-1} = h * 2/3 * f(t_{n+1})
        # 对应系数 alpha (y_{n+1}, y_n, y_{n-1}): [1, -4/3, 1/3]
        # beta: 2/3
        coeffs_alpha = np.array([1.0, -4.0 / 3.0, 1.0 / 3.0])
        beta_val = 2.0 / 3.0

        # 起步阶段 (RK4 提供前 k-1=1 步值，即 y[0], y[1])
        for n_start in range(k - 1):
            tn = t[n_start]
            yn = y[n_start]
            
            q1 = f(tn, yn)
            q2 = f(tn + h / 2.0, yn + h * q1 / 2.0)
            q3 = f(tn + h / 2.0, yn + h * q2 / 2.0)
            q4 = f(tn + h, yn + h * q3)
            
            y[n_start + 1] = yn + h * (q1 + 2*q2 + 2*q3 + q4) / 6.0

        # BDF 主推进阶段
        for n in range(k - 1, N):
            # 计算已知项 (历史值移项到右侧，注意系数对应 y[n], y[n-1])
            known = coeffs_alpha[1] * y[n] + coeffs_alpha[2] * y[n - 1]
            
            z = y[n]  # 初始猜测
            
            # Newton 迭代求解隐式方程：z - known - h*beta*f(t_{n+1}, z) = 0
            for iteration in range(max_iter):
                G = z - known - h * beta_val * f(t[n + 1], z)
                
                dG = 1.0 - h * beta_val * fy(t[n + 1], z)
                
                if abs(dG) < 1e-15:
                    break
                
                delta = G / dG
                z_new = z - delta
                
                if abs(z_new - z) < tol:
                    z = z_new
                    break
                z = z_new
            
            y[n + 1] = z
    
    elif method == "improved_euler":  
        for n in range(N):
            tn, yn = t[n], y[n]
            k1 = f(tn, yn)  
            k2 = f(tn+h, yn+h*k1) 
            y[n+1] = yn + h*(k1+k2)/2 



    return t, y  # 返回时间和解，便于后续误差对齐

# ==================== 3. 主程序：测试与验证 ====================
def main():
    # 【配置参数】步长列表
    h_list = [0.1, 0.05, 0.025, 0.0125, 0.00625]
    
    # 【方法选择】
    methods = ["trapezoid", "BDF2","improved_euler"]

    print("=" * 80)
    print("数值方法求解与误差验证报告")
    print("ODE: y' + e^(10t)y = e^(10t)sin(t) + cos(t), y(0)=0, t in (0, 1]")
    print("=" * 80)

    for method in methods:
        print(f"\n=== 方法：{method} ===")
        all_max_errors = []
        
        # 遍历步长列表
        for h in h_list:
            try:
                # 【调用主算法得到结果】
                t_num, y_num = Solve_Step(h, method)
                
                # 【计算绝对误差 (点态 + 最大)】
                # 精确解直接传入时间数组
                exact_vals = y_exact(t_num)
                abs_errors = np.abs(y_num - exact_vals)
                max_err = np.max(abs_errors)
                
                all_max_errors.append(max_err)
                
                print(f"步长 h={h:.5f} | Max Absolute Error: {max_err:.2e}")

            except Exception as e:
                print(f"步长 h={h:.5f} 发生异常：{e}")
        
        # 【计算收敛阶】(取相邻步长的误差比)
        if len(all_max_errors) >= 2:
            p_values = []
            for i in range(len(all_max_errors)-1):
                err_curr, err_prev = all_max_errors[i], all_max_errors[i+1]
                # h_current / h_previous = 0.5
                p = np.log(err_curr / err_prev) / np.log(2.0)
                p_values.append(p)
            
            print(f"收敛阶估计: [ {', '.join([f'{p:.4f}' for p in p_values])} ]")

        print("-" * 80)

    # 【最终比较分析】
    print("\n=== 结果比较与分析 ===")
    print("梯形法 (Trapezoid): A-稳定，二阶精度。对于刚性方程 e^(10t)，表现稳健，误差随 h² 衰减。")
    print("BDF2 格式：A-稳定，二阶精度。起步值需 RK4 提供，牛顿迭代求解隐式方程，计算量略大但精度相当。")
    print("观察：若梯形法收敛阶接近 2.0，BDF2 也接近 2.0，则说明两者在该问题下均发挥预期性能。")

if __name__ == "__main__":
    main()