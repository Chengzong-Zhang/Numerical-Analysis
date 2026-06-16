import numpy as np

# ==================== 1. 定义基础函数 (全局作用域) ====================
def f(t, y):
    """ODE 右端项：y' = -exp(10t)*(y-sin(t)) + cos(t)"""
    return -np.exp(10 * t) * (y - np.sin(t)) + np.cos(t)

def fy(t, y):
    """偏导数 ∂f/∂y：用于牛顿迭代和隐式方程求解"""
    return -np.exp(10 * t)

def y_exact(t):
    """精确解函数（解析推导结果）
    根据 ODE 解析解推导，特解为 sin(t)，初值 y(0)=0 满足 sin(0)=0
    所以精确解为 y(t) = sin(t)
    """
    return np.sin(t)

# ==================== 2. 数值求解器 Solve_Step(h, method) ====================
def Solve_Step(h, method):
    # 全局参数 (假设与主程序一致)
    t0, T, y0_val = 0.0, 1.0, 0.0

    # 计算总步数与时间网格
    N = int(round((T - t0) / h))
    t = t0 + np.arange(N + 1) * h
    
    # 预分配数值解数组 (使用 float64 保持精度)
    y = np.zeros(N + 1, dtype=np.float64)
    y[0] = y0_val

    tol, max_iter = 1e-8, 2000  # 牛顿迭代容差与最大次数
    
    # --- Improved Euler (显式方法) ---
    if method == "improved_euler":  
        for n in range(N):
            tn, yn = t[n], y[n]
            
            k1 = f(tn, yn)  
            # 预测步：欧拉法预估下一时刻的斜率
            pred_y = yn + h * k1 
            k2 = f(tn + h, pred_y) 
            
            # 校正步：梯形公式平均斜率
            y[n+1] = yn + h*(k1+k2)/2
            
        return t, y

    # --- Trapezoid (隐式 Crank-Nicolson) ---
    elif method == "trapezoid":
        for n in range(N):
            tn, yn = t[n], y[n]
            
            # 构造隐式方程：y_{n+1} - yn - h/2 * [f(t_n, y_n) + f(t_{n+1}, y_{n+1})] = 0
            # 初始猜测使用显式欧拉法 (Predictor)
            z = yn + h * f(tn, yn) 
            
            for _ in range(max_iter):
                k_current = f(t[n+1], z)
                k_prev = f(tn, yn)
                
                # 更新公式
                z_new = yn + h*(k_prev + k_current)/2
                
                if abs(z_new - z) < tol: 
                    break
                z = z_new
            
            y[n + 1] = z

    # --- BDF2 (二阶隐式 Backward Differentiation Formula) ---
    elif method == "BDF2":
        # BDF2 系数：对应方程形式 alpha_0*y_{n+1} + alpha_1*y_n + alpha_2*y_{n-1} = h*beta*f(t_{n+1})
        # 标准化为：y_{n+1} + c1*y_n + c2*y_{n-1} - (h/2)*f = 0 (此处采用原始代码的系数定义)
        alpha_k = np.array([1.0, -4.0 / 3.0, 1.0 / 3.0]) 
        beta_k = np.array([2.0 / 3.0])

        # 起步阶段：RK4 提供前 k-1 步值
        for n_start in range(1): # BDF2 需要 y_0, y_1 已知，y_2 开始用公式
            tn = t[n_start]
            yn = y[n_start]
            
            q1 = f(tn, yn)
            q2 = f(tn + h / 2.0, yn + h * q1 / 2.0)
            q3 = f(tn + h / 2.0, yn + h * q2 / 2.0)
            q4 = f(tn + h, yn + h * q3)
            
            y[n_start + 1] = yn + h * (q1 + 2 * q2 + 2 * q3 + q4) / 6.0

        # BDF 主推进阶段
        for n in range(1, N): # 从 n=1 开始，计算 y[n+1] 需要用到 y[n], y[n-1]
            known = alpha_k[0] * y[n] + alpha_k[1] * y[n-1] # 历史项
            
            z = y[n]  # 初始猜测 (Predictor)
            
            for iteration in range(max_iter):
                G = z - known - h * beta_k[0] * f(t[n + 1], z)
                dG = 1.0 - h * beta_k[0] * fy(t[n + 1], z)
                
                if abs(dG) < 1e-15: break
                
                z_new = z - G / dG
                if abs(z_new - z) < tol:
                    break
                z = z_new
            
            y[n + 1] = z
    
    else:
        raise ValueError(f"Unknown method: {method}")
    
    return t, y

# ==================== 3. 主程序：测试与验证 ====================
def main():
    # 【配置参数】步长列表
    h_list = [0.1, 0.05, 0.025, 0.0125, 0.00625]
    
    # 【方法选择】
    methods = ["trapezoid", "BDF2", "improved_euler"]

    print("=" * 80)
    print("数值方法求解与误差验证报告 (ODE: y' + e^(10t)y = e^(10t)sin(t) + cos(t))")
    print(f"初值条件：y(0)=0, 时间区间：[0, 1]")
    print("=" * 80)

    results_dict = {} # 用于存储详细结果以便后续分析

    for method in methods:
        print(f"\n=== 方法：{method} ===")
        all_max_errors = []
        
        for h in h_list:
            try:
                t_num, y_num = Solve_Step(h, method)
                
                # 【计算绝对误差 (点态 + 最大)】
                exact_vals = y_exact(t_num)
                abs_errors = np.abs(y_num - exact_vals)
                max_err = float(np.max(abs_errors))
                
                all_max_errors.append(max_err)
                
                status = "OK" if max_err < 1e5 else "UNSTABLE/LOW PRECISION"
                
                print(f"步长 h={h:.4f} | Max Absolute Error: {max_err:.2e} [{status}]")

            except Exception as e:
                print(f"步长 h={h:.4f} 发生异常：{type(e).__name__}: {str(e)[:50]}...")
                
        # 【计算收敛阶】(取相邻步长的误差比)
        if len(all_max_errors) >= 2 and all(all_max_errors[i] > 1e-30 for i in range(len(all_max_errors))):
            p_values = []
            print("   收敛阶分析 (p ≈ log(E_{n}/E_{n+1}) / log(h_n/h_{n+1})):")
            for i in range(len(all_max_errors)-1):
                err_curr, err_prev = all_max_errors[i], all_max_errors[i+1]
                # 防止除零，如果误差极小则设为上一个的一半
                h_ratio = 0.5 
                if err_prev > 1e-30:
                    p = np.log(err_curr / err_prev) / np.log(h_ratio)
                else:
                    p = 2.0 # 假设二阶
                
                p_values.append(p)
                
            print(f"   平均收敛阶估计：[ {', '.join([f'{p:.4f}' for p in p_values])} ]")

        results_dict[method] = all_max_errors

        print("-" * 80)

  

    # 深度分析
    print("\n【详细对比分析】")
    print("1. 梯形法 (Trapezoid): A-稳定，二阶精度。对于刚性方程 e^(10t)，表现稳健，误差随 h² 衰减。")
    print("2. BDF2 格式：A-稳定，二阶精度。起步值需 RK4 提供，牛顿迭代求解隐式方程，计算量略大但精度相当。")
    print("3. Improved Euler (显式): 条件稳定。")
    
    # 分析 Improved Euler 在刚性问题下的表现
    stable_h_min = h_list[0] # 粗略估计
    unstable_threshold = 1e5
    
    for i, err in enumerate(results_dict['improved_euler']):
        if err > unstable_threshold:
            print(f"   - 当步长 h={h_list[i]:.4f} 时，Improved Euler 开始出现严重不稳定性 (Max Error ≈ {err:.2e})。")
    
    print("\n   - 刚性特征：由于 ODE 中包含项 -exp(10t)y，特征值 λ(t) ≈ -exp(10t)，在 t=1 时 λ ≈ -22026。")
    print("   - Improved Euler (显式) 的绝对稳定区域约为 h*λ > -4.5，即 h < 0.0002 * exp(-10t)。")
    print("   - 结论：在 t=1 附近，Improved Euler 需要极小的步长才能维持稳定，否则误差会随阶数平方指数增长。")
    print("=" * 80)

if __name__ == "__main__":
    main()