import numpy as np
from scipy.special import legendre

def compute_coefficients(f, d, max_degree):
    # 生成均匀采样点
    num_points = 1000  # 调整此数以平衡精度和效率
    x = np.linspace(-d, d, num_points)
    t = x / d
    
    M = max_degree + 1
    weights = [0.0] * M
    
    for n in range(M):
        weight_numerator = d
        weight_denominator = 2 * n + 1
        if weight_denominator == 0:
            # 当n= -0.5，但n从0开始，所以这行可能不会执行
            weights[n] = float('inf')
        else:
            weights[n] = weight_numerator / weight_denominator
    
    A = np.zeros((num_points, M))
    
    for n in range(M):
        if 2 * n + 1 == 0:
            # 处理n=-0.5的情况，但这里n从0开始，所以这可能不会触发
            continue
        p = legendre(n)
        A[:, n] = p(t)
    
    dx = (2 * d) / num_points
    
    # 计算观测向量b
    b = np.zeros(num_points)
    for i in range(num_points):
        f_val = f(x[i])
        b[i] = f_val * dx
    
    a = A.T @ b  # 矩阵乘法：A^T @ b 是一个M维向量
    c = a / weights  # 每个系数除以对应的权重
    
    return c