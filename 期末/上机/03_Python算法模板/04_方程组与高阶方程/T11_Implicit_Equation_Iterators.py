# 隐式 ODE 方法所需的不动点迭代与 Newton-Raphson 迭代。  # 数值分析：用于求解每一步产生的非线性方程。
import math  # 标准库 math 提供 cos 与 sin。

G = lambda z: z-math.cos(z)  # 【必须替换】把隐式一步公式整理成 G(z)=0。
dG = lambda z: 1+math.sin(z)  # 【必须替换】G 对 z 的导数；向量问题应换成 Jacobian 矩阵。
g = lambda z: math.cos(z)  # 【必须替换】若能整理成 z=g(z)，供不动点迭代使用。
z0, tol, max_iter = 0.5, 1e-12, 100  # 【必须替换】初始猜测、容差、迭代上限。
z_fixed = z0  # 初始化不动点迭代当前值。
for _ in range(max_iter):  # 不动点迭代局部满足 |g'(z*)|<1 时通常收敛。
    z_new = g(z_fixed)  # 计算 z^{k+1}=g(z^k)。
    if abs(z_new-z_fixed) < tol: break  # 相邻迭代差作为实用停止准则。
    z_fixed = z_new  # 更新当前值并继续迭代。
fixed_result = z_new  # 保存最终不动点迭代结果。
z_newton = z0  # 初始化 Newton 迭代当前值。
for _ in range(max_iter):  # Newton 在单根附近通常二次收敛。
    delta = G(z_newton)/dG(z_newton)  # 标量 Newton 修正量；向量问题使用线性方程求解，不要显式求逆。
    z_newton_new = z_newton-delta  # Newton 更新 z^{k+1}=z^k-G/G'。
    if abs(z_newton_new-z_newton) < tol: break  # 修正量足够小时停止。
    z_newton = z_newton_new  # 更新当前值并继续迭代。
print(f"fixed-point = {fixed_result:.15g}, residual = {abs(G(fixed_result)):.3e}")  # residual 检查是否真正满足方程。
print(f"Newton      = {z_newton_new:.15g}, residual = {abs(G(z_newton_new)):.3e}")  # 输出 Newton 结果与残差。
