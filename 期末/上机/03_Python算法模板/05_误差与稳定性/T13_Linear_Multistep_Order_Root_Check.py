# 一般线性多步法：阶数条件、特征根条件与绝对稳定检查。  # 对应教材基于 Taylor 展开的待定系数法和根条件。
import math  # 标准库 math 用于阶乘。
import numpy as np  # NumPy 用于系数向量、多项式根和向量化计算。

alpha = np.array([0,-1,1], dtype=float)  # 【必须替换】按 [alpha_0,...,alpha_k] 填写；默认是二步 Adams-Bashforth 左端。
beta = np.array([-1/2,3/2,0], dtype=float)  # 【必须替换】按 [beta_0,...,beta_k] 填写；显式法最后一项通常为 0。
max_order_to_check = 8  # 【必须替换】最多检查到多少阶。
assert alpha.size == beta.size, "alpha 与 beta 长度必须相同。"  # 同一个线性多步公式必须具有同样多的系数。
k = alpha.size-1  # k 是多步法步数；数学下标从 0 到 k。
j = np.arange(k+1, dtype=float)  # 创建节点编号数组 [0,1,...,k]。
c = np.zeros(max_order_to_check+1, dtype=float)  # c[q] 对应教材中的阶条件系数 c_q。
c[0] = np.sum(alpha)  # c_0=sum(alpha_j)；相容方法必须满足 c_0=0。
for q in range(1, max_order_to_check+1):  # 依次计算 c_1,c_2,...。
    c[q] = np.sum(alpha*j**q)/math.factorial(q)-np.sum(beta*j**(q-1))/math.factorial(q-1)  # Taylor 展开后的通用阶条件。
nonzero = np.flatnonzero(np.abs(c)>1e-12)  # flatnonzero 返回所有不满足零条件的位置。
if nonzero.size == 0:  # 检查范围内所有 c_q 都近似为零。
    print(f"至少为 {max_order_to_check} 阶；请增大 max_order_to_check 继续检查。")  # 当前只能给出阶数下界。
else:  # 找到了第一个非零阶条件系数。
    first = int(nonzero[0])  # 取第一个非零系数下标。
    order = first-1  # c_0,...,c_p 为零而 c_{p+1} 非零时，方法为 p 阶。
    print(f"判定阶数 p={order}，首个非零系数 c_{first}={c[first]:.6g}。")  # 输出阶数判断。
rho_roots = np.roots(alpha[::-1])  # NumPy roots 要求最高次系数在前，因此用 [::-1] 反转系数。
inside_unit = np.all(np.abs(rho_roots) <= 1+1e-10)  # 根条件第一部分：所有根都在闭单位圆内。
unit_roots = rho_roots[np.abs(np.abs(rho_roots)-1)<1e-8]  # 筛出位于单位圆上的根。
unit_simple = all(np.sum(np.abs(rho_roots-root)<1e-8)==1 for root in unit_roots)  # 根条件第二部分：单位圆上的根必须单重。
root_condition = inside_unit and unit_simple  # 同时满足两部分才满足特征根条件。
for root in rho_roots:  # 遍历输出特征多项式的所有根。
    print(f"rho root={root}, modulus={abs(root):.10f}")  # 打印每个根及其模。
print(f"特征根条件是否满足：{root_condition}")  # 相容加根条件可推出收敛。
z = -1.0  # 【必须替换】测试 z=lambda*h；衰减问题通常 z<0。
stability_roots = np.roots((alpha-z*beta)[::-1])  # 绝对稳定特征方程 rho(xi)-z*sigma(xi)=0。
print(f"z={z:.6g} 时最大根模={np.max(np.abs(stability_roots)):.10f}；严格小于 1 表示绝对稳定。")  # 输出固定 z 的稳定性判断依据。
