# 上机考试模板索引

> 每个模板顶部都有详细注释，找到对应题型后，只需修改"① 修改这里"标注的参数区域即可直接使用。

---

## 平方逼近（对应 3-上机）

| 文件 | 题型 | 核心公式 | 关键函数 |
|------|------|----------|---------|
| [T01_最佳平方逼近_幂函数基.m](T01_最佳平方逼近_幂函数基.m) | 幂函数基 $\varphi_k=x^k$ 的最佳平方逼近 | `G = hilb(n+1)` → `c = G \ f_vec` | `hilb`, `integral`, `cond` |
| [T02_最佳平方逼近_Legendre正交基.m](T02_最佳平方逼近_Legendre正交基.m) | Legendre 正交基的最佳平方逼近 | $c_k=(2k+1)\int_0^1 f\cdot P_k(2x-1)dx$，无需解方程组 | `integral`，含 `legendre_on_01` 局部函数 |
| [T03_最小二乘拟合_多项式基.m](T03_最小二乘拟合_多项式基.m) | 离散数据多项式最小二乘拟合 | Vandermonde 设计矩阵 `A` → `c = A \ y` | `cond`, `mean`, `repmat` |
| [T04_最小二乘拟合_混合基函数.m](T04_最小二乘拟合_混合基函数.m) | 混合基 $\{1, x, \cos x, \sin x\}$ 等 | `A_mix = [ones, x, cos(x), sin(x)]` → `A \ y` | `strjoin`, `stem` |
| [T12_条件数分析与Chebyshev基.m](T12_条件数分析与Chebyshev基.m) | 多项式拟合病态性分析 + Chebyshev 基改善 | $T_k(t)=\cos(k\arccos t)$，递推：$T_{k+1}=2tT_k-T_{k-1}$ | `cond`，含 `chebyshev_T` 局部函数 |
| [T13_FFT信号压缩.m](T13_FFT信号压缩.m) | FFT 频域压缩（信号 + 图像） | `X = fft(x)` → 截断 → `x = real(ifft(X))` | `fft`, `ifft`, `fft2`, `ifft2`, `fftshift` |

---

## 数值积分（对应 4-上机）

| 文件 | 题型 | 核心公式 | 关键代码片段 |
|------|------|----------|------------|
| [T05_复合梯形公式.m](T05_复合梯形公式.m) | 复合梯形，误差随步长分析 | $T_n=\frac{h}{2}[f(a)+f(b)+2\sum f(x_i)]$，$O(h^2)$ | `h/2*(fx(1)+fx(end)+2*sum(fx(2:end-1)))` |
| [T06_复合Simpson公式.m](T06_复合Simpson公式.m) | 复合 Simpson，与梯形对比 | $S_m=\frac{h}{3}[\cdots+4\sum_\text{奇}+2\sum_\text{偶}]$，$O(h^4)$ | `h/3*(fx(1)+fx(end)+4*sum(fx(2:2:end-1))+2*sum(fx(3:2:end-2)))` |
| [T07_Romberg积分法.m](T07_Romberg积分法.m) | Romberg 外推，高精度积分 | $T_{m,j}=\frac{4^{j-1}T_{m,j-1}-T_{m-1,j-1}}{4^{j-1}-1}$，$O(h^{2m})$ | 变步长递推 + 二重循环外推 |
| [T08_自适应Simpson积分法.m](T08_自适应Simpson积分法.m) | 自适应 Simpson，递归细分 | 误差估计 $\approx\|S_2-S_{ab}\|/15$，满足则接受 | 含 `simp1` + `adapt_simp` 递归局部函数 |
| [T09_广义积分_无穷上限处理.m](T09_广义积分_无穷上限处理.m) | 无穷积分：截断 + 变量替换 | 截断法 + 变换 $t=e^{-x}$ 映射到 $[0,1]$ | `arrayfun`（处理非向量化函数），`integral(f,0,Inf)` |
| [T10_蒲丰投针_MonteCarlo估计pi.m](T10_蒲丰投针_MonteCarlo估计pi.m) | 蒲丰投针估计 $\pi$ | $P=2l/({\pi d})$，误差 $O(1/\sqrt{N})$ | `rand`, `sum`, `cumsum`, `normpdf` |
| [T11_MonteCarlo积分_两种方法.m](T11_MonteCarlo积分_两种方法.m) | 随机投点法 + 样本均值法 | 均值法：$I\approx(b-a)\cdot\frac{1}{N}\sum f(x_i)$ | `rand`, `mean`, `std`, `histogram` |

---

## 迁移使用步骤

每个模板的结构都一样，按以下步骤修改：

```
① 找到文件顶部的 "===== ① 修改这里 =====" 区域
② 替换被积函数 f、区间 [a,b]、精确值 I_exact（或数据 x_data/y_data）
③ 调整参数（n 次数、N 点数、m 子区间数等）
④ 直接运行，结果和图形自动生成
```

---

## 三条最重要的记忆点

| 易错点 | 正确做法 |
|--------|---------|
| 被积函数定义 | 必须用逐元素运算：`@(x) 4 ./ (1+x.^2)` 而非 `4/(1+x^2)` |
| Simpson 节点下标 | 奇数节点（系数4）：`fx(2:2:end-1)`；偶数内点（系数2）：`fx(3:2:end-2)` |
| 幂函数基 vs Legendre 基 | 幂函数基 → `hilb(n+1)` 极度病态；Legendre 基 → Gram 矩阵对角，条件数=1 |
