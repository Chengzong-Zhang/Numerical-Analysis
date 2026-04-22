%% =========================================================
%  problem3_poly_fitting_condition.m
%  功能：利用例3.5.3数据求9次多项式拟合，分析条件数与均方误差
%
%  题目三个子问题：
%  (1) 计算9次多项式拟合的均方误差和法方程组系数矩阵 A 的条件数
%  (2) 解释：若 f 是 Lagrange 插值多项式，理论 RMS=0 但数值 RMS≠0 的原因
%  (3) 将基函数换为 Chebyshev 多项式 {T_k(x)}，观察条件数和均方误差的变化
%
%  【数值分析背景】
%  数据点 m=10，多项式次数 n=9：参数个数=数据点数=10
%  设计矩阵 A（10×10）是方阵 → 法方程组 A^T A c = A^T y 变为方阵方程组
%  理论上9次多项式可精确过所有10个数据点（等价于9次Lagrange插值多项式）
%  → 理论 RMS = 0，但幂函数基的 Vandermonde 矩阵高度病态，数值 RMS ≠ 0
%
%  【Chebyshev 多项式（教材§3.3(一)）】
%  T_k(t) = cos(k·arccos(t))，t∈[-1,1]
%  递推：T_0=1, T_1=t, T_{k+1}(t) = 2t·T_k(t) - T_{k-1}(t)
%  关键性质：在所有首一 n 次多项式中，T_n(t)/2^{n-1} 的无穷范数最小（教材§3.3(一)）
%  作为最小二乘基函数：设计矩阵条件数远小于 Vandermonde 矩阵
%  =========================================================

clear; close all; clc;

fprintf('==============================================\n');
fprintf('  第三章上机作业3：9次多项式拟合与条件数分析\n');
fprintf('==============================================\n\n');

%% =========================================================
%  数据（例3.5.3）
%  =========================================================

x_data = [0.000, 0.895, 1.641, 2.512, 3.542, 4.054, 4.602, 5.063, 5.354, 5.617];
y_data = [1.000, 1.803, 3.680, 7.320, 13.59, 17.41, 22.19, 24.89, 26.55, 29.77];

m = length(x_data);   % m = 10 个数据点
n = m - 1;            % n = 9 次多项式（参数个数 = 数据点数）

fprintf('数据点数 m = %d，多项式次数 n = %d\n', m, n);
fprintf('参数个数（n+1）= %d = 数据点数 m → 理论上等价于 Lagrange 插值\n\n', n+1);

% 绘图用的密集点
x_plot = linspace(min(x_data), max(x_data), 500);

%% =========================================================
%  第(1)部分：幂函数基 {1, x, x^2, ..., x^9}
%  计算均方误差与法方程系数矩阵 A 的条件数
%  =========================================================

fprintf('===== 第(1)部分：幂函数基下的9次多项式拟合 =====\n\n');

% ---- 构造 Vandermonde 设计矩阵 A（10×10 方阵）----
% A(i, k+1) = x_i^k，i=1,...,10; k=0,...,9
% 数值分析含义：当 m=n+1 时，A 是方阵（Vandermonde 矩阵）
% Vandermonde 矩阵：各行是数据点 x_i 的幂次向量 [1, x_i, x_i^2, ..., x_i^n]
% 当 x_i 各不相同时 Vandermonde 矩阵可逆，但高度病态
A_mono = zeros(m, n+1);
for k = 0:n
    % x_data(:).^k：列向量化后逐元素取 k 次方
    A_mono(:, k+1) = x_data(:).^k;
end

% ---- 条件数分析 ----
% cond(A)：矩阵 2-范数条件数（最大奇异值/最小奇异值）
% 数值分析含义：
%   设计矩阵 A 的条件数 = κ(A)
%   法方程矩阵 A^T A 的条件数 = κ(A^T A) = κ(A)^2（条件数平方！）
%   解的相对误差 ≤ κ(A) × (机器精度 eps ≈ 2.2e-16)
kappa_A = cond(A_mono);
kappa_ATA = cond(A_mono' * A_mono);

fprintf('Vandermonde 设计矩阵 A（%d×%d）：\n', m, n+1);
fprintf('  cond(A)         = %.6e\n', kappa_A);
fprintf('  cond(A^T A)     = %.6e（法方程矩阵条件数）\n', kappa_ATA);
fprintf('  机器精度 eps    = %.4e\n', eps);
% eps：MATLAB 预定义常量，双精度机器精度 ≈ 2.2204e-16
fprintf('  解的相对误差上界 ≈ cond(A) × eps = %.4e\n\n', kappa_A * eps);

% ---- 求解 9 次多项式拟合 ----
% A 是 10×10 方阵，A \ y 退化为方阵方程组（LU分解），而非最小二乘
% 理论上 Ac = y 有精确解（RMS=0），但数值上由于病态性，结果不精确
c_mono = A_mono \ y_data(:);
% c_mono(k+1) = a_k 是 x^k 的系数

% ---- 计算均方根误差 ----
y_fit_mono = A_mono * c_mono;     % 在数据点处的拟合值（应近似等于 y_data）
residual_mono = y_fit_mono - y_data(:);  % 残差向量
% 均方根误差 RMS = sqrt((1/m) Σ residual_i^2) = ||residual||_2 / sqrt(m)
rms_mono = sqrt(mean(residual_mono.^2));
% mean(v.^2) = (1/m) Σ v_i^2（均方值）；sqrt 取开方得 RMS

fprintf('9次多项式拟合结果（幂函数基）：\n');
fprintf('  均方根误差（RMS）= %.4e\n', rms_mono);
fprintf('  （理论值应为 0，但因极高条件数，数值误差放大导致 RMS ≠ 0）\n\n');

% ---- 在绘图点上的拟合值 ----
A_mono_plot = zeros(length(x_plot), n+1);
for k = 0:n
    A_mono_plot(:, k+1) = x_plot(:).^k;
end
fit_mono = A_mono_plot * c_mono;

%% =========================================================
%  第(2)部分：理论分析——为何 RMS 理论为 0 但数值不为 0？
%  =========================================================

fprintf('===== 第(2)部分：理论分析 =====\n\n');

fprintf('【定理】9次多项式对10个数据点的最小二乘拟合等价于 Lagrange 插值\n\n');
fprintf('原因分析：\n');
fprintf('  1. 数据点数 m = 10，多项式次数 n = 9，参数个数 = 10\n');
fprintf('  2. 设计矩阵 A 是 10×10 Vandermonde 方阵\n');
fprintf('  3. 当数据点 x_i 互不相同时，Vandermonde 矩阵可逆\n');
fprintf('  4. 方程组 A·c = y 有唯一解（不再是超定方程组，直接求解）\n');
fprintf('  5. 解 c 使得 φ*(x_i) = y_i 精确成立 → RMS = 0（理论值）\n');
fprintf('  6. 这与过全部10个点的9次 Lagrange 插值多项式完全相同\n\n');

fprintf('【数值现象】为何计算得到的 RMS ≠ 0？\n');
fprintf('  1. 设计矩阵条件数 cond(A) = %.4e\n', kappa_A);
fprintf('  2. 机器精度 eps = %.4e\n', eps);
fprintf('  3. 数值求解 A·c = y 时，解的误差 ≤ cond(A) × eps ≈ %.4e\n', kappa_A * eps);
fprintf('  4. 舍入误差通过极高条件数被放大，使 A·c ≠ y 精确\n');
fprintf('  5. 因此残差 ||A·c - y|| ≠ 0，计算得 RMS = %.4e\n\n', rms_mono);

% ---- 对比：直接用 Lagrange 插值（理论上与9次拟合相同）----
fprintf('对比验证：Lagrange 插值多项式（理论上与9次拟合等价）\n');
% 利用 MATLAB 的 polyfit 函数做插值（n=9次多项式过所有10点）
% polyfit(x, y, n)：多项式拟合，返回从高次到低次排列的系数向量
% 语法：p = polyfit(x, y, n)
% 注意：polyfit 内部使用正交化 QR 方法，比直接用 Vandermonde 矩阵更稳定
p_lag = polyfit(x_data, y_data, 9);
% p_lag：长度为10的向量，p_lag(1)是x^9的系数，p_lag(end)是常数项

% polyval(p, x)：多项式求值
% 语法：y = polyval(p, x)，p 是从高次到低次排列的系数向量
y_lag_data = polyval(p_lag, x_data);
rms_lag = sqrt(mean((y_lag_data - y_data).^2));
fprintf('  polyfit（内置Lagrange插值）RMS = %.4e\n', rms_lag);
fprintf('  （polyfit 使用更稳定的正交化算法，RMS 更接近 0）\n\n');

%% =========================================================
%  第(3)部分：Chebyshev 多项式基 {T_k(x)}，k=0,...,9
%  观察条件数和均方误差的改善情况
%  =========================================================

fprintf('===== 第(3)部分：Chebyshev 基 {T_k : k=0,...,9} =====\n\n');
fprintf('理论背景（教材§3.3(一)）：\n');
fprintf('  Chebyshev 多项式 T_k(t) = cos(k·arccos(t))，t∈[-1,1]\n');
fprintf('  递推关系（教材式3.3.8）：T_0=1, T_1=t\n');
fprintf('  T_{k+1}(t) = 2t·T_k(t) - T_{k-1}(t)\n');
fprintf('  关键性质（教材式3.3.14）：在首一 n 次多项式中，T_n/2^{n-1} 无穷范数最小\n');
fprintf('  作为最小二乘基：设计矩阵的列接近正交 → 条件数远低于 Vandermonde\n\n');

% ---- 将 x 数据映射到 [-1, 1] ----
a = min(x_data);   % 数据 x 的最小值
b = max(x_data);   % 数据 x 的最大值
% 线性映射 t = (2x - (a+b)) / (b-a)，将 [a,b] 映射到 [-1,1]
t_data = (2*x_data - (a+b)) / (b-a);
% 验证：min(t_data) ≈ -1, max(t_data) ≈ 1
fprintf('x 范围 [%.4f, %.4f] 映射到 t 范围 [%.4f, %.4f]（标准[-1,1]）\n\n', ...
    a, b, min(t_data), max(t_data));

% ---- 构造 Chebyshev 设计矩阵 A（10×10）----
% A(i, k+1) = T_k(t_i) = cos(k·arccos(t_i))
% 数值分析含义：用 Chebyshev 多项式替换幂函数作为基，降低设计矩阵条件数
A_cheby = zeros(m, n+1);
for k = 0:n
    % T_k(t) = cos(k · arccos(t))（Chebyshev 多项式的定义式，教材式3.3.5）
    % acos(v)：对向量 v 逐元素取反余弦（结果在 [0,π] 内）
    % k * acos(...)：对 arccos 结果乘以 k（即 k 倍角）
    % cos(...)：取余弦得 T_k(t)
    A_cheby(:, k+1) = cos(k * acos(t_data(:)));
    % 注意：t_data 的元素必须在 [-1,1] 内，acos 才有实数值
end

% ---- 条件数对比 ----
kappa_Acheby = cond(A_cheby);
kappa_ATAcheby = cond(A_cheby' * A_cheby);

fprintf('Chebyshev 设计矩阵 A（%d×%d）：\n', m, n+1);
fprintf('  cond(A_Cheby)         = %.6e\n', kappa_Acheby);
fprintf('  cond(A_Cheby^T A_Cheby) = %.6e\n', kappa_ATAcheby);

fprintf('\n条件数改善对比：\n');
fprintf('  幂函数基 cond(A)      = %.6e\n', kappa_A);
fprintf('  Chebyshev 基 cond(A)  = %.6e\n', kappa_Acheby);
fprintf('  条件数比值（改善倍数）= %.4e\n\n', kappa_A / kappa_Acheby);

% ---- 求解 Chebyshev 基下的9次多项式拟合 ----
c_cheby = A_cheby \ y_data(:);
% 同样是10×10方阵，理论上 RMS=0，Chebyshev 基使其数值上更接近 0

% ---- 计算均方根误差 ----
y_fit_cheby = A_cheby * c_cheby;
rms_cheby = sqrt(mean((y_fit_cheby - y_data(:)).^2));
fprintf('Chebyshev 基拟合 RMS = %.4e\n', rms_cheby);
fprintf('幂函数基拟合 RMS     = %.4e\n', rms_mono);
fprintf('RMS 改善倍数         = %.4e\n\n', rms_mono / max(rms_cheby, 1e-20));

% ---- 在绘图点上的拟合值（Chebyshev 基）----
t_plot = (2*x_plot - (a+b)) / (b-a);   % 绘图点映射到 [-1,1]
A_cheby_plot = zeros(length(x_plot), n+1);
for k = 0:n
    A_cheby_plot(:, k+1) = cos(k * acos(t_plot(:)));
end
fit_cheby = A_cheby_plot * c_cheby;

%% =========================================================
%  汇总对比表
%  =========================================================

fprintf('==============================================\n');
fprintf('  9次多项式拟合：两种基函数全面对比\n');
fprintf('==============================================\n');
fprintf('%-25s  %-14s  %-14s  %-14s\n', '基函数', 'cond(A)', 'cond(A^T A)', 'RMS 误差');
fprintf('%s\n', repmat('-', 1, 74));
fprintf('%-25s  %-14.4e  %-14.4e  %-14.4e\n', ...
    '幂函数基 {x^k}', kappa_A, kappa_ATA, rms_mono);
fprintf('%-25s  %-14.4e  %-14.4e  %-14.4e\n', ...
    'Chebyshev 基 {T_k(t)}', kappa_Acheby, kappa_ATAcheby, rms_cheby);
fprintf('\n理论 RMS = 0（m=n+1 时的精确插值情形）\n');
fprintf('数值 RMS 不为0的根本原因：条件数极大导致舍入误差被放大\n\n');

%% =========================================================
%  绘图
%  =========================================================

figure('Name', '9次多项式拟合：条件数与均方误差分析', 'NumberTitle', 'off', ...
    'Position', [50, 50, 1400, 500]);

% ---- 左图：幂函数基拟合 ----
subplot(1, 3, 1);
plot(x_data, y_data, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', '数据点');
hold on;
plot(x_plot, fit_mono, 'b-', 'LineWidth', 1.5, 'DisplayName', '幂函数基9次拟合');
xlabel('x');  ylabel('y');
title({sprintf('幂函数基 {x^k}'), ...
       sprintf('cond(A)=%.2e', kappa_A), ...
       sprintf('RMS=%.2e', rms_mono)}, 'FontSize', 9);
legend('Location', 'northwest', 'FontSize', 8);
grid on;  hold off;

% ---- 中图：Chebyshev 基拟合 ----
subplot(1, 3, 2);
plot(x_data, y_data, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', '数据点');
hold on;
plot(x_plot, fit_cheby, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Chebyshev基9次拟合');
xlabel('x');  ylabel('y');
title({sprintf('Chebyshev 基 {T_k(t)}'), ...
       sprintf('cond(A)=%.2e', kappa_Acheby), ...
       sprintf('RMS=%.2e', rms_cheby)}, 'FontSize', 9);
legend('Location', 'northwest', 'FontSize', 8);
grid on;  hold off;

% ---- 右图：残差对比（绝对值，数据点处）----
subplot(1, 3, 3);
% semilogy：y 轴为对数刻度（适合显示跨越多个数量级的量）
semilogy(1:m, abs(y_fit_mono  - y_data(:)), 'b-o', 'LineWidth', 1.2, ...
    'MarkerSize', 6, 'DisplayName', sprintf('幂函数基(RMS=%.2e)', rms_mono));
hold on;
semilogy(1:m, abs(y_fit_cheby - y_data(:)), 'r-s', 'LineWidth', 1.2, ...
    'MarkerSize', 6, 'DisplayName', sprintf('Chebyshev基(RMS=%.2e)', rms_cheby));
xlabel('数据点编号 i');
ylabel('|φ*(x_i) - y_i|（对数刻度）');
title('各数据点处的绝对残差对比', 'FontSize', 10);
legend('Location', 'best', 'FontSize', 8);
grid on;  hold off;

sgtitle('9次多项式拟合（m=n+1=10）：幂函数基 vs Chebyshev 基', 'FontSize', 12);

fprintf('图形已生成，请查看弹出的图形窗口。\n\n');
fprintf('程序运行完毕。\n');
