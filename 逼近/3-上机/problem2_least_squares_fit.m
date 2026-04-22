%% =========================================================
%  problem2_least_squares_fit.m
%  功能：利用例3.5.3的数据求最小线性二乘拟合
%
%  题目两个部分：
%  (1) 在多项式函数空间 P_n 中求最小二乘拟合（测试 n=1,2,3,4）
%  (2) 改用函数集 {1, x, cos(x), sin(x)} 作为基函数进行拟合
%
%  【数值分析背景：线性最小二乘法（教材§3.5.1）】
%  给定 m 个数据点 (x_i, y_i)，i=1,...,m，用函数：
%    φ(x) = Σ_{k=0}^n c_k φ_k(x)
%  拟合数据，使均方误差（损失函数）最小（教材式3.5.1）：
%    L(c) = (1/m) Σ_{i=1}^m [φ(x_i) - y_i]^2 → min
%
%  设计矩阵 A（m×(n+1) 矩阵）：A_{ik} = φ_{k-1}(x_i)
%  将极值条件化为线性方程组（法方程组，教材式3.5.2→3.5.4）：
%    Φ c = b，其中 Φ_{jk} = (φ_j, φ_k) = Σ_i φ_j(x_i)φ_k(x_i)（离散内积）
%           b_j = (φ_j, y) = Σ_i φ_j(x_i) y_i
%  等价矩阵形式：(A^T A) c = A^T y
%
%  MATLAB 实现：用 A\y 直接求解（QR分解，比先算 A^T A 更稳定）
%
%  例3.5.3数据（10个点）：
%    x = [0.000, 0.895, 1.641, 2.512, 3.542, 4.054, 4.602, 5.063, 5.354, 5.617]
%    y = [1.000, 1.803, 3.680, 7.320, 13.59, 17.41, 22.19, 24.89, 26.55, 29.77]
%  =========================================================

clear; close all; clc;

fprintf('==============================================\n');
fprintf('  第三章上机作业2：最小二乘数据拟合\n');
fprintf('  利用例3.5.3数据\n');
fprintf('==============================================\n\n');

%% =========================================================
%  数据（例3.5.3）
%  =========================================================

% 原始测量数据（10个点）
% 数值分析含义：实验数据含测量误差，不能用插值（会保留所有误差），
% 应用最小二乘法找最能反映数据规律的函数
x_data = [0.000, 0.895, 1.641, 2.512, 3.542, 4.054, 4.602, 5.063, 5.354, 5.617];
y_data = [1.000, 1.803, 3.680, 7.320, 13.59, 17.41, 22.19, 24.89, 26.55, 29.77];

% length(v)：向量 v 的元素个数（最长维的大小）
m = length(x_data);   % m = 10，数据点个数
fprintf('例3.5.3 数据：m = %d 个数据点\n', m);
fprintf('x 范围：[%.3f, %.3f]\n', min(x_data), max(x_data));
fprintf('y 范围：[%.3f, %.3f]\n\n', min(y_data), max(y_data));
% min(v)：向量 v 的最小元素；max(v)：最大元素

% 用于绘图的密集点（不参与计算，只用于可视化）
% linspace(a, b, N)：在 [a,b] 内均匀生成 N 个点（含端点）
x_plot = linspace(min(x_data), max(x_data), 500);

%% =========================================================
%  第(1)部分：多项式空间 P_n 中的最小二乘拟合
%  基函数 φ_k(x) = x^k，k=0,...,n
%  对 n=1,2,3,4 分别计算
%  =========================================================

fprintf('===== 第(1)部分：多项式函数空间 P_n 中的最小二乘拟合 =====\n\n');
fprintf('模型：φ(x) = a_0 + a_1·x + a_2·x^2 + ... + a_n·x^n\n\n');
fprintf('理论背景（教材式3.5.2）：\n');
fprintf('  设计矩阵 A（m×(n+1)）：A(i,k+1) = x_i^k\n');
fprintf('  法方程组：(A^T A) c = A^T y\n');
fprintf('  均方误差：L = (1/m)||Ac - y||_2^2\n\n');

n_tests = [1, 2, 3, 4];          % 测试的多项式次数
rms_poly = zeros(size(n_tests));  % 存储各 n 的均方根误差
fit_poly = cell(size(n_tests));   % 存储各 n 在绘图点上的拟合值

for idx = 1 : length(n_tests)
    n = n_tests(idx);
    fprintf('--- n = %d 次多项式拟合 ---\n', n);

    % ---- 构造设计矩阵 A（m×(n+1) 矩阵）----
    % A(i, k+1) = x_i^k，i=1,...,m（数据点下标）; k=0,...,n（基函数下标）
    % 数值分析含义：A 的第 k+1 列是第 k 个基函数 φ_k(x)=x^k 在所有数据点处的取值向量
    % 这是一个 Vandermonde 型矩阵（各列是 x 的不同次幂）
    A = zeros(m, n+1);   % zeros(m,n)：预分配 m×n 全零矩阵
    for k = 0:n
        % x_data(:).^k：x_data 先转为列向量（(:)操作），再逐元素取 k 次方
        % .^ 是逐元素幂运算符（区别于矩阵幂 ^）
        A(:, k+1) = x_data(:).^k;
        % A(:, k+1)：取 A 的第 k+1 整列（: 表示所有行）
    end

    % ---- 求解最小二乘问题 ----
    % 方法：A \ y（左除，当 A 不是方阵时 MATLAB 自动用 QR 分解求最小二乘解）
    % 语法：c = A \ b
    % 数值分析：QR分解法比先算 A^T A 再求逆更稳定（避免平方条件数 κ(A^T A)=κ(A)^2）
    % 与法方程组 (A^T A)c = A^T y 等价，但数值精度更高
    c = A \ y_data(:);
    % y_data(:)：确保 y_data 是列向量（(:) 将任意形状变为列向量）
    % c 是系数向量，c(k+1) = a_k 是 x^k 的系数

    % ---- 输出设计矩阵的条件数 ----
    % cond(A)：设计矩阵 A 的 2-范数条件数
    % 数值分析意义：法方程组矩阵 A^T A 的条件数 = cond(A)^2（条件数平方！）
    % 这说明用法方程组比直接用 A 的 QR 分解更病态
    kappa_A   = cond(A);
    kappa_ATA = cond(A' * A);
    % A'：矩阵 A 的转置（若 A 含复数则为共轭转置；本题 A 为实矩阵）
    fprintf('  设计矩阵 A（%d×%d）条件数 cond(A) = %.4e\n', m, n+1, kappa_A);
    fprintf('  法方程矩阵 A^TA 条件数 = %.4e（≈ cond(A)^2）\n', kappa_ATA);

    % ---- 计算在数据点处的拟合值 ----
    y_fit_data = A * c;
    % A * c：矩阵-向量乘法，结果是拟合函数在各数据点处的值
    % y_fit_data(i) = Σ_k c(k+1) * x_i^k = φ*(x_i)

    % ---- 均方根误差（RMS Error） ----
    % RMS = sqrt((1/m) Σ(φ*(x_i)-y_i)^2) = ||Ac-y||_2 / sqrt(m)
    % 数值分析含义：衡量拟合函数与数据的平均偏差（均方误差的开方）
    % 注意：均方误差 L = RMS^2，教材用的是 L，这里输出 RMS 更直观
    rms_poly(idx) = sqrt(mean((y_fit_data - y_data(:)).^2));
    % mean(v)：向量 v 的算术平均值 = (1/m)Σv_i
    % .^2：逐元素平方；(v1-v2).^2：对应元素之差的平方
    fprintf('  均方根误差（RMS）= %.6f\n', rms_poly(idx));
    fprintf('  均方误差（L）= %.6f\n\n', rms_poly(idx)^2);

    % ---- 在绘图密集点上计算拟合值 ----
    % 为绘图构造密集点上的设计矩阵（同样结构，行数变为 length(x_plot)）
    A_plot = zeros(length(x_plot), n+1);
    for k = 0:n
        A_plot(:, k+1) = x_plot(:).^k;
    end
    fit_poly{idx} = A_plot * c;
    % fit_poly{idx}：cell 数组用 {} 索引，存储第 idx 个拟合结果（列向量）
end

%% =========================================================
%  第(2)部分：函数集 {1, x, cos(x), sin(x)} 拟合
%  这是一个"非多项式"线性逼近问题（函数不再是多项式）
%  =========================================================

fprintf('===== 第(2)部分：基函数集 {1, x, cos(x), sin(x)} 拟合 =====\n\n');
fprintf('理论背景（教材§3.5.1 推广形式）：\n');
fprintf('  基函数：φ_0=1, φ_1=x, φ_2=cos(x), φ_3=sin(x)\n');
fprintf('  拟合模型：f(x) = a_0 + a_1·x + a_2·cos(x) + a_3·sin(x)\n');
fprintf('  线性无关性：这四个函数在数据点集上满足 Haar 条件（教材定理3.5.1）\n');
fprintf('  法方程组同样为 (A^T A)c = A^T y，只是 A 的列不再是幂函数\n\n');

% ---- 构造混合基函数设计矩阵 A（m×4 矩阵）----
% 每列对应一个基函数在所有数据点处的取值
% A(:,1) = φ_0(x_i) = 1        （全1列向量）
% A(:,2) = φ_1(x_i) = x_i      （x 值向量）
% A(:,3) = φ_2(x_i) = cos(x_i) （逐元素余弦）
% A(:,4) = φ_3(x_i) = sin(x_i) （逐元素正弦）
A_mix = [ones(m, 1), x_data(:), cos(x_data(:)), sin(x_data(:))];
% ones(m,1)：m×1 全1列向量（对应常数基函数 φ_0=1）
% x_data(:)：列向量化（确保是列向量）
% cos(v)：对向量 v 逐元素取余弦（MATLAB 内置函数，支持向量化）
% [A, B, C, D]：水平拼接四个列向量为 m×4 矩阵

fprintf('设计矩阵 A 的结构：\n');
fprintf('  A = [1, x_i, cos(x_i), sin(x_i)]（%d×4 矩阵）\n', m);
kappa_mix = cond(A_mix);
fprintf('  条件数 cond(A) = %.4e\n', kappa_mix);

% ---- 用左除求最小二乘解 ----
c_mix = A_mix \ y_data(:);
% A 是 10×4 的过定方程组（超定，m>n），\ 使用 QR 分解求最小范数最小二乘解

fprintf('\n拟合系数：\n');
fprintf('  a_0（常数项 1 的系数）   = %10.6f\n', c_mix(1));
fprintf('  a_1（x 的系数）         = %10.6f\n', c_mix(2));
fprintf('  a_2（cos(x) 的系数）    = %10.6f\n', c_mix(3));
fprintf('  a_3（sin(x) 的系数）    = %10.6f\n', c_mix(4));
fprintf('\n拟合函数：f(x) = %.4f + %.4f·x + %.4f·cos(x) + %.4f·sin(x)\n', ...
    c_mix(1), c_mix(2), c_mix(3), c_mix(4));
% ... 是 MATLAB 的续行符，表示一条语句跨多行书写

% ---- 计算拟合误差 ----
y_fit_mix = A_mix * c_mix;    % 在数据点处的拟合值
rms_mix = sqrt(mean((y_fit_mix - y_data(:)).^2));  % 均方根误差
fprintf('\n均方根误差（RMS）= %.6f\n', rms_mix);
fprintf('均方误差（L）    = %.6f\n\n', rms_mix^2);

% ---- 在绘图点上的拟合值 ----
A_mix_plot = [ones(length(x_plot),1), x_plot(:), cos(x_plot(:)), sin(x_plot(:))];
fit_mix = A_mix_plot * c_mix;

%% =========================================================
%  汇总对比表
%  =========================================================

fprintf('==============================================\n');
fprintf('  各拟合方法均方根误差（RMS）汇总对比\n');
fprintf('==============================================\n');
fprintf('%-35s  %-14s  %-14s\n', '拟合方法', 'RMS 误差', '均方误差 L');
fprintf('%s\n', repmat('-', 1, 68));

for idx = 1:length(n_tests)
    n = n_tests(idx);
    fprintf('%-35s  %-14.6f  %-14.6f\n', ...
        sprintf('多项式 P_%d（%d 个参数）', n, n+1), ...
        rms_poly(idx), rms_poly(idx)^2);
end
fprintf('%-35s  %-14.6f  %-14.6f\n', ...
    '混合基 {1, x, cos(x), sin(x)}', rms_mix, rms_mix^2);
fprintf('\n');

% 对比教材例3.5.3的参考结果（n=2次多项式）
fprintf('教材例3.5.3参考值（P_2拟合）：\n');
fprintf('  a_2=0.7542712, a_1=0.98752785, a_0=0.53554452\n');
fprintf('  L = 0.17507558（均方误差）\n\n');

%% =========================================================
%  绘图
%  =========================================================

% ---- 图1：多项式拟合结果（2×2 子图）----
figure('Name', '多项式最小二乘拟合（例3.5.3数据）', 'NumberTitle', 'off', ...
    'Position', [50, 50, 1200, 800]);

for idx = 1:length(n_tests)
    n = n_tests(idx);
    % subplot(2, 2, idx)：2行2列子图，第 idx 个
    subplot(2, 2, idx);

    % 绘制原始数据点（红色圆圈）
    % 'ro'：红色(r)圆圈(o)标记；'MarkerSize'：标记大小；'MarkerFaceColor'：填充颜色
    plot(x_data, y_data, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', ...
        'DisplayName', '数据点');
    hold on;

    % 绘制拟合曲线（蓝色实线）
    plot(x_plot, fit_poly{idx}, 'b-', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('P_%d 拟合曲线', n));

    xlabel('x');  ylabel('y');
    title({sprintf('多项式拟合 P_%d（次数 n=%d，参数个数 %d）', n, n, n+1), ...
           sprintf('RMS = %.4f', rms_poly(idx))}, 'FontSize', 9);
    legend('Location', 'northwest', 'FontSize', 8);
    grid on;
    hold off;
end
% sgtitle：整个 figure 的总标题
sgtitle('例3.5.3数据：多项式空间 P_n 中的最小二乘拟合', 'FontSize', 13);
print(gcf, 'p2_fig1', '-dpng', '-r150');   % 保存图1到磁盘（150 DPI）

% ---- 图2：所有方法汇总对比 ----
figure('Name', '最小二乘拟合方法汇总对比', 'NumberTitle', 'off', ...
    'Position', [100, 100, 1200, 500]);

subplot(1, 2, 1);
% 绘制数据点
plot(x_data, y_data, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'DisplayName', '数据点');
hold on;

% 颜色列表（按顺序分配给各 n）
colors = {'b', 'r', 'm', 'g'};
for idx = 1:length(n_tests)
    plot(x_plot, fit_poly{idx}, [colors{idx}, '-'], 'LineWidth', 1.5, ...
        'DisplayName', sprintf('P_%d 拟合(RMS=%.3f)', n_tests(idx), rms_poly(idx)));
end
% 也绘制混合基拟合
plot(x_plot, fit_mix, 'c--', 'LineWidth', 2, ...
    'DisplayName', sprintf('混合基(RMS=%.3f)', rms_mix));

xlabel('x');  ylabel('y');
title('各拟合方法曲线对比', 'FontSize', 11);
legend('Location', 'northwest', 'FontSize', 8);
grid on;  hold off;

subplot(1, 2, 2);
% 绘制各方法 RMS 误差的条形图
all_rms = [rms_poly, rms_mix];
all_labels = {'P_1', 'P_2', 'P_3', 'P_4', '混合基'};
% bar(x, y)：绘制竖直条形图
bar(1:5, all_rms, 'FaceColor', [0.3, 0.6, 0.9]);
% 'FaceColor'：条形颜色（RGB 归一化值）
set(gca, 'XTickLabel', all_labels, 'XTick', 1:5);
% set(gca, ...)：设置当前坐标轴(gca=get current axes)属性
% 'XTickLabel'：设置 x 轴刻度标签；'XTick'：设置刻度位置
xlabel('拟合方法');
ylabel('均方根误差（RMS）');
title('各方法 RMS 误差对比', 'FontSize', 11);
% 在每个条形上方显示数值
for k = 1:5
    % text(x, y, str)：在坐标 (x,y) 处显示字符串
    % 'HorizontalAlignment','center'：水平居中对齐
    % 'VerticalAlignment','bottom'：垂直底部对齐（数字在条形顶部上方）
    text(k, all_rms(k), sprintf('%.4f', all_rms(k)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9);
end
grid on;

sgtitle('例3.5.3数据最小二乘拟合：多项式基 vs 混合基', 'FontSize', 13);
print(gcf, 'p2_fig2', '-dpng', '-r150');   % 保存图2到磁盘（150 DPI）

fprintf('图形已保存：p2_fig1.png, p2_fig2.png\n\n');
fprintf('程序运行完毕。\n');
