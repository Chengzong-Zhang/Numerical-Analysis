%% =========================================================
%  模板T03：最小二乘数据拟合 —— 多项式基（Vandermonde 设计矩阵）
%
%  【问题描述】
%  给定 m 个数据点 (xᵢ, yᵢ)，i=1,...,m
%  用 n 次多项式 φ(x) = a₀ + a₁x + a₂x² + ... + aₙxⁿ 拟合，
%  使均方误差最小：
%    L = Σᵢ [yᵢ - φ(xᵢ)]² → min（这是离散最小二乘）
%
%  【与最佳平方逼近的区别】
%  - 最佳平方逼近：针对已知连续函数，用积分定义内积，求精确逼近
%  - 最小二乘拟合：针对离散数据点，用求和定义内积，允许不过点
%    （数据有误差时不应用插值，用最小二乘更合理）
%
%  【矩阵形式】
%  设计矩阵 A（m×(n+1)）：A(i, k+1) = xᵢᵏ（第 i 行第 k 列是 xᵢ 的 k 次方）
%  最小二乘问题等价于：min ||Ac - y||₂²
%  正规方程（法方程组）：(AᵀA)c = Aᵀy
%  MATLAB 中用 A \ y 直接求解（QR 分解，比先算 AᵀA 更稳定）
%  =========================================================

clear; close all; clc;

%% ===== ① 修改这里：输入你的数据和问题参数 =====
% 数据点（改成题目给的数据）
x_data = [0.000, 0.895, 1.641, 2.512, 3.542, 4.054, 4.602, 5.063, 5.354, 5.617];
y_data = [1.000, 1.803, 3.680, 7.320, 13.59, 17.41, 22.19, 24.89, 26.55, 29.77];

n_list = [1, 2, 3, 4];   % 要测试的多项式次数（可以只写一个，如 n_list = 2）
%% =====================================================

m = length(x_data);   % 数据点个数
fprintf('数据点个数 m = %d\n', m);
fprintf('x 范围：[%.3f, %.3f]\n', min(x_data), max(x_data));
fprintf('y 范围：[%.3f, %.3f]\n\n', min(y_data), max(y_data));

% 绘图用的密集点（不参与计算，只用于可视化拟合曲线）
x_plot = linspace(min(x_data), max(x_data), 500);

% 存储各次数的结果
rms_list   = zeros(size(n_list));   % 均方根误差
cond_list  = zeros(size(n_list));   % 条件数
fit_curves = cell(size(n_list));    % 绘图用的拟合曲线值

for idx = 1 : length(n_list)
    n = n_list(idx);
    fprintf('===== n=%d 次多项式拟合（共 %d 个参数）=====\n', n, n+1);

    %% ---- 第一步：构造设计矩阵 A（m×(n+1) 矩阵）----
    % A(i, k+1) = xᵢᵏ，其中 i=1,...,m 是数据点下标，k=0,...,n 是幂次
    % 第 k+1 列就是基函数 φ_k(x) = x^k 在所有数据点处的取值
    % 这种矩阵结构叫做 Vandermonde 矩阵
    A = zeros(m, n+1);   % 预分配 m×(n+1) 全零矩阵
    for k = 0:n
        A(:, k+1) = x_data(:) .^ k;
        % x_data(:) 将 x_data 变成列向量（: 拉平所有元素）
        % .^ k 对列向量每个元素求 k 次方（逐元素幂）
        % A(:, k+1) 表示 A 的第 k+1 整列（: 表示所有行）
    end

    %% ---- 第二步：求解最小二乘问题 ----
    % A \ y（左除）：当 A 是 m×(n+1) 且 m > n+1 时（超定方程组），
    % MATLAB 自动用 QR 分解求最小二乘解（最小化 ||Ac-y||₂）
    % 这比先算 (A'*A) 再求逆更数值稳定（避免条件数平方）
    c = A \ y_data(:);
    % y_data(:) 确保 y_data 是列向量
    % c 是系数向量：c(1)=a₀, c(2)=a₁, ..., c(n+1)=aₙ

    fprintf('  拟合系数：');
    for k = 0:n
        fprintf('a_%d=%.4f  ', k, c(k+1));
    end
    fprintf('\n');

    %% ---- 第三步：计算条件数 ----
    kappa_A   = cond(A);        % 设计矩阵 A 的条件数
    kappa_ATA = cond(A' * A);   % 法方程矩阵 AᵀA 的条件数（≈ κ(A)²）
    % A'：矩阵 A 的转置（实矩阵的转置就是各行列互换）
    cond_list(idx) = kappa_A;
    fprintf('  cond(A) = %.4e，cond(AᵀA) = %.4e（注意：约是 cond(A) 的平方！）\n', ...
        kappa_A, kappa_ATA);

    %% ---- 第四步：计算拟合误差 ----
    y_fit = A * c;   % 在数据点处的拟合值（矩阵-向量乘法）
    % y_fit(i) = Σ_k c(k+1)*xᵢᵏ = φ*(xᵢ)

    % 均方根误差（Root Mean Square Error）：衡量拟合质量的核心指标
    % RMS = sqrt( (1/m)·Σᵢ(φ*(xᵢ) - yᵢ)² ) = ||Ac-y||₂ / sqrt(m)
    rms = sqrt(mean((y_fit - y_data(:)).^2));
    % mean(v) 求向量 v 的算术平均值；.^2 逐元素平方
    rms_list(idx) = rms;
    fprintf('  均方根误差 RMS = %.6f\n', rms);
    fprintf('  均方误差   L   = %.6f（= RMS²）\n\n', rms^2);

    %% ---- 第五步：计算绘图用的拟合曲线 ----
    % 为绘图点构造同样结构的设计矩阵（行数变为 length(x_plot)）
    A_plot = zeros(length(x_plot), n+1);
    for k = 0:n
        A_plot(:, k+1) = x_plot(:) .^ k;
    end
    fit_curves{idx} = A_plot * c;   % 拟合曲线在绘图点处的值
    % fit_curves 是 cell 数组，用 {} 索引，每个元素是一个向量
end

%% ---- 第六步：结果汇总 ----
fprintf('============== 汇总 ==============\n');
fprintf('%-15s  %-12s  %-12s\n', '拟合方法', 'RMS 误差', 'cond(A)');
fprintf('%s\n', repmat('-', 1, 43));   % repmat('-', 1, 43) 生成 43 个 '-'
for idx = 1:length(n_list)
    fprintf('%-15s  %-12.6f  %-12.4e\n', ...
        sprintf('P_%d (%d参数)', n_list(idx), n_list(idx)+1), ...
        rms_list(idx), cond_list(idx));
end
fprintf('\n');

%% ---- 第七步：绘图 ----
% 图1：各次数拟合曲线对比
figure('Name', '最小二乘多项式拟合', 'NumberTitle', 'off', 'Position', [50,50,1200,800]);

for idx = 1:length(n_list)
    n = n_list(idx);
    subplot(2, ceil(length(n_list)/2), idx);   % 自动计算子图布局

    % 绘制原始数据点
    plot(x_data, y_data, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', ...
        'DisplayName', '数据点');
    % 'ro'：红色(r)圆圈(o)；'MarkerSize'：标记大小；'MarkerFaceColor'：填充色
    hold on;
    % 绘制拟合曲线
    plot(x_plot, fit_curves{idx}, 'b-', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('P_%d 拟合', n));
    hold off;
    xlabel('x'); ylabel('y');
    title(sprintf('n=%d 次多项式拟合\nRMS=%.4f, cond(A)=%.2e', ...
        n, rms_list(idx), cond_list(idx)));
    legend('Location', 'northwest');
    grid on;
end
sgtitle('多项式最小二乘拟合（各次数对比）');

% 图2：RMS 误差条形图
figure('Name', 'RMS误差对比', 'NumberTitle', 'off');
bar(1:length(n_list), rms_list, 'FaceColor', [0.3, 0.6, 0.9]);
% bar(x, y) 绘制竖直条形图；'FaceColor'：条形填充色（RGB归一化值）
labels = arrayfun(@(n) sprintf('P_%d', n), n_list, 'UniformOutput', false);
% arrayfun 对数组每个元素调用函数；'UniformOutput', false 表示返回 cell 数组
set(gca, 'XTickLabel', labels, 'XTick', 1:length(n_list));
% set(gca, ...) 设置当前坐标轴属性；gca = get current axes
for k = 1:length(n_list)
    text(k, rms_list(k)*1.05, sprintf('%.4f', rms_list(k)), ...
        'HorizontalAlignment', 'center', 'FontSize', 9);
    % text(x, y, str) 在坐标 (x,y) 处显示文字
end
xlabel('拟合方法'); ylabel('RMS 误差');
title('各次数多项式拟合的 RMS 误差对比');
grid on;
