%% =========================================================
%  模板T04：最小二乘数据拟合 —— 混合基函数
%
%  【问题描述】
%  给定数据点 (xᵢ, yᵢ)，用非多项式的基函数组合拟合：
%    φ(x) = a₀φ₀(x) + a₁φ₁(x) + a₂φ₂(x) + a₃φ₃(x)
%
%  常见混合基函数例子（按题目要求修改）：
%    例1：{1, x, cos(x), sin(x)}  →  φ(x) = a₀ + a₁x + a₂cos(x) + a₃sin(x)
%    例2：{1, x, x², eˣ}          →  φ(x) = a₀ + a₁x + a₂x² + a₃eˣ
%    例3：{1, x, ln(x), 1/x}      →  φ(x) = a₀ + a₁x + a₂ln(x) + a₃/x
%
%  【方法与多项式拟合完全一样！】
%  只是设计矩阵 A 的每列变成对应基函数在数据点处的值，
%  然后同样用 A \ y 求解最小二乘问题。
%
%  【优势】
%  合适的基函数能用更少的参数达到更好的拟合效果
%  =========================================================

clear; close all; clc;

%% ===== ① 修改这里：输入你的数据 =====
x_data = [0.000, 0.895, 1.641, 2.512, 3.542, 4.054, 4.602, 5.063, 5.354, 5.617];
y_data = [1.000, 1.803, 3.680, 7.320, 13.59, 17.41, 22.19, 24.89, 26.55, 29.77];
%% =====================================

m = length(x_data);   % 数据点个数

%% ===== ② 修改这里：定义你的基函数 =====
% 方式：在下面的 cell 数组中，每个元素是一个匿名函数（一个基函数）
% 对应名称写在 basis_names 里（用于打印和绘图）

basis_funcs = {
    @(x) ones(size(x)),   % φ₀ = 1（常数）
    @(x) x,               % φ₁ = x（线性）
    @(x) cos(x),          % φ₂ = cos(x)
    @(x) sin(x)           % φ₃ = sin(x)
};
basis_names = {'1', 'x', 'cos(x)', 'sin(x)'};

% ---- 其他常用组合，把上面那段替换成下面任意一段 ----
%
% 多项式+指数：
% basis_funcs = {@(x) ones(size(x)), @(x) x, @(x) x.^2, @(x) exp(x)};
% basis_names = {'1', 'x', 'x^2', 'e^x'};
%
% 只用三角函数（周期数据适用）：
% basis_funcs = {@(x) ones(size(x)), @(x) cos(x), @(x) sin(x), @(x) cos(2*x)};
% basis_names = {'1', 'cos(x)', 'sin(x)', 'cos(2x)'};
%% =========================================

nb = length(basis_funcs);   % 基函数个数（即拟合参数个数）
fprintf('混合基函数最小二乘拟合\n');
fprintf('基函数：{%s}\n', strjoin(basis_names, ', '));
% strjoin(cell, sep) 将 cell 数组的字符串用分隔符连接成一个字符串
fprintf('数据点个数 m=%d，参数个数 nb=%d\n\n', m, nb);

%% ---- 第一步：构造设计矩阵 A（m×nb 矩阵）----
% A(i, k) = φ_{k-1}(xᵢ)，即第 k 列是第 k 个基函数在所有数据点处的值
A = zeros(m, nb);   % 预分配
for k = 1:nb
    A(:, k) = basis_funcs{k}(x_data(:));
    % basis_funcs{k} 取第 k 个基函数（cell 数组用 {} 索引）
    % (x_data(:)) 对列向量形式的 x_data 调用该基函数
    % 结果赋给 A 的第 k 列（A(:, k) 表示第 k 列的所有行）
end

fprintf('设计矩阵 A（%d×%d）条件数：cond(A) = %.4e\n\n', m, nb, cond(A));

%% ---- 第二步：求解最小二乘问题 ----
c = A \ y_data(:);
% 当 m > nb（超定）时，\ 自动求最小二乘解 min ||Ac - y||₂
% c(k) 是第 k 个基函数的系数

fprintf('拟合系数：\n');
for k = 1:nb
    fprintf('  a_%d（%s 的系数）= %.8f\n', k-1, basis_names{k}, c(k));
end

% 拼接拟合函数字符串（便于展示）
fprintf('\n拟合函数：φ(x) = ');
for k = 1:nb
    if k == 1
        fprintf('%.4f·%s', c(k), basis_names{k});
    else
        if c(k) >= 0
            fprintf(' + %.4f·%s', c(k), basis_names{k});
        else
            fprintf(' - %.4f·%s', abs(c(k)), basis_names{k});
        end
    end
end
fprintf('\n\n');

%% ---- 第三步：计算误差 ----
y_fit_data = A * c;   % 在数据点处的拟合值
rms = sqrt(mean((y_fit_data - y_data(:)).^2));   % 均方根误差
fprintf('均方根误差 RMS = %.6f\n', rms);
fprintf('均方误差   L   = %.6f\n\n', rms^2);

%% ---- 第四步：与多项式拟合对比（可选，改为 false 跳过） ----
compare_with_poly = true;

if compare_with_poly
    fprintf('--- 与等参数数目的多项式拟合对比 ---\n');
    n_poly = nb - 1;   % 多项式次数与参数个数相同
    A_poly = zeros(m, nb);
    for k = 0:n_poly
        A_poly(:, k+1) = x_data(:).^k;
    end
    c_poly = A_poly \ y_data(:);
    y_fit_poly = A_poly * c_poly;
    rms_poly = sqrt(mean((y_fit_poly - y_data(:)).^2));
    fprintf('  多项式 P_%d 拟合 RMS = %.6f\n', n_poly, rms_poly);
    fprintf('  混合基函数拟合  RMS = %.6f\n', rms);
    if rms < rms_poly
        fprintf('  结论：混合基函数效果更好（RMS减小 %.2f%%）\n\n', (rms_poly-rms)/rms_poly*100);
    else
        fprintf('  结论：多项式效果更好\n\n');
    end
end

%% ---- 第五步：绘图 ----
x_plot = linspace(min(x_data), max(x_data), 500);

% 在绘图密集点上构造设计矩阵并计算拟合曲线
A_plot = zeros(length(x_plot), nb);
for k = 1:nb
    A_plot(:, k) = basis_funcs{k}(x_plot(:));
end
y_fit_plot = A_plot * c;   % 拟合曲线在密集点的值

figure('Name', '混合基函数最小二乘拟合', 'NumberTitle', 'off', 'Position', [50,50,1100,450]);

subplot(1, 2, 1);
plot(x_data, y_data, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'DisplayName', '数据点');
hold on;
plot(x_plot, y_fit_plot, 'r-', 'LineWidth', 2, 'DisplayName', '混合基拟合');
if compare_with_poly
    % 在绘图点计算多项式拟合曲线
    A_poly_plot = zeros(length(x_plot), nb);
    for k = 0:n_poly
        A_poly_plot(:, k+1) = x_plot(:).^k;
    end
    y_poly_plot = A_poly_plot * c_poly;
    plot(x_plot, y_poly_plot, 'b--', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('P_%d 多项式拟合', n_poly));
end
hold off;
xlabel('x'); ylabel('y');
title(sprintf('拟合结果对比\n混合基RMS=%.4f', rms));
legend('Location', 'northwest');
grid on;

subplot(1, 2, 2);
% 绘制残差（数据点处的拟合误差）
stem(x_data, y_fit_data - y_data(:), 'filled', 'MarkerSize', 5);
% stem：茎叶图，适合展示离散残差
xlabel('x'); ylabel('拟合值 - 实际值');
title(sprintf('残差图（RMS=%.4f）', rms));
grid on;
yline(0, 'k--');   % 在 y=0 处画水平虚线作参考线

sgtitle(sprintf('混合基函数 {%s} 最小二乘拟合', strjoin(basis_names, ', ')));
