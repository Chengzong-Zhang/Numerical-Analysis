%% =========================================================
%  模板T12：多项式拟合条件数分析 + Chebyshev 多项式基
%
%  【问题背景】
%  用 n 次多项式拟合 m 个数据点：
%  (a) 当 m = n+1 时，最小二乘退化为多项式插值（RMS理论=0，但数值上不为0）
%  (b) 幂函数基（Vandermonde 矩阵）随着 n 增大极度病态
%  (c) 改用 Chebyshev 多项式基可显著改善条件数
%
%  【Chebyshev 多项式】
%  定义：T_k(t) = cos(k·arccos(t))，t∈[-1,1]
%  三项递推：
%    T_0(t) = 1
%    T_1(t) = t
%    T_{k+1}(t) = 2t·T_k(t) - T_{k-1}(t)
%
%  优势：在 [-1,1] 上近似正交，设计矩阵条件数远小于 Vandermonde 矩阵
%  使用前需将数据映射到 [-1,1]：t = 2(x-x_min)/(x_max-x_min) - 1
%  =========================================================

clear; close all; clc;

%% ===== ① 修改这里：定义数据 =====
% 使用例3.5.3的10个数据点（可替换为题目数据）
x_data = [0.000, 0.895, 1.641, 2.512, 3.542, 4.054, 4.602, 5.063, 5.354, 5.617];
y_data = [1.000, 1.803, 3.680, 7.320, 13.59, 17.41, 22.19, 24.89, 26.55, 29.77];
n = 9;   % 拟合多项式次数（当 n+1 = m 时是插值问题）
%% =====================================

m = length(x_data);
fprintf('条件数分析：%d 个数据点，%d 次多项式拟合\n\n', m, n);

%% ===================================================
%  第(1)部分：幂函数基（Vandermonde 矩阵）
%  ===================================================
fprintf('===== 幂函数基（Vandermonde 矩阵）=====\n\n');

A_mono = zeros(m, n+1);   % 设计矩阵
for k = 0:n
    A_mono(:, k+1) = x_data(:).^k;
end

kappa_A   = cond(A_mono);
kappa_ATA = cond(A_mono' * A_mono);
fprintf('cond(A)   = %.4e\n', kappa_A);
fprintf('cond(AᵀA) = %.4e（≈ cond(A)²！用法方程组更病态）\n\n', kappa_ATA);

c_mono = A_mono \ y_data(:);   % QR 分解求最小二乘解
y_fit_mono = A_mono * c_mono;
rms_mono = sqrt(mean((y_fit_mono - y_data(:)).^2));
fprintf('均方根误差 RMS = %.2e\n', rms_mono);
if rms_mono < 1e-6
    fprintf('【注意】RMS 理论上应为0（m=n+1 是插值），\n');
    fprintf('但数值上不为0（%.2e），这就是数值病态的直接体现！\n\n', rms_mono);
end

%% ===================================================
%  第(2)部分：为什么 m=n+1 时 RMS 理论为 0 但数值不为 0？
%  ===================================================
fprintf('===== 理论分析 =====\n\n');
fprintf('当参数个数 n+1 = 数据点个数 m 时：\n');
fprintf('  - 最小二乘问题退化为恰定方程组，即多项式插值\n');
fprintf('  - 理论上解是唯一的，残差为零（RMS=0）\n');
fprintf('  - 但 Vandermonde 矩阵的条件数约 %.2e\n', kappa_A);
fprintf('  - 机器精度 eps ≈ %.2e\n', eps);
fprintf('  - 数值求解误差约 ≈ cond(A) × eps ≈ %.2e\n', kappa_A*eps);
fprintf('  - 这就是数值 RMS 不为零的原因（舍入误差被放大）\n\n');

%% ===================================================
%  第(3)部分：Chebyshev 多项式基
%  ===================================================
fprintf('===== Chebyshev 多项式基 =====\n\n');

% 首先将 x 数据映射到 [-1,1]（Chebyshev 的标准区间）
x_min = min(x_data);
x_max = max(x_data);
t_data = 2*(x_data - x_min)/(x_max - x_min) - 1;
% 映射公式：t = 2(x - x_min)/(x_max - x_min) - 1
% 当 x=x_min 时 t=-1，x=x_max 时 t=1

fprintf('数据 x 映射到 t∈[-1,1]：t = 2(x-%.3f)/%.3f - 1\n\n', x_min, x_max-x_min);

% 构造 Chebyshev 基设计矩阵
A_cheb = zeros(m, n+1);
for k = 0:n
    A_cheb(:, k+1) = chebyshev_T(k, t_data(:));
    % chebyshev_T(k, t)：第 k 个 Chebyshev 多项式（定义在文件末尾）
end

kappa_cheb   = cond(A_cheb);
kappa_cheb_ATA = cond(A_cheb' * A_cheb);
fprintf('Chebyshev 基设计矩阵：\n');
fprintf('  cond(A_cheb)   = %.4e（远小于 Vandermonde！）\n', kappa_cheb);
fprintf('  cond(AᵀA)      = %.4e\n\n', kappa_cheb_ATA);

% 改善比例
fprintf('条件数改善比例：幂函数基/Chebyshev基 = %.2e 倍\n\n', kappa_A/kappa_cheb);

% 求解
c_cheb = A_cheb \ y_data(:);
y_fit_cheb = A_cheb * c_cheb;
rms_cheb = sqrt(mean((y_fit_cheb - y_data(:)).^2));
fprintf('Chebyshev 基 RMS = %.2e\n\n', rms_cheb);

%% ---- 汇总对比 ----
fprintf('==============================================\n');
fprintf('  条件数与误差汇总\n');
fprintf('==============================================\n');
fprintf('%-20s  %-15s  %-15s\n', '方法', 'cond(A)', 'RMS 误差');
fprintf('%s\n', repmat('-', 1, 52));
fprintf('%-20s  %-15.4e  %-15.4e\n', '幂函数基（Vandermonde）', kappa_A,    rms_mono);
fprintf('%-20s  %-15.4e  %-15.4e\n', 'Chebyshev 基',           kappa_cheb, rms_cheb);
fprintf('\n');
fprintf('结论：Chebyshev 基的条件数比幂函数基小 %.0e 倍，\n', kappa_A/kappa_cheb);
fprintf('      数值稳定性大幅提升，与 Legendre 基类似\n\n');

%% ---- 绘图对比 ----
x_plot = linspace(x_min, x_max, 500);
t_plot = 2*(x_plot - x_min)/(x_max - x_min) - 1;   % 对应的 t 值

% 幂函数基拟合曲线
A_mono_plot = zeros(length(x_plot), n+1);
for k = 0:n
    A_mono_plot(:, k+1) = x_plot(:).^k;
end
y_mono_plot = A_mono_plot * c_mono;

% Chebyshev 基拟合曲线
A_cheb_plot = zeros(length(x_plot), n+1);
for k = 0:n
    A_cheb_plot(:, k+1) = chebyshev_T(k, t_plot(:));
end
y_cheb_plot = A_cheb_plot * c_cheb;

figure('Name', '幂函数基 vs Chebyshev 基', 'NumberTitle', 'off', 'Position', [50,50,1200,500]);

subplot(1, 2, 1);
plot(x_data, y_data, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'DisplayName', '数据点');
hold on;
plot(x_plot, y_mono_plot, 'b-', 'LineWidth', 1.5, 'DisplayName', ...
    sprintf('幂函数基（cond=%.2e）', kappa_A));
plot(x_plot, y_cheb_plot, 'r--', 'LineWidth', 1.5, 'DisplayName', ...
    sprintf('Chebyshev基（cond=%.2e）', kappa_cheb));
hold off;
xlabel('x'); ylabel('y');
title('拟合曲线对比');
legend('Location', 'northwest');
grid on;

subplot(1, 2, 2);
% 残差图（在数据点处的误差）
stem(x_data, abs(y_fit_mono - y_data(:)), 'b', 'DisplayName', ...
    sprintf('幂函数基残差（RMS=%.2e）', rms_mono));
hold on;
stem(x_data, abs(y_fit_cheb - y_data(:)), 'r', 'DisplayName', ...
    sprintf('Chebyshev基残差（RMS=%.2e）', rms_cheb));
% stem：茎叶图，适合显示离散误差
hold off;
xlabel('x'); ylabel('|残差|');
title('残差绝对值对比');
legend;
grid on;

sgtitle(sprintf('n=%d 次多项式拟合：幂函数基 vs Chebyshev 基条件数对比', n));


%% =========================================================
%  局部函数（必须放在文件末尾）
%  =========================================================

function vals = chebyshev_T(k, t)
    % 计算第 k 个 Chebyshev 多项式 T_k(t) 在 t 处的值
    %
    % 定义：T_k(t) = cos(k·arccos(t))，t∈[-1,1]
    % 三项递推（数值更稳定）：
    %   T_0(t) = 1
    %   T_1(t) = t
    %   T_{k+1}(t) = 2t·T_k(t) - T_{k-1}(t)
    %
    % 输入：k - 阶数，t - [-1,1] 内的点（可为向量）
    % 输出：vals - T_k(t) 在各点的值

    if k == 0
        vals = ones(size(t));       % T_0(t) = 1
    elseif k == 1
        vals = t;                   % T_1(t) = t
    else
        % 三项递推（和 Legendre 类似，但系数不同）
        T_prev = ones(size(t));     % T_0(t) = 1
        T_curr = t;                 % T_1(t) = t
        for j = 1 : k-1
            T_next = 2*t.*T_curr - T_prev;   % T_{j+1} = 2t·T_j - T_{j-1}
            % 注意：系数是2，比 Legendre 的 (2j+1)/(j+1) 更简单！
            T_prev = T_curr;
            T_curr = T_next;
        end
        vals = T_curr;
    end
end
