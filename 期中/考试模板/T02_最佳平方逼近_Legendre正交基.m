% =========================================================
%  模板T02：最佳平方逼近 —— Legendre 正交基
%
%  【问题描述】
%  求函数 f(x) 在区间 [0,1] 上关于 Legendre 正交多项式基的
%  n 次最佳平方逼近：
%    φ*(x) = Σ_{k=0}^n c_k · P_k(2x-1)
%
%  【为什么用 Legendre 基？】
%  Legendre 多项式 P_k(t) 在 [-1,1] 上满足正交性：
%    ∫_{-1}^1 P_m(t)·P_k(t) dt = 2/(2k+1) · δ_{mk}
%  变量替换 t=2x-1 后在 [0,1] 上：
%    ∫_0^1 P_m(2x-1)·P_k(2x-1) dx = δ_{mk}/(2k+1)
%
%  因此 Gram 矩阵是对角阵！条件数=1，数值完全稳定！
%  系数公式直接给出（无需解方程组）：
%    c_k = (f, φ_k) / (φ_k, φ_k) = (2k+1) · ∫_0^1 f(x)·P_k(2x-1) dx
%
%  【三项递推关系（重要！）】
%    P_0(t) = 1
%    P_1(t) = t
%    P_{k+1}(t) = (2k+1)/(k+1)·t·P_k(t) - k/(k+1)·P_{k-1}(t)
%  =========================================================

clear; close all; clc;

%% ===== ① 修改这里：定义你的问题 =====
f = @(x) exp(x);    % 被逼近函数，改成题目要求的函数
n = 5;               % 逼近次数（修改这里）
% 注意：本模板针对 [0,1] 区间；其他区间需要调整变量替换
%% =====================================

fprintf('最佳平方逼近（Legendre 正交基）：n=%d 次，区间[0,1]\n\n', n);

%% ---- 第一步：计算各 Legendre 展开系数 ----
% 由于 Gram 矩阵是对角阵，系数可逐项独立计算，无需解方程组！
%
% c_k = (f, φ_k) / (φ_k, φ_k)
%      = ∫_0^1 f(x)·P_k(2x-1) dx  /  (1/(2k+1))
%      = (2k+1) · ∫_0^1 f(x)·P_k(2x-1) dx
%
% 这是最大的优势：每个系数独立计算，不受其他系数影响！

c = zeros(n+1, 1);   % 预分配系数向量，c(k+1) 对应 P_k 的系数
for k = 0:n
    % integral 计算内积 (f, φ_k) = ∫_0^1 f(x)·P_k(2x-1) dx
    % legendre_on_01(k, x)：计算 P_k(2x-1) 在点 x 处的值（见文件末尾的局部函数）
    inner_fk = integral(@(x) f(x) .* legendre_on_01(k, x), 0, 1);

    norm_sq_k = 1 / (2*k + 1);   % (φ_k, φ_k) = ∫_0^1 P_k²(2x-1) dx = 1/(2k+1)

    c(k+1) = inner_fk / norm_sq_k;   % c_k = 内积 / 模长平方
    % 等价写法：c(k+1) = (2*k+1) * inner_fk;

    fprintf('  c_%d = %12.8f\n', k, c(k+1));
end
fprintf('\n');
fprintf('Gram 矩阵条件数 = 1（正交基，完全稳定！）\n\n');

%% ---- 第二步：在密集点上求逼近多项式的值 ----
x_plot = linspace(0, 1, 500);
y_exact = f(x_plot);

% φ*(x) = Σ_{k=0}^n c_k · P_k(2x-1)
phi_plot = zeros(size(x_plot));   % 初始化为零向量
for k = 0:n
    % legendre_on_01(k, x_plot) 计算 P_k(2x-1) 在所有绘图点的值（向量）
    % c(k+1) * ... 是标量乘以向量（MATLAB 自动广播）
    phi_plot = phi_plot + c(k+1) * legendre_on_01(k, x_plot);
end

%% ---- 第三步：用公式计算 L2 误差（利用 Bessel 不等式，无需积分） ----
% ||f - φ*||₂² = ||f||₂² - Σ_{k=0}^n c_k² · (φ_k, φ_k)
%              = ||f||₂² - Σ_{k=0}^n c_k² / (2k+1)
%
% 这个公式说明：L2 误差随 n 增大单调减小（Bessel 不等式）

% ||f||₂² = ∫_0^1 [f(x)]² dx（需要数值计算，对于 e^x 有解析值 (e²-1)/2）
f_L2sq = integral(@(x) f(x).^2, 0, 1);   % 用数值积分（通用做法）
% 注意：如果题目给的是 e^x，解析值是 (exp(2)-1)/2，更精确

norms_sq = 1 ./ (2*(0:n)' + 1);   % 列向量 [1, 1/3, 1/5, ..., 1/(2n+1)]
% (0:n)' 是行向量 [0,1,...,n] 的转置（列向量）；./ 逐元素除

L2_err_sq = f_L2sq - sum(c.^2 .* norms_sq);
% c.^2 各系数逐元素平方；.* 对应元素相乘；sum(...) 求和
L2_err_sq = max(0, L2_err_sq);   % 防止因舍入误差变成极小负数
L2_err = sqrt(L2_err_sq);

% 最大误差
max_err = max(abs(y_exact - phi_plot));

fprintf('误差分析：\n');
fprintf('  L2  误差  ||f-φ*||₂ = %.6e\n', L2_err);
fprintf('  最大误差 ||f-φ*||∞ = %.6e\n\n', max_err);

%% ---- 第四步：绘图 ----
figure('Name', '最佳平方逼近（Legendre 正交基）', 'NumberTitle', 'off');

subplot(1, 2, 1);
plot(x_plot, y_exact, 'k-', 'LineWidth', 2, 'DisplayName', 'f(x) 精确值');
hold on;
plot(x_plot, phi_plot, 'r--', 'LineWidth', 1.5, 'DisplayName', sprintf('Legendre 逼近 n=%d', n));
hold off;
xlabel('x'); ylabel('y');
title(sprintf('n=%d 次最佳平方逼近（Legendre 正交基）', n));
legend('Location', 'northwest');
grid on;

subplot(1, 2, 2);
% semilogy：y轴用对数刻度，便于观察误差量级
% 这里改用普通坐标，直接看误差曲线
plot(x_plot, abs(y_exact - phi_plot), 'r-', 'LineWidth', 1.5);
xlabel('x'); ylabel('|f(x) - φ*(x)|');
title(sprintf('逐点误差（L2误差=%.2e）', L2_err));
grid on;

sgtitle(sprintf('Legendre 正交基：f(x) 在 [0,1] 的 %d 次最佳平方逼近', n));

%% ---- 附加：与幂函数基的对比（如果题目要求对比） ----
fprintf('========================================\n');
fprintf('与幂函数基（T01）的对比：\n');
fprintf('  Legendre 基条件数 = 1（完全稳定）\n');
fprintf('  幂函数基（Hilbert矩阵）条件数 = %.4e（极度病态）\n', cond(hilb(n+1)));
fprintf('  两者逼近精度理论上相同，但数值稳定性天壤之别\n');
fprintf('========================================\n');


%% =========================================================
%  局部函数（必须放在脚本文件末尾）
%  =========================================================

function vals = legendre_on_01(k, x)
    % 计算第 k 次 Legendre 多项式在 [0,1] 上的值
    % 即计算 P_k(t) 在 t = 2x-1 处的值（将 [0,1] 映射到 [-1,1]）
    %
    % 【三项递推公式】（避免高阶微分，数值高效）
    %   P_0(t) = 1
    %   P_1(t) = t
    %   P_{j+1}(t) = (2j+1)/(j+1)·t·P_j(t) - j/(j+1)·P_{j-1}(t)
    %
    % 输入：k - 阶数（非负整数），x - [0,1] 内的点（可为向量）
    % 输出：vals - P_k(2x-1) 在各点的值（与 x 形状相同）

    t = 2*x - 1;   % 变量替换：x∈[0,1] → t∈[-1,1]
    % 当 x 为向量时，t 也是同形状的向量（MATLAB 广播规则）

    if k == 0
        vals = ones(size(x));   % P_0(t) = 1（常数）
    elseif k == 1
        vals = t;               % P_1(t) = t = 2x-1
    else
        % 用三项递推从 P_0, P_1 推到 P_k
        p_prev = ones(size(x));   % P_0(t) = 1，与 x 等大小的向量
        p_curr = t;               % P_1(t) = t
        for j = 1 : k-1
            % 递推：P_{j+1} = [(2j+1)·t·P_j - j·P_{j-1}] / (j+1)
            p_next = ((2*j+1) * t .* p_curr - j * p_prev) / (j+1);
            % .* 是逐元素乘法（t 和 p_curr 都是向量）
            p_prev = p_curr;   % 滚动更新：前一项
            p_curr = p_next;   % 滚动更新：当前项
        end
        vals = p_curr;   % 循环结束后 p_curr 就是 P_k(t)
    end
end