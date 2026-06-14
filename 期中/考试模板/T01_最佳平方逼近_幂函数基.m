%% =========================================================
%  模板T01：最佳平方逼近 —— 幂函数基（Hilbert矩阵）
%
%  【问题描述】
%  求函数 f(x) 在区间 [a,b] 上的 n 次最佳平方逼近多项式：
%    φ*(x) = c₀ + c₁x + c₂x² + ... + cₙxⁿ
%  使得 ∫_a^b [f(x) - φ*(x)]² dx 最小
%
%  【核心公式】
%  法方程组：G·c = f_vec
%    G(i,j) = (φᵢ, φⱼ) = ∫_a^b xⁱ · xʲ dx = (b^{i+j+1} - a^{i+j+1}) / (i+j+1)
%    特别地，当 [a,b] = [0,1] 时，G 就是 Hilbert 矩阵，G(i,j) = 1/(i+j-1)
%    f_vec(j) = (f, φⱼ) = ∫_a^b f(x) · xʲ dx
%
%  【注意】幂函数基的 Gram 矩阵是 Hilbert 矩阵，极度病态，
%          条件数随 n 指数增长（n=5: ~10^12，n=10: ~10^13）
%  =========================================================

clear; close all; clc;

%% ===== ① 修改这里：定义你的问题 =====
f = @(x) exp(x);    % 被逼近函数，改成题目要求的函数
a = 0;               % 积分区间左端点
b = 1;               % 积分区间右端点
n = 5;               % 逼近多项式次数（修改这里）
%% =====================================

fprintf('最佳平方逼近（幂函数基）：n=%d 次多项式，区间[%g,%g]\n\n', n, a, b);

%% ---- 第一步：构造 Gram 矩阵 G（Hilbert 矩阵型） ----
% G 是 (n+1)×(n+1) 的对称正定矩阵
% G(i,j) = ∫_a^b x^(i-1) · x^(j-1) dx = ∫_a^b x^(i+j-2) dx
%         = [x^(i+j-1)/(i+j-1)]_a^b = (b^(i+j-1) - a^(i+j-1)) / (i+j-1)
%
% 当 [a,b]=[0,1] 时：G(i,j) = 1/(i+j-1)，正好是 Hilbert 矩阵，用 hilb(n+1)
% 当 [a,b] 不是 [0,1] 时：需要手动计算每个元素

if a == 0 && b == 1
    % 特殊情形：[0,1] 区间，直接用 MATLAB 内置 Hilbert 矩阵
    G = hilb(n+1);
    % hilb(m) 生成 m×m Hilbert 矩阵，H(i,j) = 1/(i+j-1)
else
    % 一般区间 [a,b]：手动计算
    G = zeros(n+1, n+1);   % 预分配 (n+1)×(n+1) 矩阵
    for i = 0:n
        for j = 0:n
            % G(i+1, j+1) = ∫_a^b x^i · x^j dx = ∫_a^b x^(i+j) dx
            %             = (b^(i+j+1) - a^(i+j+1)) / (i+j+1)
            G(i+1, j+1) = (b^(i+j+1) - a^(i+j+1)) / (i+j+1);
        end
    end
end

% 输出条件数（衡量病态程度）
kappa = cond(G);
% cond(A) 返回矩阵 A 的 2-范数条件数 = σ_max / σ_min
% 含义：求解 Gc=f_vec 时，系数误差的放大倍数上限
% 如果 κ ≈ 10^k，则求解过程约损失 k 位有效数字
fprintf('Gram 矩阵（%d×%d）条件数：cond(G) = %.4e\n', n+1, n+1, kappa);
fprintf('（损失约 %.0f 位有效数字，机器精度 eps≈2.2e-16）\n\n', log10(kappa));
% log10(x) 以10为底的对数，用来估计损失几位十进制有效数字

%% ---- 第二步：构造右端向量 f_vec ----
% f_vec(j+1) = (f, φⱼ) = ∫_a^b f(x) · x^j dx，j = 0, 1, ..., n
f_vec = zeros(n+1, 1);   % 预分配列向量，避免循环中动态扩容
for j = 0:n
    % integral(fun, a, b)：自适应数值积分（内部用 Gauss-Kronrod 规则）
    % @(x) f(x) .* x.^j：匿名函数，x.^j 是逐元素幂（支持向量输入）
    f_vec(j+1) = integral(@(x) f(x) .* x.^j, a, b);
end

%% ---- 第三步：求解法方程组 G·c = f_vec ----
% A \ b（MATLAB 左除运算符）：求解线性方程组 Ax = b
% 比 inv(G)*f_vec 更稳定（内部用 Cholesky 或 LU 分解）
% c(k+1) 是 x^k 的系数，即 φ*(x) = c(1) + c(2)*x + c(3)*x² + ...
c = G \ f_vec;

fprintf('逼近系数（从低次到高次）：\n');
for k = 0:n
    fprintf('  c_%d（x^%d 的系数）= %12.8f\n', k, k, c(k+1));
end
fprintf('\n');

%% ---- 第四步：在密集点上求逼近多项式的值（Horner 嵌套法） ----
% 生成绘图用的密集点（不参与计算，只用于可视化）
x_plot = linspace(a, b, 500);   % 500 个等间距点，足够光滑

% Horner 嵌套乘法（秦九韶算法）求多项式值
% φ*(x) = c(1) + c(2)x + c(3)x² + ... + c(n+1)x^n
%       = (...((c(n+1)·x + c(n))·x + c(n-1))·x + ... + c(1))
% 只需 n 次乘法和 n 次加法，比直接求幂更快、更稳定
phi_plot = c(n+1) * ones(size(x_plot));  % 初始化为最高次系数
% ones(size(x_plot))：生成与 x_plot 等大小的全1数组，用于标量×向量
for k = n:-1:1
    % n:-1:1 生成从 n 递减到 1 的整数序列（步长 -1）
    phi_plot = phi_plot .* x_plot + c(k);
    % .* 是逐元素乘法（x_plot 是向量，每个元素分别乘）
    % + c(k) 是标量加到向量每个元素（MATLAB 自动广播）
end

%% ---- 第五步：计算误差 ----
y_exact = f(x_plot);   % 精确函数值

% L2 误差（均方误差）：||f - φ*||₂ = sqrt(∫_a^b [f(x)-φ*(x)]² dx)
% 用数值积分计算，而不是在绘图点上近似（更准确）
L2_err = sqrt(integral(@(x) (f(x) - horner_eval(c, x)).^2, a, b));
% horner_eval：在任意点集上用 Horner 法求多项式值，定义在文件末尾

% 最大误差（无穷范数）：||f - φ*||∞ = max|f(x) - φ*(x)|
max_err = max(abs(y_exact - phi_plot));
% abs(v) 向量逐元素取绝对值；max(v) 向量最大元素

fprintf('误差分析：\n');
fprintf('  L2  误差  ||f-φ*||₂ = %.6e\n', L2_err);
fprintf('  最大误差 ||f-φ*||∞ = %.6e\n\n', max_err);

%% ---- 第六步：绘图 ----
figure('Name', '最佳平方逼近（幂函数基）', 'NumberTitle', 'off');
% figure 创建新图形窗口；'Name' 设置窗口标题

subplot(1, 2, 1);   % 1行2列的第1个子图
plot(x_plot, y_exact, 'k-', 'LineWidth', 2, 'DisplayName', 'f(x) 精确值');
% 'k-'：黑色(k)实线(-); 'LineWidth'：线宽; 'DisplayName'：图例文字
hold on;   % 保持当前图，后续绘图叠加
plot(x_plot, phi_plot, 'b--', 'LineWidth', 1.5, 'DisplayName', sprintf('φ*(x) n=%d 逼近', n));
hold off;
xlabel('x'); ylabel('y');
title(sprintf('n=%d 次最佳平方逼近（幂函数基）\ncond(G)=%.2e', n, kappa));
legend('Location', 'northwest');   % 'Location'：图例位置
grid on;   % 显示网格线

subplot(1, 2, 2);   % 第2个子图：误差曲线
plot(x_plot, abs(y_exact - phi_plot), 'r-', 'LineWidth', 1.5);
xlabel('x'); ylabel('|f(x) - φ*(x)|');
title(sprintf('逐点误差（最大误差 = %.2e）', max_err));
grid on;

sgtitle(sprintf('最佳平方逼近（幂函数基）：f(x) 在 [%g,%g] 上的 %d 次逼近', a, b, n));
% sgtitle：整个 figure 的总标题（R2018b 以上版本支持）


%% =========================================================
%  局部函数（必须放在脚本文件末尾）
%  =========================================================

function vals = horner_eval(c, x)
    % 用 Horner 嵌套法（秦九韶算法）在点集 x 上求多项式值
    % 输入：c - 系数向量（c(k+1) 是 x^k 的系数），x - 求值点（可为向量）
    % 输出：vals - 多项式在各点的值
    n = length(c) - 1;
    vals = c(n+1) * ones(size(x));   % 从最高次系数初始化
    for k = n:-1:1
        vals = vals .* x + c(k);    % 逐元素：每个 x 分量分别计算
    end
end
