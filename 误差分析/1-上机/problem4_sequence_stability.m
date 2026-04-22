%% =========================================================
%  problem4_sequence_stability.m
%  功能：考察由递推公式定义的实数序列的数值稳定性
%  =========================================================
%
%  【问题描述】
%  考察实数序列 {x_n}，由以下递推关系定义：
%    x_{n+1} = (13/3) · x_n - (4/3) · x_{n-1}，  n = 1, 2, 3, ...
%
%  分别考察两组初值：
%    情形一：x_0 = 1，x_1 = 1/3
%    情形二：x_0 = 1，x_1 = 4
%
%  【数学分析：特征方程与通解】
%  这是一个二阶线性常系数递推方程，求通解的方法类似常微分方程。
%
%  设 x_n = r^n 代入递推关系，得到特征方程：
%    r^{n+1} = (13/3) r^n - (4/3) r^{n-1}
%  两边除以 r^{n-1}（r≠0）：
%    r² = (13/3) r - (4/3)
%    3r² = 13r - 4
%    3r² - 13r + 4 = 0
%
%  用求根公式：r = (13 ± √(169 - 48)) / 6 = (13 ± √121) / 6 = (13 ± 11) / 6
%    r₁ = (13 + 11)/6 = 24/6 = 4
%    r₂ = (13 - 11)/6 =  2/6 = 1/3
%
%  通解：x_n = A · 4^n + B · (1/3)^n  （A、B 由初值确定）
%
%  【情形一：x_0 = 1，x_1 = 1/3】
%  代入初值求 A、B：
%    n=0：A + B = 1
%    n=1：4A + B/3 = 1/3
%  由第二式：12A + B = 1
%  代入第一式：12A + (1-A) = 1 → 11A = 0 → A = 0，B = 1
%  精确解：x_n = (1/3)^n  （单调趋近 0）
%
%  【数值稳定性分析——情形一】
%  理想情况下 A = 0，但浮点计算中 x_0, x_1 有舍入误差，
%  实际上等价于 A ≠ 0 的微小扰动：A ≈ ε（ε 为机器精度量级）。
%  扰动后的解：x_n ≈ (1/3)^n + ε · 4^n
%  由于 4^n 指数增长，而 (1/3)^n 指数衰减，随着 n 增大：
%    误差 ε · 4^n → ∞（尽管 ε 极小），而真值 → 0
%  → 算法不稳定！（4 是增长特征根，是"寄生解"）
%
%  【情形二：x_0 = 1，x_1 = 4】
%  代入初值求 A、B：
%    n=0：A + B = 1
%    n=1：4A + B/3 = 4
%  由第二式：12A + B = 12
%  代入第一式：12A + (1-A) = 12 → 11A = 11 → A = 1，B = 0
%  精确解：x_n = 4^n  （单调增长）
%
%  【数值稳定性分析——情形二】
%  扰动后的解：x_n ≈ 4^n + ε · (1/3)^n
%  由于 (1/3)^n → 0，扰动项相对于主项 4^n 趋近 0：
%    相对误差 = ε(1/3)^n / 4^n = ε(1/12)^n → 0 迅速
%  → 算法稳定！（主项 4^n 远大于误差项，错误被"淹没"）
%  =========================================================

clear; close all; clc;

fprintf('==============================================\n');
fprintf('  题目4：递推序列的数值稳定性分析\n');
fprintf('==============================================\n\n');

% 递推系数
c1 = 13/3;   % x_n 的系数
c2 = 4/3;    % x_{n-1} 的系数（递推公式中带负号：x_{n+1} = c1*x_n - c2*x_{n-1}）

fprintf('递推关系：x_{n+1} = (13/3)·x_n - (4/3)·x_{n-1}\n');
fprintf('特征方程：3r² - 13r + 4 = 0\n');
fprintf('特征根：r₁ = 4，r₂ = 1/3\n');
fprintf('通解形式：x_n = A·4^n + B·(1/3)^n\n\n');

N = 60;   % 递推步数（足够多以看出不稳定的增长）


%% =========================================================
%  第一部分：情形一（x₀=1, x₁=1/3）——算法不稳定
%  =========================================================
%
%  精确解：x_n = (1/3)^n
%  浮点计算中存在舍入误差，等效于 A ≠ 0 的微小扰动，
%  导致 4^n 分量随 n 增大而主导结果。
%  =========================================================

fprintf('--- 情形一：x₀ = 1，x₁ = 1/3 ---\n\n');

% ---- 解析求 A, B ----
x0_case1 = 1;
x1_case1 = 1/3;
% 线性方程组：[1  1 ; 4  1/3] * [A; B] = [x0; x1]
% 用 MATLAB 的 \ 运算符（左除）求解线性方程组 Ax = b
coeff_matrix = [1, 1; 4, 1/3];  % 系数矩阵（行：n=0 和 n=1 的方程）
rhs_case1    = [x0_case1; x1_case1]; % 右端项
AB_case1 = coeff_matrix \ rhs_case1; % 求解 [A; B]
A1 = AB_case1(1);
B1 = AB_case1(2);
fprintf('情形一的系数：A = %.6e，B = %.6e\n', A1, B1);
fprintf('理论上 A = 0，B = 1，精确解 x_n = (1/3)^n\n\n');

% ---- 计算精确解（用 MATLAB 浮点运算计算 (1/3)^n，仅供参考） ----
n_axis = 0:N;
x_exact_case1 = (1/3).^n_axis;  % .^：逐元素幂运算（作用于数组中每个元素）

% ---- 正向递推（浮点计算，含舍入误差） ----
x_computed_case1 = zeros(1, N+1);  % 预分配计算结果数组
x_computed_case1(1) = x0_case1;   % x_0 = 1
x_computed_case1(2) = x1_case1;   % x_1 = 1/3

for n = 2:N
    % 递推：x_{n+1} = (13/3)·x_n - (4/3)·x_{n-1}
    % 注意 MATLAB 下标从 1 开始，所以 x_{n} 对应 x_computed(n+1)
    x_computed_case1(n+1) = c1 * x_computed_case1(n) - c2 * x_computed_case1(n-1);
end

% ---- 打印结果对比 ----
fprintf('%-5s  %-22s  %-22s  %-15s\n', 'n', '递推计算值', '精确值 (1/3)^n', '相对误差');
fprintf('%s\n', repmat('-', 1, 68));

for n = [0, 1, 2, 5, 10, 15, 20, 25, 30, 35, 40]
    if n > N; break; end
    xc = x_computed_case1(n+1);
    xe = x_exact_case1(n+1);
    if xe ~= 0
        rel_err = abs(xc - xe) / abs(xe);
    else
        rel_err = Inf;
    end
    fprintf('%4d  %-22.10e  %-22.10e  %-15.4e\n', n, xc, xe, rel_err);
end

% 找序列首次出现负值的位置（精确值始终为正）
fprintf('\n诊断：\n');
neg_idx = find(x_computed_case1 < 0, 1);  % find：返回满足条件的第一个下标
if ~isempty(neg_idx)
    fprintf('  n = %d 时首次出现负值（真值 = %.2e > 0）→ 算法已失效！\n', ...
        neg_idx - 1, x_exact_case1(neg_idx));
end

% 找序列首次超过 1 的位置（应单调减小）
big_idx = find(x_computed_case1 > 1 & n_axis > 0, 1);
if ~isempty(big_idx)
    fprintf('  n = %d 时计算值超过 1（真值 = %.2e → 0）→ 发散！\n', ...
        big_idx - 1, x_exact_case1(big_idx));
end
fprintf('\n');


%% =========================================================
%  第二部分：情形二（x₀=1, x₁=4）——算法稳定
%  =========================================================
%
%  精确解：x_n = 4^n
%  扰动来自 (1/3)^n 分量，但其相对于 4^n 趋近 0，算法稳定。
%  =========================================================

fprintf('--- 情形二：x₀ = 1，x₁ = 4 ---\n\n');

% ---- 解析求 A, B ----
x0_case2 = 1;
x1_case2 = 4;
rhs_case2 = [x0_case2; x1_case2];
AB_case2 = coeff_matrix \ rhs_case2;
A2 = AB_case2(1);
B2 = AB_case2(2);
fprintf('情形二的系数：A = %.6e，B = %.6e\n', A2, B2);
fprintf('理论上 A = 1，B = 0，精确解 x_n = 4^n\n\n');

% ---- 精确解 ----
x_exact_case2 = 4.^n_axis;  % 4^n（对数组逐元素计算）

% ---- 正向递推 ----
x_computed_case2 = zeros(1, N+1);
x_computed_case2(1) = x0_case2;
x_computed_case2(2) = x1_case2;

for n = 2:N
    x_computed_case2(n+1) = c1 * x_computed_case2(n) - c2 * x_computed_case2(n-1);
end

% ---- 打印结果对比 ----
fprintf('%-5s  %-22s  %-22s  %-15s\n', 'n', '递推计算值', '精确值 4^n', '相对误差');
fprintf('%s\n', repmat('-', 1, 68));

for n = [0, 1, 2, 5, 10, 15, 20, 25, 30]
    if n > N; break; end
    xc = x_computed_case2(n+1);
    xe = x_exact_case2(n+1);
    if xe ~= 0
        rel_err = abs(xc - xe) / abs(xe);
    else
        rel_err = 0;
    end
    fprintf('%4d  %-22.10e  %-22.10e  %-15.4e\n', n, xc, xe, rel_err);
end

fprintf('\n情形二递推结果与精确值 4^n 始终吻合，验证了算法的数值稳定性。\n\n');


%% =========================================================
%  第三部分：误差增长的理论解释与验证
%  =========================================================
%
%  情形一误差分析：
%    计算值 = (1/3)^n + δ·4^n（δ 为初始扰动，ε₀ 量级）
%    真值   = (1/3)^n
%    误差   = δ·4^n（指数增长）
%    相对误差 = δ·4^n / (1/3)^n = δ·12^n（以 12 为底指数增长！）
%
%  情形二误差分析：
%    计算值 = 4^n + δ·(1/3)^n（δ 为初始扰动）
%    真值   = 4^n
%    相对误差 = δ·(1/3)^n / 4^n = δ·(1/12)^n → 0（指数衰减！）
%  =========================================================

fprintf('--- 第三部分：误差增长的理论分析 ---\n\n');

% 估算初始扰动 δ：由于 x_1 = 1/3 = 0.333... 在双精度中有舍入误差
delta = abs(1/3 - 0.333333333333333);  % 1/3 的舍入误差（近似估算）
% 更准确地：直接用计算值与精确值的初始差
delta_real = abs(x_computed_case1(2) - (1/3));
fprintf('情形一：初始 x₁ = 1/3 的舍入误差 δ ≈ %.2e\n', delta_real);

fprintf('\n理论误差预测（情形一）：\n');
fprintf('  相对误差 ≈ δ·12^n\n');
fprintf('\n%-5s  %-20s  %-20s\n', 'n', '预测相对误差 δ·12^n', '实际相对误差');
fprintf('%s\n', repmat('-', 1, 48));

% 取初始误差的精确估计（包含 x_0 和 x_1 两步的累积）
% 实际上直接比较 computed vs exact 即可
for n = [5, 10, 15, 20, 25, 30]
    if n > N; break; end
    xe = x_exact_case1(n+1);
    xc = x_computed_case1(n+1);
    if xe ~= 0
        actual_rel = abs(xc - xe) / abs(xe);
    else
        actual_rel = Inf;
    end
    predicted_rel = eps * 12^n;  % 用机器精度 eps 作为初始扰动估计
    fprintf('%4d  %-20.4e  %-20.4e\n', n, predicted_rel, actual_rel);
end

fprintf('\n（预测值仅为量级估计，误差与实际初始扰动的精确值有关）\n\n');


%% =========================================================
%  第四部分：图形可视化
%  =========================================================

fprintf('--- 第四部分：图形可视化 ---\n\n');

% 限制显示范围（情形一超出范围后数值极大，取有意义的前几十步）
n_show = min(50, N);  % 显示前 n_show+1 个点

figure('Name', '题目4：递推序列稳定性', 'NumberTitle', 'off', 'Position', [100 100 1050 600]);

% ---- 左上：情形一——计算值 vs 精确值（绝对值对数坐标） ----
subplot(2, 2, 1);

n_plot = 0:n_show;
% 取绝对值后用对数坐标（因为正向递推可能产生负值）
valid1 = x_computed_case1(1:n_show+1) ~= 0;  % 避免 log(0)

semilogy(n_plot, abs(x_exact_case1(1:n_show+1)),    'g-o', ...
    'MarkerSize', 4, 'LineWidth', 1.5, 'DisplayName', '精确值 (1/3)^n');
hold on;
semilogy(n_plot, abs(x_computed_case1(1:n_show+1)), 'r--s', ...
    'MarkerSize', 4, 'LineWidth', 1.5, 'DisplayName', '递推计算值');

xlabel('n');
ylabel('|x_n|（对数刻度）');
title('情形一：x₀=1, x₁=1/3（不稳定）');
legend('Location', 'NorthEast');
grid on;

% ---- 右上：情形二——计算值 vs 精确值 ----
subplot(2, 2, 2);

semilogy(n_plot, x_exact_case2(1:n_show+1),    'g-o', ...
    'MarkerSize', 4, 'LineWidth', 1.5, 'DisplayName', '精确值 4^n');
hold on;
semilogy(n_plot, abs(x_computed_case2(1:n_show+1)), 'b--^', ...
    'MarkerSize', 4, 'LineWidth', 1.5, 'DisplayName', '递推计算值');

xlabel('n');
ylabel('|x_n|（对数刻度）');
title('情形二：x₀=1, x₁=4（稳定）');
legend('Location', 'NorthWest');
grid on;

% ---- 左下：情形一的相对误差增长 ----
subplot(2, 2, 3);

rel_err_case1 = zeros(1, n_show+1);
for i = 1:n_show+1
    xe = x_exact_case1(i);
    xc = x_computed_case1(i);
    if xe ~= 0
        rel_err_case1(i) = abs(xc - xe) / abs(xe);
    else
        rel_err_case1(i) = NaN;  % NaN：Not a Number，跳过该点
    end
end

% 绘制实际相对误差
semilogy(n_plot, rel_err_case1, 'r-s', 'MarkerSize', 4, 'LineWidth', 1.5, ...
    'DisplayName', '实际相对误差');
hold on;

% 绘制理论预测：eps × 12^n
theory_rel_err = eps * 12.^n_plot;
semilogy(n_plot, theory_rel_err, 'k--', 'LineWidth', 1.5, ...
    'DisplayName', '理论预测 eps×12^n');

% 机器精度参考线
yline(eps, 'b:', 'LineWidth', 1.2, 'DisplayName', 'eps');

xlabel('n');
ylabel('相对误差（对数刻度）');
title('情形一相对误差：以 12^n 速率增长');
legend('Location', 'NorthWest');
grid on;

% ---- 右下：情形二的相对误差（应趋于 0）----
subplot(2, 2, 4);

rel_err_case2 = zeros(1, n_show+1);
for i = 1:n_show+1
    xe = x_exact_case2(i);
    xc = x_computed_case2(i);
    if xe ~= 0
        rel_err_case2(i) = abs(xc - xe) / abs(xe);
    end
end
% 将恰好为 0 的误差替换为 eps（避免 log(0)）
rel_err_case2(rel_err_case2 == 0) = eps * 0.1;

semilogy(n_plot, rel_err_case2, 'b-^', 'MarkerSize', 4, 'LineWidth', 1.5, ...
    'DisplayName', '实际相对误差');
hold on;

% 理论预测：相对误差 ≈ eps × (1/12)^n（衰减）
theory_rel_err2 = eps * (1/12).^n_plot;
semilogy(n_plot, theory_rel_err2, 'k--', 'LineWidth', 1.5, ...
    'DisplayName', '理论预测 eps×(1/12)^n');

yline(eps, 'r:', 'LineWidth', 1.2, 'DisplayName', 'eps');

xlabel('n');
ylabel('相对误差（对数刻度）');
title('情形二相对误差：迅速衰减（稳定）');
legend('Location', 'NorthEast');
grid on;

sgtitle('递推序列 x_{n+1}=(13/3)x_n-(4/3)x_{n-1} 的数值稳定性分析');

fprintf('图形已生成，请查看弹出的图形窗口。\n');


%% =========================================================
%  结论总结
%  =========================================================
fprintf('\n==============================================\n');
fprintf('               结论总结\n');
fprintf('==============================================\n\n');

fprintf('递推关系 x_{n+1} = (13/3)x_n - (4/3)x_{n-1} 的特征根为 r₁=4, r₂=1/3。\n\n');

fprintf('情形一（x₀=1, x₁=1/3）：\n');
fprintf('  - 精确解 x_n = (1/3)^n（r₂ 分量，趋近 0）\n');
fprintf('  - 浮点误差激发 r₁=4 对应的分量 δ·4^n（"寄生解"）\n');
fprintf('  - 相对误差以 12^n 速率增长 → 算法不稳定！\n\n');

fprintf('情形二（x₀=1, x₁=4）：\n');
fprintf('  - 精确解 x_n = 4^n（r₁ 分量，快速增长）\n');
fprintf('  - 浮点误差激发 r₂=1/3 对应的分量 δ·(1/3)^n（快速衰减）\n');
fprintf('  - 相对误差以 (1/12)^n 速率衰减 → 算法稳定！\n\n');

fprintf('一般规律：\n');
fprintf('  若算法的"主解"对应较小的特征根（|r₂|=1/3），\n');
fprintf('  而误差会激发较大特征根（|r₁|=4）的增长，则算法不稳定。\n');
fprintf('  若算法的"主解"对应较大的特征根，误差的小特征根分量\n');
fprintf('  相对主项趋近 0，则算法稳定。\n\n');

fprintf('程序运行完毕。\n');