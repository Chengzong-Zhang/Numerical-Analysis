%% =========================================================
%  problem3_integral_recursion.m
%  功能：验证用递推公式计算积分 y_n = ∫₀¹ x^n eˣ dx 的算法不稳定性
%  =========================================================
%
%  【积分的定义与递推关系】
%  定义积分序列：
%    y_n = ∫₀¹ x^n eˣ dx，  n = 0, 1, 2, ...
%
%  初始值（n=0 时的精确值）：
%    y_0 = ∫₀¹ eˣ dx = [eˣ]₀¹ = e - 1 ≈ 1.71828...
%
%  由分部积分法推导递推公式：
%    ∫₀¹ x^n eˣ dx = [x^n eˣ]₀¹ - n∫₀¹ x^{n-1} eˣ dx
%                  = e - n·y_{n-1}
%  即：y_n = e - n·y_{n-1}  （n = 1, 2, 3, ...）
%  改写索引：y_{n+1} = e - (n+1)·y_n  （即题目中的形式）
%
%  【积分的真实值界估计】
%  由于 x ∈ [0,1] 时 1 ≤ eˣ ≤ e，有：
%    1/(n+1) ≤ y_n ≤ e/(n+1)
%  即 y_n > 0 且随 n 增大而趋近 0。
%
%  【为什么正向递推不稳定？——误差分析】
%  设 ε_n 为第 n 步的累积误差（ỹ_n - y_n，ỹ_n 是计算值）。
%  由递推公式 y_{n+1} = e - (n+1)·y_n，误差满足：
%    ε_{n+1} = -(n+1)·ε_n
%  因此：
%    |ε_n| = n! · |ε_0|
%  即使初始误差 ε_0 极小（如机器精度量级），经过 n 步后误差被放大 n! 倍！
%  对于 n = 20，n! ≈ 2.4 × 10¹⁸，相当于将 eps ≈ 2.2e-16 放大到 O(1)——完全失真。
%
%  【稳定的改进：反向递推】
%  由 y_{n+1} = e - (n+1)·y_n 解出 y_n：
%    y_n = (e - y_{n+1}) / (n+1)
%  误差传播：ε_n = -ε_{n+1} / (n+1)
%  即每步误差缩小 (n+1) 倍，从大 n 往小 n 递推时误差急剧衰减 → 数值稳定！
%  起点：取足够大的 N（如 N=40），令 ỹ_N = 0（或 e/(N+1) 的近似），
%  误差 ε_N 最多为 e/(N+1) ≈ 0.068（N=40 时），
%  反向传播后 ε_0 ≈ ε_N / N! ≈ 0，对低阶项精度无影响。
%  =========================================================

clear; close all; clc;

fprintf('==============================================\n');
fprintf('  题目3：积分递推公式的数值稳定性验证\n');
fprintf('==============================================\n\n');

% 自然常数 e（MATLAB 内置：exp(1)）
e_val = exp(1);  % e = 2.71828...
fprintf('使用自然常数 e = exp(1) = %.15f\n\n', e_val);

% 计算阶数范围（验证 n = 0 到 N_forward）
N_forward = 25;   % 正向递推到第 N_forward 阶
N_backward = 40;  % 反向递推从第 N_backward 阶开始（足够大，使起点误差可忽略）


%% =========================================================
%  第一部分：用 MATLAB 数值积分计算参考真值
%  =========================================================
%
%  integral(fun, a, b)：用自适应数值积分计算 ∫_a^b fun(x) dx
%  此函数基于高精度算法（如 Gauss-Kronrod 公式），结果可视为"精确值"参考。
%  =========================================================

fprintf('--- 第一部分：参考真值（MATLAB 数值积分）---\n\n');

y_true = zeros(1, N_forward + 1);  % 预分配：y_true(k) 存储 y_{k-1} 的真值（MATLAB 下标从1起）

for n = 0:N_forward
    % 定义被积函数：x^n * e^x
    % @(x) 是匿名函数（anonymous function）语法，定义一个临时函数
    % x.^n 用逐元素幂运算（.^），使函数接受向量输入（integral 要求）
    integrand = @(x) x.^n .* exp(x);

    % integral：MATLAB 内置自适应数值积分函数
    % 'RelTol', 1e-15：相对误差容限（要求极高精度）
    % 'AbsTol', 1e-20：绝对误差容限
    y_true(n+1) = integral(integrand, 0, 1, 'RelTol', 1e-15, 'AbsTol', 1e-20);
end

fprintf('参考真值（前10阶）：\n');
for n = 0:9
    fprintf('  y_%2d = %.15f\n', n, y_true(n+1));
end
fprintf('  ...\n\n');

% 验证：真值应满足 y_n ∈ [1/(n+1), e/(n+1)]
fprintf('验证真值界：y_n ∈ [1/(n+1), e/(n+1)]\n');
for n = [5, 10, 15, 20]
    lb = 1/(n+1);   % 下界：eˣ ≥ 1
    ub = e_val/(n+1); % 上界：eˣ ≤ e
    in_range = (y_true(n+1) >= lb) && (y_true(n+1) <= ub);
    fprintf('  y_%2d = %.6e  ∈ [%.4e, %.4e]? %s\n', ...
        n, y_true(n+1), lb, ub, string(in_range));
end
fprintf('\n');


%% =========================================================
%  第二部分：正向递推（Forward Recursion）—— 验证不稳定性
%  =========================================================
%
%  递推公式：y_{n+1} = e - (n+1) · y_n
%  初始值：y_0 = e - 1（精确值）
%
%  注意：y_0 = e - 1 本身有舍入误差 ε_0 ≈ eps（机器精度量级），
%  这个微小误差在后续步骤被 (n+1)! 放大，导致结果完全错误。
%  =========================================================

fprintf('--- 第二部分：正向递推（不稳定）---\n\n');

y_forward = zeros(1, N_forward + 1);  % 正向递推计算值

% 初始值：y_0 = e - 1（精确，但存在 O(eps) 的舍入误差）
y_forward(1) = e_val - 1;

% 正向递推：y_{n+1} = e - (n+1) * y_n
for n = 0 : N_forward - 1
    % n+2 是 MATLAB 下标，对应 y_{n+1}
    y_forward(n+2) = e_val - (n+1) * y_forward(n+1);
    % 每步误差放大因子 = n+1，即 ε_{n+1} = (n+1) · ε_n（绝对值）
end

% 显示正向递推结果与真值的比较
fprintf('%-5s  %-20s  %-20s  %-15s  %-15s\n', ...
    'n', '正向递推计算值', '参考真值', '绝对误差', '相对误差');
fprintf('%s\n', repmat('-', 1, 78));

for n = 0:N_forward
    y_f = y_forward(n+1);
    y_r = y_true(n+1);
    abs_err = abs(y_f - y_r);
    if y_r ~= 0
        rel_err = abs_err / abs(y_r);
    else
        rel_err = Inf;
    end
    % 只在关键节点打印（前5阶 + 5的倍数），否则表格太长
    if n <= 4 || mod(n, 5) == 0
        fprintf('%4d  %-20.10e  %-20.10e  %-15.4e  %-15.4e\n', ...
            n, y_f, y_r, abs_err, rel_err);
    end
end

% 找正向递推第一次出现负值或绝对误差超过真值的位置
fprintf('\n正向递推发散分析：\n');
for n = 0:N_forward
    if y_forward(n+1) < 0
        fprintf('  n = %d 时首次出现负值（真值为正），算法已完全失效！\n', n);
        break;
    end
end

% 理论预测：误差放大为 n! * eps
fprintf('\n理论误差预测（n! × eps × y_0）：\n');
factorial_n = 1;
for n = 1:min(20, N_forward)
    factorial_n = factorial_n * n;
    predicted = factorial_n * eps * abs(y_forward(1));
    fprintf('  n=%2d: 预测 |ε_n| ≈ %.2e,  实际 |ε_n| = %.2e\n', ...
        n, predicted, abs(y_forward(n+1) - y_true(n+1)));
    if predicted > 10 && n > 5
        fprintf('  （误差已超过 10，后续不再显示）\n');
        break;
    end
end
fprintf('\n');


%% =========================================================
%  第三部分：反向递推（Backward Recursion）—— 验证稳定性
%  =========================================================
%
%  递推公式：y_n = (e - y_{n+1}) / (n+1)
%  起点：令 ỹ_{N_backward} = 0（近似，真值约为 e/(N+1) ≈ 0.067）
%
%  误差传播：ε_n = -ε_{n+1} / (n+1)
%  反向传播时误差每步缩小，从 n=N_backward 到 n=0，误差缩小 N_backward! 倍。
%  因此即使起点误差很大，最终 y_0 的误差也极小。
%  =========================================================

fprintf('--- 第三部分：反向递推（稳定）---\n\n');

y_backward = zeros(1, N_backward + 1);  % 预分配反向递推数组

% 起点：ỹ_{N_backward} = 0（有意取简单近似，误差约 e/(N_backward+1)）
y_backward(N_backward + 1) = 0;
starting_error = y_true(min(N_backward+1, length(y_true)));
fprintf('反向起点：ỹ_%d = 0（真值约 e/%d ≈ %.4e，起点误差：%.4e）\n\n', ...
    N_backward, N_backward+1, e_val/(N_backward+1), e_val/(N_backward+1));

% 反向递推：y_n = (e - y_{n+1}) / (n+1)
for n = N_backward - 1 : -1 : 0
    % n+1 是 MATLAB 下标，对应 y_n
    % 由 y_{n+1} = e - (n+1)*y_n 解出 y_n：
    y_backward(n+1) = (e_val - y_backward(n+2)) / (n+1);
end

% 显示反向递推结果与真值的比较
fprintf('%-5s  %-20s  %-20s  %-15s\n', 'n', '反向递推计算值', '参考真值', '相对误差');
fprintf('%s\n', repmat('-', 1, 63));

for n = 0:min(20, N_forward)
    y_b = y_backward(n+1);
    y_r = y_true(n+1);
    if y_r ~= 0
        rel_err = abs(y_b - y_r) / abs(y_r);
    else
        rel_err = 0;
    end
    if n <= 4 || mod(n, 5) == 0
        fprintf('%4d  %-20.15f  %-20.15f  %-15.4e\n', n, y_b, y_r, rel_err);
    end
end
fprintf('\n反向递推结果与参考真值高度吻合，验证了反向递推的数值稳定性。\n\n');


%% =========================================================
%  第四部分：图形可视化
%  =========================================================

fprintf('--- 第四部分：图形可视化 ---\n\n');

n_axis = 0:N_forward;  % x 轴：阶数 n

figure('Name', '题目3：积分递推公式稳定性', 'NumberTitle', 'off', 'Position', [100 100 1000 500]);

% ---- 上图：三种方法的计算值对比 ----
subplot(2, 1, 1);

% 参考真值（取对数绝对值，便于在大范围内比较）
semilogy(n_axis, abs(y_true(1:N_forward+1)), 'g-o', ...
    'MarkerSize', 5, 'LineWidth', 1.5, 'DisplayName', '参考真值（数值积分）');
hold on;

% 正向递推值（绝对值，因为可能为负）
semilogy(n_axis, abs(y_forward), 'r--s', ...
    'MarkerSize', 5, 'LineWidth', 1.5, 'DisplayName', '正向递推（不稳定）');

% 反向递推值（仅显示前 N_forward+1 项）
semilogy(n_axis, abs(y_backward(1:N_forward+1)), 'b-.^', ...
    'MarkerSize', 5, 'LineWidth', 1.5, 'DisplayName', '反向递推（稳定）');

xlabel('n（积分阶数）');
ylabel('|y_n|（对数刻度）');
title('三种方法计算 y_n = ∫₀¹ xⁿeˣdx 的结果');
legend('Location', 'NorthEast');
grid on;

% ---- 下图：正向递推的误差增长 ----
subplot(2, 1, 2);

% 计算正向递推的绝对误差
abs_err_forward = abs(y_forward - y_true(1:N_forward+1));

% 避免 log(0) 问题：将 0 替换为极小值
abs_err_forward(abs_err_forward == 0) = eps;

semilogy(n_axis, abs_err_forward, 'r-s', ...
    'MarkerSize', 5, 'LineWidth', 1.5, 'DisplayName', '正向递推绝对误差');
hold on;

% 理论误差曲线：|ε_n| ≈ n! × eps
theory_err = zeros(1, N_forward+1);
theory_err(1) = eps * abs(y_forward(1));
for n = 1:N_forward
    theory_err(n+1) = theory_err(n) * n;  % 乘以 n（对应 |ε_n| = n! × ε_0）
end
semilogy(n_axis, theory_err, 'k--', 'LineWidth', 1.5, 'DisplayName', '理论误差 n!×ε₀');

% 机器精度参考线
yline(eps, 'b:', 'LineWidth', 1.2, 'DisplayName', 'eps（机器精度）');
% yline：在指定 y 值处画水平线

xlabel('n（积分阶数）');
ylabel('正向递推绝对误差（对数刻度）');
title('正向递推误差按 n! 增长（验证不稳定性）');
legend('Location', 'NorthWest');
grid on;

sgtitle('积分 y_n = ∫₀¹ xⁿeˣdx 递推公式的稳定性分析');

fprintf('图形已生成，请查看弹出的图形窗口。\n');
fprintf('\n程序运行完毕。\n');