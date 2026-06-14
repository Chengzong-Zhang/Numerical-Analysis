%% =========================================================
%  模板T05：复合梯形公式
%
%  【公式】
%  将 [a,b] 分成 n 等份，步长 h = (b-a)/n，节点 xᵢ = a + i·h
%
%    T_n = h/2 · [f(x₀) + f(xₙ) + 2·Σᵢ₌₁ⁿ⁻¹ f(xᵢ)]
%        = h/2 · [f(a) + f(b) + 2·(内点之和)]
%
%  【误差】
%  截断误差 = -(b-a)h²/12 · f''(ξ)，精度阶数 O(h²)
%  即步长 h 减半，误差缩小为 1/4
%
%  【最优步长】
%  h 太小时，舍入误差 O(eps/h) 主导，误差反而增大
%  梯形法最优步长：h_opt ≈ eps^(1/3) ≈ 6×10⁻⁶
%  =========================================================

clear; close all; clc;

%% ===== ① 修改这里：定义你的积分 =====
f = @(x) 4 ./ (1 + x.^2);   % 被积函数（逐元素运算！）
a = 0;                        % 积分下限
b = 1;                        % 积分上限
I_exact = pi;                 % 精确值（用于计算误差；不知道时用 integral 估算）
%% =====================================

% 如果不知道精确值，用下面这行估算（误差约 1e-12）：
% I_exact = integral(f, a, b, 'RelTol', 1e-12);

fprintf('复合梯形公式\n');
fprintf('积分区间：[%g, %g]，精确值：%.10f\n\n', a, b, I_exact);

%% ---- 用法一：给定子区间数 n，直接计算 ----
n = 1000;   % 子区间数（越大越精确，但有最优上限）
h = (b - a) / n;            % 步长
x = linspace(a, b, n+1);   % 生成 n+1 个等间距节点（含两端点）
% linspace(a, b, N) 在 [a,b] 内生成 N 个均匀点：a, a+h, a+2h, ..., b
fx = f(x);                  % 在所有节点处计算函数值（向量化）

% 套公式：端点系数1，内点系数2，再乘 h/2
T = h/2 * (fx(1) + fx(end) + 2*sum(fx(2:end-1)));
% fx(1)：第一个元素（端点 f(a)）
% fx(end)：最后一个元素（端点 f(b)）
% fx(2:end-1)：从第2个到倒数第2个（内点），共 n-1 个
% sum(v)：向量 v 所有元素之和

fprintf('n=%d（h=%.2e）：T_n = %.10f，误差 = %.4e\n', n, h, T, abs(T - I_exact));

%% ---- 用法二：测试不同 n（观察收敛阶数）----
fprintf('\n--- 改变 n，观察误差随 h 的变化 ---\n');
fprintf('%-10s  %-10s  %-20s  %-15s\n', 'n', 'h', 'T_n 结果', '误差 |T_n - I|');
fprintf('%s\n', repmat('-', 1, 58));

n_test = [10, 100, 1000, 10000, 100000];
for n_k = n_test
    h_k = (b - a) / n_k;
    x_k = linspace(a, b, n_k + 1);
    fx_k = f(x_k);
    T_k = h_k/2 * (fx_k(1) + fx_k(end) + 2*sum(fx_k(2:end-1)));
    err_k = abs(T_k - I_exact);
    fprintf('%-10d  %-10.2e  %-20.15f  %-15.4e\n', n_k, h_k, T_k, err_k);
end
fprintf('\n');
fprintf('观察：n 增大10倍（h缩小10倍），误差应缩小约100倍（因为阶数是 O(h²)）\n\n');

%% ---- 用法三：分析误差随 h 的变化（双对数图） ----
h_vals = logspace(0, -8, 80);   % h 从 1 到 1e-8，对数间隔
% logspace(a, b, n)：生成从 10^a 到 10^b 的 n 个对数均匀间隔的数

n_vals = ceil((b - a) ./ h_vals);   % 对应的 n（向上取整）
err_vals = zeros(size(h_vals));

for k = 1:length(h_vals)
    n_k = n_vals(k);
    h_k = (b - a) / n_k;
    x_k = linspace(a, b, n_k + 1);
    fx_k = f(x_k);
    T_k = h_k/2 * (fx_k(1) + fx_k(end) + 2*sum(fx_k(2:end-1)));
    err_vals(k) = abs(T_k - I_exact);
end

% 找最优步长（误差最小处）
[min_err, idx_min] = min(err_vals);
% [val, idx] = min(v)：同时返回最小值和其下标

fprintf('最优步长分析：\n');
fprintf('  最小误差 = %.4e，对应 h = %.4e，n = %d\n', ...
    min_err, (b-a)/n_vals(idx_min), n_vals(idx_min));
fprintf('  理论最优 h ≈ eps^(1/3) = %.4e\n\n', eps^(1/3));
fprintf('  结论：h < %.2e 后，舍入误差主导，再减小 h 误差反而增大\n\n', (b-a)/n_vals(idx_min));

% 绘图
figure('Name', '复合梯形公式：误差随步长变化', 'NumberTitle', 'off');
actual_h = (b - a) ./ n_vals;   % 实际步长（因为 n 取整，略有差异）
loglog(actual_h, err_vals, 'b-', 'LineWidth', 1.5, 'DisplayName', '实际误差');
% loglog：双对数坐标图（两轴都用对数刻度）
% 在双对数图上，O(h^p) 的误差曲线是斜率为 p 的直线
hold on;
% 绘制理论 O(h²) 参考线（斜率=2）
h_ref = logspace(0, -5, 50);
loglog(h_ref, 0.1*h_ref.^2, 'k--', 'LineWidth', 1, 'DisplayName', 'O(h²) 参考斜率');
% 标注最优步长
xline((b-a)/n_vals(idx_min), 'r:', '最优步长', 'LineWidth', 1.5);
% xline(x, 线型, 标签)：在 x 处画竖线
hold off;
xlabel('步长 h（对数刻度）');
ylabel('|误差|（对数刻度）');
title('复合梯形公式：误差 vs 步长（双对数图）');
legend('Location', 'NorthEast');
grid on;
axis([1e-8, 1, 1e-16, 10]);   % 设置坐标轴范围
