%% =========================================================
%  模板T11：Monte Carlo 数值积分 —— 两种方法
%
%  【问题描述】
%  用随机采样估计定积分 I = ∫_a^b f(x) dx
%
%  【方法一：随机投点法（Hit-or-Miss）】
%  在矩形区域 [a,b]×[0, f_max] 内均匀撒点，
%  统计落在曲线 f(x) 下方的比例：
%    I ≈ (b-a)·f_max · (命中次数/总点数)
%  方差：σ² = (b-a)²·f_max² · P(1-P) / N，其中 P=I/((b-a)·f_max)
%
%  【方法二：样本均值法（Sample Mean，更高效！）】
%  在 [a,b] 均匀采样 N 个点 x₁,...,xN，对函数值求平均：
%    I = ∫_a^b f(x) dx = (b-a)·E[f(X)] ≈ (b-a)·(1/N)·Σf(xᵢ)
%  方差：σ² = (b-a)²·Var[f(X)] / N
%
%  【对比】
%  样本均值法方差更小（通常小6倍以上），用相同随机数得到更精确的结果。
%  推荐优先使用样本均值法！
%
%  【收敛速度】
%  两种方法误差都是 O(1/√N)，与维数无关（高维积分的优势）。
%  =========================================================

clear; close all; clc;

%% ===== ① 修改这里：定义你的积分 =====
f = @(x) exp(-x);            % 被积函数
a = 0;                        % 积分下限
b = 1;                        % 积分上限
I_exact = 1 - exp(-1);       % 精确值（用 integral 估算也可以）
N = 1e6;                      % 随机点数（越大越精确）
%% =====================================

% 如果不知道精确值：
% I_exact = integral(f, a, b, 'RelTol', 1e-12);

fprintf('Monte Carlo 数值积分\n');
fprintf('∫_%.0f^%.0f f(x) dx，精确值 = %.10f\n\n', a, b, I_exact);

%% ===================================================
%  方法一：随机投点法（Hit-or-Miss）
%  ===================================================
fprintf('--- 方法一：随机投点法（Hit-or-Miss）---\n');

% 需要 f(x) 在 [a,b] 上的最大值（用于确定矩形上界）
% 简单估算：在密集点上取最大值
x_check = linspace(a, b, 10000);
f_max = max(f(x_check)) * 1.05;   % 多留5%的余量
% 如果能解析确定最大值更好（如 f(x)=e^{-x} 在 [0,1] 上最大值是 f(0)=1）

fprintf('f_max = %.4f（矩形高度）\n', f_max);
fprintf('矩形面积 = (b-a)·f_max = %.4f\n\n', (b-a)*f_max);

% 生成 N 个均匀随机点（向量化，一次性生成所有）
x_rand = a + (b-a) * rand(1, N);   % x 坐标：均匀分布在 [a,b]
y_rand = f_max * rand(1, N);       % y 坐标：均匀分布在 [0, f_max]
% rand(1, N)：生成 N 个 [0,1] 均匀随机数

% 判断每个点是否落在曲线下方（即 y ≤ f(x)）
hit = (y_rand <= f(x_rand));   % 逻辑向量（true=命中，false=未命中）
% f(x_rand)：对向量 x_rand 每个元素求 f 值（向量化）
% <= 是逐元素比较，返回同大小的逻辑向量

N_hit = sum(hit);         % 命中次数（sum 对逻辑向量求和）
P_hat = N_hit / N;        % 命中概率的估计值
I_hit = (b-a) * f_max * P_hat;   % 积分估计

fprintf('命中次数 %d / %d（命中率 %.4f）\n', N_hit, N, P_hat);
fprintf('结果 = %.10f，误差 = %.4e\n', I_hit, abs(I_hit - I_exact));
fprintf('理论标准差 = %.4e（≈误差量级）\n\n', ...
    (b-a)*f_max*sqrt(P_hat*(1-P_hat)/N));

%% ===================================================
%  方法二：样本均值法（Sample Mean）
%  ===================================================
fprintf('--- 方法二：样本均值法（Sample Mean）---\n');

% 在 [a,b] 均匀采样 N 个点
x_mean = a + (b-a) * rand(1, N);   % N 个均匀随机点

% 对函数值求平均，再乘以区间长度
% I = ∫_a^b f(x) dx = (b-a)·E[f(X)]，X~Uniform(a,b)
f_vals = f(x_mean);                     % 在 N 个随机点处求 f 值（向量化）
I_mean = (b-a) * mean(f_vals);         % 平均值乘区间长度
% mean(v)：向量 v 的算术平均值 = (1/N)·Σv_i

f_std = std(f_vals);   % f 在采样点处的样本标准差
% std(v)：样本标准差（分母为 N-1）

fprintf('f 均值 = %.6f，f 标准差 = %.6f\n', mean(f_vals), f_std);
fprintf('结果 = %.10f，误差 = %.4e\n', I_mean, abs(I_mean - I_exact));
fprintf('理论标准差 = %.4e（= (b-a)·σ_f/√N）\n\n', (b-a)*f_std/sqrt(N));

%% ---- 效率对比 ----
fprintf('--- 效率对比 ---\n');
var_hit  = (b-a)^2 * f_max^2 * P_hat*(1-P_hat);   % 投点法方差（单次估计）
var_mean = (b-a)^2 * f_std^2;                       % 均值法方差（单次估计）
fprintf('投点法方差 = %.6e\n', var_hit/N);
fprintf('均值法方差 = %.6e\n', var_mean/N);
fprintf('效率比 = %.2f（均值法方差是投点法的 %.2f 倍，越小越好）\n\n', ...
    var_hit/var_mean, var_mean/var_hit);
fprintf('结论：样本均值法更高效，方差更小！\n\n');

%% ---- 收敛性分析（误差随 N 的变化）----
fprintf('--- 收敛性分析 ---\n');
N_range = round(logspace(2, log10(N), 50));   % 从 100 到 N 的 50 个检查点
% logspace(a, b, n)：对数均匀间隔的 n 个点（10^a 到 10^b）

err_hit_arr  = zeros(size(N_range));   % 投点法误差
err_mean_arr = zeros(size(N_range));   % 均值法误差

% 对同一组随机数累积计算（用 cumsum 和 cumulative mean）
f_vals_all = f(a + (b-a)*rand(1, N));     % 均值法用的函数值
y_rand_full = f_max * rand(1, N);         % 投点法用的 y 坐标
x_rand_full = a + (b-a) * rand(1, N);     % 投点法用的 x 坐标
hit_full = (y_rand_full <= f(x_rand_full));

cum_f   = cumsum(f_vals_all);   % 函数值累积和
cum_hit = cumsum(hit_full);     % 命中次数累积和

for idx = 1:length(N_range)
    n_k = N_range(idx);
    I_mean_k = (b-a) * cum_f(n_k) / n_k;
    I_hit_k  = (b-a) * f_max * cum_hit(n_k) / n_k;
    err_hit_arr(idx)  = abs(I_hit_k  - I_exact);
    err_mean_arr(idx) = abs(I_mean_k - I_exact);
end

%% ---- 绘图 ----
figure('Name', 'Monte Carlo 积分', 'NumberTitle', 'off', 'Position', [50,50,1200,800]);

% 子图1：积分示意（随机投点）
subplot(2, 3, 1);
N_show = min(2000, N);
x_s = x_rand(1:N_show); y_s = y_rand(1:N_show);
hit_s = hit(1:N_show);
% 绘制命中点（绿色）和未命中点（红色）
plot(x_s(hit_s), y_s(hit_s), 'g.', 'MarkerSize', 2, 'DisplayName', '命中');
hold on;
plot(x_s(~hit_s), y_s(~hit_s), 'r.', 'MarkerSize', 2, 'DisplayName', '未命中');
% 绘制函数曲线
x_curve = linspace(a, b, 300);
plot(x_curve, f(x_curve), 'b-', 'LineWidth', 2, 'DisplayName', 'f(x)');
hold off;
xlabel('x'); ylabel('y');
title(sprintf('随机投点法（前%d个点）', N_show));
legend('Location', 'northeast', 'FontSize', 7);
grid on;

% 子图2：函数值分布（均值法）
subplot(2, 3, 2);
histogram(f_vals(1:min(N, 50000)), 50, 'Normalization', 'pdf', 'FaceColor', [0.3 0.6 0.9]);
hold on;
xline(mean(f_vals), 'r-', 'LineWidth', 2, 'DisplayName', sprintf('均值=%.4f', mean(f_vals)));
hold off;
xlabel('f(x) 值'); ylabel('概率密度');
title(sprintf('f(X) 的分布（X~Uniform[%g,%g]）', a, b));
legend;
grid on;

% 子图3：两种方法的误差收敛对比（双对数图）
subplot(2, 3, [3, 6]);
loglog(N_range, err_hit_arr,  'r-', 'LineWidth', 1.5, 'DisplayName', '投点法误差');
hold on;
loglog(N_range, err_mean_arr, 'b-', 'LineWidth', 1.5, 'DisplayName', '均值法误差');
% 理论 O(1/√N) 参考线
loglog(N_range, 2./sqrt(N_range), 'k--', 'LineWidth', 1, 'DisplayName', 'O(1/√N)');
hold off;
xlabel('随机点数 N（对数刻度）');
ylabel('误差（对数刻度）');
title('两种方法误差随 N 收敛（双对数图）');
legend('Location', 'northeast');
grid on;
% 双对数图上斜率=-0.5 对应 O(N^{-1/2}) = O(1/√N)

% 子图4：估计值随 N 的变化
subplot(2, 3, 4);
cum_mean_plot = (b-a) * cumsum(f_vals_all) ./ (1:N);   % 实时估计值
semilogx(1:N, cum_mean_plot, 'b-', 'LineWidth', 1.0);
hold on;
semilogx([1, N], [I_exact, I_exact], 'r--', 'LineWidth', 1.5);
hold off;
xlabel('N（对数刻度）'); ylabel('估计值');
title('样本均值法：估计值随 N 收敛');
xlim([1, N]);
legend({'估计值', '精确值'}, 'Location', 'northeast');
grid on;

% 子图5：误差汇总条形图
subplot(2, 3, 5);
final_errors = [abs(I_hit-I_exact), abs(I_mean-I_exact)];
bar(final_errors, 'FaceColor', [0.3 0.7 0.5]);
set(gca, 'XTickLabel', {'投点法', '均值法'});
set(gca, 'YScale', 'log');   % y轴对数刻度
for k = 1:2
    text(k, final_errors(k)*1.5, sprintf('%.2e', final_errors(k)), ...
        'HorizontalAlignment', 'center', 'FontSize', 9);
end
ylabel('|估计误差|（对数刻度）');
title(sprintf('两种方法最终误差对比（N=%.0e）', N));
grid on;

sgtitle(sprintf('Monte Carlo 积分：∫_%.0f^%.0f f(x)dx，N=%.0e', a, b, N));

%% ---- 汇总 ----
fprintf('==============================================\n');
fprintf('  Monte Carlo 积分结果汇总（N=%.0e）\n', N);
fprintf('==============================================\n');
fprintf('精确值     = %.10f\n', I_exact);
fprintf('投点法     = %.10f  误差 = %.4e\n', I_hit, abs(I_hit-I_exact));
fprintf('样本均值法 = %.10f  误差 = %.4e\n', I_mean, abs(I_mean-I_exact));
fprintf('MATLAB内置 = %.10f  误差 = %.4e\n', integral(f,a,b), abs(integral(f,a,b)-I_exact));
fprintf('\n推荐：优先使用样本均值法（方差更小，代码更简单）\n');