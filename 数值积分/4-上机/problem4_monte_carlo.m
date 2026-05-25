%% =========================================================
%  problem4_monte_carlo.m
%  功能：分别用随机投点法和样本均值法计算定积分
%    I = ∫₀¹ e^{-x} dx
%  并与精确值比较，同时画出积分的面积区域
%
%  【精确值】
%  I = ∫₀¹ e^{-x} dx = [-e^{-x}]₀¹ = 1 - e^{-1} ≈ 0.63212055882...
%
%  【两种 Monte Carlo 方法的对比】
%  ─────────────────────────────────────────────────────────
%  方法               原理                     方差
%  ─────────────────────────────────────────────────────────
%  随机投点法（Hit-Miss）生成均匀点判断落在曲线下方  较大（取决于目标区域形状）
%  样本均值法（Sample Mean）直接对函数值取平均       较小（等于 Var[f(X)]）
%  ─────────────────────────────────────────────────────────
%  两种方法的误差均为 O(1/√N)，但样本均值法方差更小，效率更高。
%  =========================================================

clear; close all; clc;

rng(2024);   % 固定随机种子，确保结果可复现

fprintf('==============================================\n');
fprintf('   第四题：Monte Carlo 计算 I = ∫₀¹ e^{-x} dx\n');
fprintf('==============================================\n\n');

% 被积函数 f(x) = e^{-x}
% exp(-x)：MATLAB 中自然指数函数，支持向量/矩阵输入
f = @(x) exp(-x);

% 积分精确值
I_exact = 1 - exp(-1);   % = 1 - 1/e
fprintf('精确值：I = 1 - e^{-1} = %.15f\n\n', I_exact);

% 总样本数
N_total = 1e6;   % 100万个点


%% =========================================================
%  第一部分：随机投点法（Hit-or-Miss / Rejection Sampling）
%  =========================================================
%
%  【算法原理】
%  在包含积分区域的矩形 [a,b]×[0,M] 内均匀撒点，
%  统计落在曲线 y=f(x) 下方的点的比例：
%    I ≈ (b-a)·M · (落在曲线下方的点数 / 总点数)
%
%  本题：[0,1]×[0,1]（因为 f(x)=e^{-x} ∈ (e^{-1}, 1) ⊂ [0,1]）
%  矩形面积 = 1×1 = 1
%  所以：I ≈ 命中次数 / N（矩形面积=1，简化公式）
%
%  【几何含义】
%  "命中"的点均匀散布在积分区域（曲线下方阴影部分）内，
%  积分值 = 阴影面积 = 矩形面积 × 命中比例
%
%  【方差分析】
%  设 Z = 1(Y ≤ f(X))（示性函数，命中为1，未命中为0）
%  E[Z] = I（期望 = 积分值）
%  Var[Z] = I(1-I)（Bernoulli 分布方差）
%  对 I ≈ 0.632 估计量的标准差 ≈ √(I(1-I)/N) ≈ √(0.232/N)
%  =========================================================

fprintf('--- 第一部分：随机投点法（Hit-or-Miss）---\n\n');

% 在 [0,1]×[0,1] 内均匀生成 N_total 个点 (x,y)
x_hit = rand(1, N_total);   % x 坐标，均匀分布于 [0,1)
y_hit = rand(1, N_total);   % y 坐标，均匀分布于 [0,1)
% rand(1, N)：生成 1×N 的 [0,1) 均匀随机数

% 判断点是否落在曲线下方：y ≤ f(x) = e^{-x}
% 返回逻辑向量（true/false），sum() 统计 true 的个数
hits = sum(y_hit <= f(x_hit));

% 积分估计：命中比例 × 矩形面积（此处矩形面积=1）
I_hitormiss = hits / N_total;

fprintf('随机投点法（N=%.0e）：\n', N_total);
fprintf('  命中次数：%d / %d（比例 %.6f）\n', hits, N_total, hits/N_total);
fprintf('  I 估计值：%.10f\n', I_hitormiss);
fprintf('  绝对误差：%.6e\n', abs(I_hitormiss - I_exact));
fprintf('  理论标准差：%.6e\n\n', sqrt(I_exact*(1-I_exact)/N_total));


%% =========================================================
%  第二部分：样本均值法（Sample Mean Method）
%  =========================================================
%
%  【算法原理】
%  利用积分的概率解释：
%    I = ∫_a^b f(x)·(1/(b-a)) dx · (b-a) = E[f(X)] · (b-a)
%  其中 X ∼ Uniform[a,b]。
%  故 I ≈ (b-a)/N · ∑_{i=1}^N f(Xᵢ)，Xᵢ ∼ Uniform[a,b]
%
%  本题 a=0, b=1，故：I ≈ (1/N) ∑ f(Xᵢ) = (1/N) ∑ e^{-Xᵢ}
%
%  【与投点法的方差对比】
%  样本均值法：Var[f(X)] = E[f²] - (E[f])²
%    = ∫₀¹ e^{-2x}dx - I² = (1-e^{-2})/2 - (1-e^{-1})² ≈ 0.0398
%  标准差 ≈ √(Var[f]/N) ≈ √(0.0398/N)
%
%  投点法：Var[Z] = I(1-I) ≈ 0.232
%  标准差 ≈ √(0.232/N)
%
%  结论：样本均值法方差约为投点法的 0.0398/0.232 ≈ 1/6，效率高约6倍
%  =========================================================

fprintf('--- 第二部分：样本均值法（Sample Mean）---\n\n');

% 在 [0,1] 上均匀生成 N_total 个采样点 x
x_sample = rand(1, N_total);

% 积分估计：样本 f 值的平均值（× (b-a) = 1）
I_samplemean = mean(f(x_sample));
% mean(v)：计算向量 v 的算术平均值，等价于 sum(v)/length(v)

% 样本方差（用于验证理论值）
var_sample = var(f(x_sample));
% var(v)：计算向量 v 的样本方差（分母为 N-1）

fprintf('样本均值法（N=%.0e）：\n', N_total);
fprintf('  I 估计值：  %.10f\n', I_samplemean);
fprintf('  绝对误差：  %.6e\n', abs(I_samplemean - I_exact));
fprintf('  样本标准差：%.6e\n', sqrt(var_sample/N_total));
fprintf('  理论标准差：%.6e\n\n', sqrt(((1-exp(-2))/2 - I_exact^2)/N_total));


%% =========================================================
%  第三部分：两种方法的收敛分析（误差随 N 变化）
%  =========================================================

fprintf('--- 第三部分：收敛性对比分析 ---\n\n');

N_list = round(logspace(2, 6, 50));   % 从 100 到 10^6

err_hit = zeros(1, length(N_list));   % 投点法误差
err_smp = zeros(1, length(N_list));   % 样本均值法误差

% 预先生成所有随机数，累积使用（避免不同 N 之间的随机波动干扰比较）
x_all_h = rand(1, N_total);   % 投点法 x
y_all_h = rand(1, N_total);   % 投点法 y
x_all_s = rand(1, N_total);   % 样本均值法 x
% 预计算命中标志和函数值
hit_all   = (y_all_h <= f(x_all_h));   % 逻辑向量：第 i 个点是否命中
fval_all  = f(x_all_s);                % 样本函数值

% cumsum：计算累积和，使得不同 N 只需查表而无需重新求和
cumhit    = cumsum(hit_all);           % 累积命中次数
cumfval   = cumsum(fval_all);          % 累积函数值之和

for k = 1:length(N_list)
    N_k = N_list(k);
    I_h = cumhit(N_k) / N_k;          % 投点法估计
    I_s = cumfval(N_k) / N_k;         % 样本均值法估计
    err_hit(k) = abs(I_h - I_exact);
    err_smp(k) = abs(I_s - I_exact);
end

% 理论标准差曲线
std_hit = sqrt(I_exact * (1 - I_exact) ./ N_list);   % 投点法理论标准差
var_f   = (1 - exp(-2))/2 - I_exact^2;               % Var[f(X)]
std_smp = sqrt(var_f ./ N_list);                       % 样本均值法理论标准差

fprintf('效率比较（理论方差比）：\n');
fprintf('  投点法方差：  Var[Z] = I(1-I) = %.6f\n', I_exact*(1-I_exact));
fprintf('  样本均值法方差：Var[f(X)] = %.6f\n', var_f);
fprintf('  样本均值法效率是投点法的 %.1f 倍\n\n', ...
    I_exact*(1-I_exact)/var_f);


%% =========================================================
%  第四部分：综合结果汇总
%  =========================================================

fprintf('==============================================\n');
fprintf('           综合结果汇总\n');
fprintf('==============================================\n\n');

fprintf('%-25s  %-18s  %-12s\n', '方法', '估计值', '误差');
fprintf('%s\n', repmat('-', 1, 58));
fprintf('%-25s  %-18.10f  %-12.4e\n', '随机投点法', I_hitormiss, abs(I_hitormiss-I_exact));
fprintf('%-25s  %-18.10f  %-12.4e\n', '样本均值法', I_samplemean, abs(I_samplemean-I_exact));
fprintf('%-25s  %-18.10f  %-12s\n', '精确值 1-e^{-1}', I_exact, '0');
fprintf('\n');


%% =========================================================
%  第五部分：图形可视化
%  =========================================================

figure('Name', '第四题：Monte Carlo 计算积分', 'NumberTitle', 'off', ...
       'Position', [50 50 1100 800]);


% ---- 子图1：积分区域面积示意图 ----
subplot(2, 3, 1);
x_curve = linspace(0, 1, 300);   % 曲线绘制点
y_curve = f(x_curve);

% 绘制曲线 y = e^{-x} 下的面积（填充区域）
fill([x_curve, fliplr(x_curve)], [y_curve, zeros(1,300)], ...
     [0.65 0.85 0.95], 'EdgeColor', 'b', 'LineWidth', 2, ...
     'DisplayName', '积分区域（I = 阴影面积）');
% fill：绘制多边形并填充颜色
% 第一组 x 正向（0→1），第二组 x 反向（1→0），构成闭合多边形
% zeros(1,300)：下边界（x轴）对应的 y=0
hold on;
plot(x_curve, y_curve, 'b-', 'LineWidth', 2, 'DisplayName', 'y = e^{-x}');
% 在图上标注精确积分值
area_text = sprintf('面积 = %.6f', I_exact);
text(0.35, 0.3, area_text, 'FontSize', 12, 'Color', 'b', 'FontWeight', 'bold');
% text(x, y, str)：在坐标 (x,y) 处放置文字标注
xlabel('x');
ylabel('y = e^{-x}');
title('积分区域（曲线下阴影面积）');
xlim([0, 1]);  ylim([0, 1.1]);
legend('Location', 'NorthEast', 'FontSize', 8);
grid on;


% ---- 子图2：随机投点法示意图（显示前5000个点）----
subplot(2, 3, 2);
N_show = 5000;   % 只显示5000个点，避免画面太密

x_show = x_all_h(1:N_show);
y_show = y_all_h(1:N_show);
hit_show = hit_all(1:N_show);

% 画曲线下的填充区域（参考）
fill([x_curve, fliplr(x_curve)], [y_curve, zeros(1,300)], ...
     [0.9 0.95 0.98], 'EdgeColor', 'none');
hold on;

% 命中点（蓝色）和未命中点（红色）
% x_show(~hit_show)：逻辑索引，选出未命中的 x 坐标
plot(x_show(~hit_show), y_show(~hit_show), 'r.', 'MarkerSize', 2, ...
     'DisplayName', '未命中');
plot(x_show(hit_show), y_show(hit_show), 'b.', 'MarkerSize', 2, ...
     'DisplayName', '命中');
plot(x_curve, y_curve, 'k-', 'LineWidth', 1.5, 'DisplayName', 'y = e^{-x}');

xlabel('x');  ylabel('y');
title(sprintf('随机投点法示意（N=%d）\n命中率=%.4f，I≈%.4f', ...
    N_show, sum(hit_show)/N_show, sum(hit_show)/N_show));
legend('Location', 'NorthEast', 'FontSize', 7);
xlim([0,1]);  ylim([0,1]);


% ---- 子图3：样本均值法示意图（直方图显示 f(x) 的样本分布）----
subplot(2, 3, 3);
% histogram：绘制 f(x) 值的分布直方图
% x_sample 均匀分布于 [0,1]，f(x)=e^{-x} 在 [e^{-1}, 1] 上取值
histogram(f(x_all_s(1:10000)), 30, 'Normalization', 'pdf', ...
    'FaceColor', [0.4 0.7 0.4], 'EdgeColor', 'white');
hold on;

% 理论概率密度：X 均匀分布 → f(X)=e^{-X} 的 PDF
% 对于 f(x)=e^{-x}，X∼U[0,1]，f(X) 的 PDF：
% f(X) = t ⟺ X = -ln(t)，dX/dt = 1/t，P(f(X)=t) = 1/t（t∈[e^{-1}, 1]）
t_range = linspace(exp(-1), 1, 200);
pdf_fX = 1 ./ t_range;   % f(X) 的概率密度函数
plot(t_range, pdf_fX, 'r-', 'LineWidth', 2, 'DisplayName', '理论密度');

% 标注均值（即积分值）
xline(I_exact, 'b--', 'LineWidth', 1.5, 'DisplayName', sprintf('均值=%.4f',I_exact));
xlabel('f(x) = e^{-x} 的取值');
ylabel('概率密度');
title(sprintf('样本均值法：f(X) 的分布\n样本均值 = ∫e^{-x}dx ≈ %.6f', ...
    mean(f(x_all_s(1:10000)))));
legend('Location', 'NorthEast', 'FontSize', 7);


% ---- 子图4：两种方法误差随 N 的收敛对比（双对数图）----
subplot(2, 3, 4);
loglog(N_list, err_hit, 'r-', 'LineWidth', 1.5, 'DisplayName', '投点法误差');
hold on;
loglog(N_list, err_smp, 'b-', 'LineWidth', 1.5, 'DisplayName', '样本均值法误差');
loglog(N_list, std_hit, 'r--', 'LineWidth', 1.0, 'DisplayName', '投点法理论σ');
loglog(N_list, std_smp, 'b--', 'LineWidth', 1.0, 'DisplayName', '样本均值法理论σ');

% 添加 O(1/√N) 参考线
N_ref = [1e2, 1e6];
loglog(N_ref, 0.5./sqrt(N_ref), 'k:', 'LineWidth', 1.2, 'DisplayName', 'O(N^{-1/2})');

xlabel('样本数 N（对数刻度）');
ylabel('误差（对数刻度）');
title('两种方法收敛速度比较');
legend('Location', 'SouthWest', 'FontSize', 7);
grid on;


% ---- 子图5：两种方法实时估计值（前50000次）----
subplot(2, 3, 5);
N_live = 50000;   % 显示前50000次的动态过程
n_axis = 1:N_live;

% 利用已计算的累积量，直接取出前 N_live 点的累积均值
live_hit  = cumhit(1:N_live)  ./ n_axis;   % 投点法累积均值
live_samp = cumfval(1:N_live) ./ n_axis;   % 样本均值法累积均值

plot(n_axis, live_hit,  'r-', 'LineWidth', 0.8, 'DisplayName', '投点法');
hold on;
plot(n_axis, live_samp, 'b-', 'LineWidth', 0.8, 'DisplayName', '样本均值法');
yline(I_exact, 'k--', 'LineWidth', 1.5, 'DisplayName', sprintf('精确值=%.6f',I_exact));
xlabel('样本数 N');
ylabel('I 累积估计值');
title('实时收敛过程（前 50000 次）');
legend('Location', 'NorthEast', 'FontSize', 8);
grid on;
ylim([I_exact-0.05, I_exact+0.05]);   % 限制纵轴范围以便观察收敛


% ---- 子图6：多次重复实验的估计值分布 ----
subplot(2, 3, 6);
M_rep = 1000;   % 重复次数
N_each = 1000;  % 每次样本量

pi_hit_rep  = zeros(1, M_rep);   % 投点法重复结果
pi_smp_rep  = zeros(1, M_rep);   % 样本均值法重复结果

for m = 1:M_rep
    xh = rand(1, N_each);  yh = rand(1, N_each);
    h_m = sum(yh <= f(xh));
    pi_hit_rep(m) = h_m / N_each;

    xs = rand(1, N_each);
    pi_smp_rep(m) = mean(f(xs));
end

% 绘制两种方法估计值的分布直方图（叠加比较）
histogram(pi_hit_rep, 30, 'Normalization', 'pdf', ...
    'FaceColor', [1 0.6 0.6], 'EdgeColor', 'white', 'FaceAlpha', 0.6, ...
    'DisplayName', sprintf('投点法(σ=%.4f)', std(pi_hit_rep)));
hold on;
histogram(pi_smp_rep, 30, 'Normalization', 'pdf', ...
    'FaceColor', [0.6 0.7 1], 'EdgeColor', 'white', 'FaceAlpha', 0.6, ...
    'DisplayName', sprintf('样本均值法(σ=%.4f)', std(pi_smp_rep)));
xline(I_exact, 'k--', 'LineWidth', 2, 'DisplayName', '精确值');
xlabel('I 估计值');
ylabel('概率密度');
title(sprintf('估计值分布对比（M=%d次，N=%d）\n样本均值法方差更小', M_rep, N_each));
legend('Location', 'NorthWest', 'FontSize', 7);
% 可以清楚看到：样本均值法的分布更窄（方差更小），估计更稳定

sgtitle('第四题：Monte Carlo 计算 I = ∫₀¹ e^{-x}dx');

fprintf('图形已生成，请查看弹出的图形窗口。\n');
fprintf('\n程序运行完毕。\n');
