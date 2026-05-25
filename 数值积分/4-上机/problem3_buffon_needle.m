%% =========================================================
%  problem3_buffon_needle.m
%  功能：蒲丰投针（Buffon's Needle）Monte Carlo 模拟估计 π
%
%  【历史背景】
%  1777年，法国数学家 Georges-Louis Leclerc（蒲丰伯爵）提出：
%  在铺有平行线的地板上，随机投下一根针，
%  针与某条平行线相交的概率是多少？
%  这是最早的概率与几何相结合的问题之一。
%
%  【问题建模】
%  设平行线间距为 d，针的长度为 l（l ≤ d，短针情形）。
%  随机变量：
%    y  ∼ Uniform[0, d]：针的中心距最近下方平行线的距离
%    θ  ∼ Uniform[0, π]：针与平行线的夹角
%  针与平行线相交的条件：
%    针在 y 方向的"投影半长" = (l/2)·sin(θ)
%    相交 ⟺ y < (l/2)·sin(θ)  或  y > d - (l/2)·sin(θ)
%
%  【解析推导——为什么能估计 π？】
%  P(相交) = ∫₀^π ∫₀^d [y < l/2·sinθ 或 y > d-l/2·sinθ] dy·dθ / (π·d)
%           = 2/π ∫₀^π (l/2·sinθ) dθ / d
%           = 2l/(πd) · ∫₀^π sinθ dθ / π
%  实际上精确推导：
%    P = 2/(πd) ∫₀^π (l/2)sinθ dθ = 2l/(πd) · [-cosθ]₀^π / 1 ...
%    P = 2l/(πd)
%  因此：π ≈ 2l·N / (d·crossings)
%  其中 N 为投针总次数，crossings 为相交次数。
%
%  【Monte Carlo 方法核心思想】
%  大数定律：实验频率在实验次数→∞时收敛于真实概率。
%  Monte Carlo 误差：标准差 σ ∝ 1/√N，
%  即误差每减小10倍，需要将实验次数增加100倍。
%  =========================================================

clear; close all; clc;

% 设置随机数种子（使结果可复现）
% rng(seed)：随机数生成器（Random Number Generator）初始化
% 固定种子后，每次运行产生完全相同的随机序列
rng(42);

fprintf('==============================================\n');
fprintf('   第三题：蒲丰投针模拟估计 π\n');
fprintf('==============================================\n\n');


%% =========================================================
%  参数设置
%  =========================================================

l = 1.0;   % 针的长度（length of needle）
d = 2.0;   % 平行线间距（distance between parallel lines）
           % 要求 l ≤ d（短针情形，公式 P = 2l/(πd) 成立）
           % 若 l > d（长针），公式更复杂（需分情形讨论）

fprintf('参数设置：\n');
fprintf('  针长 l = %.1f，平行线间距 d = %.1f\n', l, d);
fprintf('  理论相交概率 P = 2l/(πd) = %.6f\n', 2*l/(pi*d));
fprintf('  精确参考值 π = %.10f\n\n', pi);

assert(l <= d, '需要 l ≤ d（短针情形）');  % assert：条件不满足时报错


%% =========================================================
%  第一部分：大规模模拟（一次性投入全部针）
%  =========================================================

fprintf('--- 第一部分：大规模模拟（N = 10^7）---\n\n');

N_total = 1e7;   % 总投针次数（1千万次）

% 生成随机针的中心位置和角度
% rand(1, N)：生成 1×N 的均匀分布随机数，范围 [0, 1)
y_center = d * rand(1, N_total);    % 中心到最近下方线的距离，均匀分布于 [0, d)
theta    = pi * rand(1, N_total);   % 针与平行线夹角，均匀分布于 [0, π)

% 判断相交：针的中心到下方平行线距离 y < l/2·sin(θ)，或到上方线距离 < l/2·sin(θ)
% 即：y < l/2·sinθ  或  d-y < l/2·sinθ
half_proj = (l / 2) * sin(theta);   % 针在垂直方向的投影半长（逐元素·）

% 逻辑运算：| 是逐元素"或"，< 是逐元素比较，返回逻辑数组（0或1）
cross_flag = (y_center < half_proj) | (y_center > d - half_proj);

% sum()：对逻辑数组求和，即统计 true 的个数（相交次数）
N_cross = sum(cross_flag);

% 由 P = 2l/(πd) 反推 π
pi_estimate = 2 * l * N_total / (d * N_cross);

fprintf('投针次数：N = %.0e\n', N_total);
fprintf('相交次数：%d\n', N_cross);
fprintf('实验相交频率：%.6f（理论值 %.6f）\n', N_cross/N_total, 2*l/(pi*d));
fprintf('π 估计值：%.10f\n', pi_estimate);
fprintf('绝对误差：%.6e\n', abs(pi_estimate - pi));
fprintf('相对误差：%.4f %%\n\n', abs(pi_estimate - pi) / pi * 100);


%% =========================================================
%  第二部分：误差随实验次数的收敛分析
%  =========================================================
%
%  【Monte Carlo 收敛理论】
%  若真实概率为 p，N 次实验中相交频率 p̂ = X/N（X∼Binomial(N,p)），
%  则由中心极限定理：
%    (p̂ - p) / √(p(1-p)/N) → N(0,1)
%  π 的估计量 π̂ = 2l/(d·p̂) 的标准差约为：
%    σ_{π̂} ≈ π · √(p(1-p)/N) / p = π·√((1-p)/Np)
%  即误差 σ ∝ 1/√N，Monte Carlo 收敛速度与维数无关（维数灾免疫）
%  =========================================================

fprintf('--- 第二部分：误差随实验次数的收敛 ---\n\n');

% 实验规模序列（对数均匀分布）
N_seq = round(logspace(2, 7, 60));   % 从 100 到 10^7，共60个测试点

pi_seq  = zeros(1, length(N_seq));   % 各规模对应的 π 估计值
err_seq = zeros(1, length(N_seq));   % 对应的绝对误差

% 预先生成所有随机数（避免循环内多次调用 rand，提高效率）
y_all     = d * rand(1, N_total);
theta_all = pi * rand(1, N_total);
half_all  = (l/2) * sin(theta_all);
cross_all = (y_all < half_all) | (y_all > d - half_all);  % 全部结果
cumcross  = cumsum(cross_all);   % cumsum：累积求和，cumcross(k)=前k次实验的相交总数

for k = 1:length(N_seq)
    N_k = N_seq(k);
    % 用前 N_k 次投针的累积结果估计 π（无需重新生成随机数）
    cross_k = cumcross(N_k);   % 前 N_k 次中的相交次数
    if cross_k > 0
        pi_seq(k) = 2 * l * N_k / (d * cross_k);
    else
        pi_seq(k) = Inf;   % 尚无相交，无法估计（极端情形，N 很小时偶发）
    end
    err_seq(k) = abs(pi_seq(k) - pi);
end

% 理论收敛曲线：σ ≈ π·√((1-p₀)/(N·p₀))，p₀ = 2l/(πd)
p0 = 2*l/(pi*d);   % 真实相交概率
theory_std = pi * sqrt((1-p0) ./ (N_seq * p0));

fprintf('%-12s  %-15s  %-15s\n', 'N', 'π估计值', '绝对误差');
fprintf('%s\n', repmat('-', 1, 45));
display_idx = round(linspace(1, length(N_seq), 12));   % 选12个有代表性的点显示
for k = display_idx
    fprintf('%-12d  %-15.8f  %-15.4e\n', N_seq(k), pi_seq(k), err_seq(k));
end
fprintf('\n');


%% =========================================================
%  第三部分：多次重复实验——验证误差分布
%  =========================================================
%
%  中心极限定理保证：大量重复实验时，π̂ 的分布趋近正态分布
%  本部分用 M 次独立实验（每次 N 根针）验证这一理论
%  =========================================================

fprintf('--- 第三部分：重复实验验证误差分布（中心极限定理）---\n\n');

M_repeat = 500;    % 重复实验次数
N_each   = 5000;   % 每次实验投针数

pi_repeat = zeros(1, M_repeat);   % 每次实验的 π 估计值

for m = 1:M_repeat
    y_m     = d * rand(1, N_each);
    theta_m = pi * rand(1, N_each);
    hp_m    = (l/2) * sin(theta_m);
    cross_m = sum((y_m < hp_m) | (y_m > d - hp_m));
    if cross_m > 0
        pi_repeat(m) = 2*l*N_each / (d*cross_m);
    else
        pi_repeat(m) = NaN;   % NaN（Not a Number）：MATLAB 中的缺失值标记
    end
end

% 去掉 NaN 值
pi_repeat = pi_repeat(~isnan(pi_repeat));
% ~isnan()：逻辑非+isnan，过滤掉 NaN 值

fprintf('%d 次独立实验（每次 N=%d 根针）统计：\n', M_repeat, N_each);
fprintf('  π 估计均值：  %.8f\n', mean(pi_repeat));
% mean(v)：计算向量 v 的算术平均值
fprintf('  π 真实值：    %.8f\n', pi);
fprintf('  估计标准差：  %.6f\n', std(pi_repeat));
% std(v)：计算向量 v 的样本标准差（分母为 N-1）
fprintf('  理论标准差：  %.6f\n', pi*sqrt((1-p0)/(N_each*p0)));
fprintf('\n');


%% =========================================================
%  第四部分：可视化
%  =========================================================

figure('Name', '第三题：蒲丰投针 Monte Carlo', 'NumberTitle', 'off', ...
       'Position', [50 50 1100 700]);


% ---- 子图1：针的实际模拟示意图 ----
subplot(2, 3, 1);
n_show = 200;   % 显示 200 根针（不宜太多，否则太密）

% 随机生成要显示的针
y_show     = d * rand(1, n_show);
theta_show = pi * rand(1, n_show);
hp_show    = (l/2) * sin(theta_show);
cross_show = (y_show < hp_show) | (y_show > d - hp_show);

hold on;
% 绘制平行线（假设显示4条，y=0, d, 2d, 3d）
n_lines = 4;
for line_k = 0:n_lines-1
    % plot([x1,x2], [y1,y2])：绘制从(x1,y1)到(x2,y2)的线段
    plot([0, l*3], [line_k*d, line_k*d], 'k-', 'LineWidth', 1.5);
end

% 绘制每根针
for k = 1:n_show
    % 针的两端点（根据中心坐标和角度计算）
    y_k = y_show(k) + floor(rand()*n_lines)*d;  % 随机分配到哪两条线之间
    y_k = mod(y_k, n_lines*d);                   % mod：取余（保持在显示区域内）
    x_center = l * rand();    % 随机水平位置
    dx = (l/2) * cos(theta_show(k));   % 水平投影半长
    dy = (l/2) * sin(theta_show(k));   % 垂直投影半长

    x1 = x_center - dx;  x2 = x_center + dx;
    y1 = y_k - dy;        y2 = y_k + dy;

    if cross_show(k)
        plot([x1,x2], [y1,y2], 'r-', 'LineWidth', 1.2);  % 相交针：红色
    else
        plot([x1,x2], [y1,y2], 'b-', 'LineWidth', 0.8);  % 不相交：蓝色
    end
end
axis equal;
xlim([0, l*3]);
ylim([0, n_lines*d]);
title(sprintf('投针示意（红=相交，蓝=不相交）\nn=%d，相交%d根', n_show, sum(cross_show)));
xlabel('x');  ylabel('y');


% ---- 子图2：π 估计值随实验次数的收敛 ----
subplot(2, 3, 2);
semilogx(N_seq, pi_seq, 'b-', 'LineWidth', 1.2, 'DisplayName', 'π 估计值');
hold on;
% 精确值水平参考线
yline(pi, 'r--', 'LineWidth', 1.5, 'DisplayName', '真实 π');
% 用阴影显示 ±2σ 置信带（约 95% 置信区间）
% fill：[x正向, x反向], [上界, 下界] 构成填充多边形
fill([N_seq, fliplr(N_seq)], ...
     [pi+2*theory_std, fliplr(pi-2*theory_std)], ...
     [1 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', '±2σ 置信带');
% 'FaceAlpha'：透明度，0=完全透明，1=不透明；'EdgeColor','none'：无边框线
xlabel('投针次数 N（对数刻度）');
ylabel('π 估计值');
title('π 估计值随实验次数的收敛');
legend('Location', 'NorthEast');
grid on;


% ---- 子图3：误差随 N 的收敛（双对数，验证 O(1/√N)）----
subplot(2, 3, 3);
loglog(N_seq, err_seq, 'b-', 'LineWidth', 1.2, 'DisplayName', '实际误差');
hold on;
loglog(N_seq, theory_std, 'r--', 'LineWidth', 1.5, 'DisplayName', '理论 σ = C/√N');
xlabel('N（对数刻度）');
ylabel('|π̂ - π|（对数刻度）');
title('收敛速度验证：误差 ∝ 1/√N');
legend('Location', 'NorthEast');
grid on;


% ---- 子图4：重复实验 π 估计值的分布直方图 ----
subplot(2, 3, 4);
histogram(pi_repeat, 30, 'Normalization', 'pdf', ...
    'FaceColor', [0.4 0.6 0.8], 'EdgeColor', 'white');
% histogram(data, nbins, ...)：绘制直方图
% 'Normalization','pdf'：纵轴归一化为概率密度（面积=1）
hold on;
% 叠加理论正态分布曲线
mu_theory = pi;
sig_theory = pi*sqrt((1-p0)/(N_each*p0));
x_norm = linspace(pi-4*sig_theory, pi+4*sig_theory, 200);
% normpdf(x, mu, sigma)：正态分布概率密度函数（需 Statistics Toolbox）
% 手动计算以避免依赖 Toolbox：
y_norm = (1/(sig_theory*sqrt(2*pi))) * exp(-0.5*((x_norm-mu_theory)/sig_theory).^2);
plot(x_norm, y_norm, 'r-', 'LineWidth', 2, 'DisplayName', '理论正态分布');
xline(pi, 'k--', 'LineWidth', 1.5, 'DisplayName', '真实 π');
xlabel('π 估计值');
ylabel('概率密度');
title(sprintf('π 估计值分布（M=%d次，N=%d）\n均值=%.4f，std=%.4f', ...
    M_repeat, N_each, mean(pi_repeat), std(pi_repeat)));
legend('Location', 'NorthEast', 'FontSize', 8);


% ---- 子图5：被积函数（几何概率密度）----
subplot(2, 3, 5);
theta_axis = linspace(0, pi, 200);
% 对于给定角度 θ，相交概率 = l·sin(θ)/d（d 上均匀的针穿过一条线）
P_theta = (l/d) * sin(theta_axis);
fill([theta_axis, fliplr(theta_axis)], [P_theta, zeros(1,200)], ...
     [0.8 0.9 0.7], 'EdgeColor', 'g', 'LineWidth', 1.2);
xlabel('θ（角度，弧度）');
ylabel('P(相交|θ) = l·sinθ/d');
title(sprintf('条件相交概率\n面积 = 2l/(πd) = %.4f', 2*l/(pi*d)));
xticks([0 pi/4 pi/2 3*pi/4 pi]);
xticklabels({'0', 'π/4', 'π/2', '3π/4', 'π'});
grid on;


% ---- 子图6：累计相交次数 vs 理论期望 ----
subplot(2, 3, 6);
n_show_seq = 1:10000;   % 显示前 10000 次投针的累积过程
% cumsum：累积和，cumcross(k) = 前k次实验的相交总次数
cum_cross_show = cumcross(n_show_seq);
cum_pi_show    = 2*l*n_show_seq ./ (d * cum_cross_show);
cum_pi_show(cum_cross_show == 0) = NaN;   % 避免除以零
plot(n_show_seq, cum_pi_show, 'b-', 'LineWidth', 0.8);
hold on;
yline(pi, 'r--', 'LineWidth', 1.5);
xlabel('投针次数 N');
ylabel('π 累积估计值');
title('π 估计值实时收敛过程（前10000次）');
grid on;

sgtitle('第三题：蒲丰投针 Monte Carlo 模拟估计 π');

fprintf('图形已生成，请查看弹出的图形窗口。\n');
fprintf('\n程序运行完毕。\n');
