%% =========================================================
%  模板T10：蒲丰投针（Buffon's Needle）—— Monte Carlo 估计 π
%
%  【数学背景】
%  平面上画一组平行线，间距为 d。
%  向平面随机投一根长度为 l 的针（l ≤ d）。
%  理论上，针与某条直线相交的概率为：
%    P = 2l / (πd)
%  因此可反推：π = 2l / (d·P)
%  用大量随机投针实验估计概率 P，从而估计 π。
%
%  【随机模拟的两个随机量】
%  设针的中心到最近直线的距离为 X（均匀分布在 [0, d/2]），
%  针与直线的夹角为 θ（均匀分布在 [0, π/2]，利用对称性）。
%  针与直线相交的条件：X ≤ l/2·sin(θ)
%
%  【收敛速度】
%  Monte Carlo 方法误差 = O(1/√N)，N 增大100倍，误差才减小10倍。
%  这比数值积分慢得多，但对高维问题无可替代。
%  =========================================================

clear; close all; clc;

%% ===== ① 修改这里：实验参数 =====
l = 1;       % 针的长度
d = 2;       % 平行线间距（要求 l ≤ d）
N = 1e6;     % 投针总次数（越大越精确，但时间也越长）
%% =====================================

fprintf('蒲丰投针 Monte Carlo 估计 π\n');
fprintf('针长 l=%g，线距 d=%g，投针次数 N=%.0e\n\n', l, d, N);

%% ---- 第一步：大规模模拟 ----
% 生成 N 根随机针（向量化，一次性生成所有随机数）
% rand(1, N)：N 个在 [0,1] 均匀分布的随机数（行向量）
x_center = (d/2) * rand(1, N);   % 针中心到最近直线的距离，均匀分布在 [0, d/2]
theta = (pi/2) * rand(1, N);     % 针与直线的夹角，均匀分布在 [0, π/2]

% 判断是否相交：X ≤ (l/2)·sin(θ)
% 结果是逻辑向量（0 或 1）
cross = (x_center <= l/2 * sin(theta));
% sin(theta)：对向量 theta 每个元素求正弦（向量化）
% <= 是逐元素比较，返回等大小的逻辑向量

% 统计相交次数并估计概率
n_cross = sum(cross);   % sum 对逻辑向量求和 = true 的个数
P_hat = n_cross / N;    % 相交概率的估计值

% 反推 π = 2l / (d·P)
pi_est = 2 * l / (d * P_hat);

fprintf('相交次数：%d / %d\n', n_cross, N);
fprintf('估计概率 P ≈ %.6f（理论值 2l/(πd) = %.6f）\n', P_hat, 2*l/(pi*d));
fprintf('π 估计值 = %.8f\n', pi_est);
fprintf('π 精确值 = %.8f\n', pi);
fprintf('绝对误差 = %.4e\n\n', abs(pi_est - pi));

%% ---- 第二步：收敛性分析（误差随 N 的变化）----
fprintf('--- 收敛性分析：误差随实验次数 N 的变化 ---\n');

% 对同一组随机数，累计计算估计值（避免重复生成随机数）
N_checkpoints = round(logspace(2, log10(N), 50));   % 50个检查点
% logspace(a, b, n)：生成从 10^a 到 10^b 的 n 个对数均匀间隔的值
pi_running = zeros(size(N_checkpoints));
err_running = zeros(size(N_checkpoints));

cumsum_cross = cumsum(cross);   % cumsum：累积和（第 k 元素是前 k 个的和）
for idx = 1:length(N_checkpoints)
    n_k = N_checkpoints(idx);
    P_k = cumsum_cross(n_k) / n_k;   % 前 n_k 次的相交概率估计
    pi_k = 2*l / (d*P_k);
    pi_running(idx) = pi_k;
    err_running(idx) = abs(pi_k - pi);
end

%% ---- 第三步：多次重复实验（验证标准差）----
fprintf('--- 重复实验验证中心极限定理 ---\n');
n_repeat = 200;      % 重复实验次数
pi_estimates = zeros(1, n_repeat);

for rep = 1:n_repeat
    xc = (d/2) * rand(1, N);       % 每次独立生成随机针
    th = (pi/2) * rand(1, N);
    cr = (xc <= l/2 * sin(th));
    P_rep = sum(cr) / N;
    pi_estimates(rep) = 2*l / (d*P_rep);
end

pi_mean = mean(pi_estimates);    % 200次估计的均值
pi_std  = std(pi_estimates);     % 200次估计的标准差
% std(v)：向量 v 的样本标准差（分母为 N-1）

fprintf('%d 次重复实验（每次 N=%.0e 针）：\n', n_repeat, N);
fprintf('  均值    = %.8f（偏差 = %.4e）\n', pi_mean, abs(pi_mean-pi));
fprintf('  标准差  = %.6e\n', pi_std);
% 理论标准差：σ_P = sqrt(P(1-P)/N)，σ_π = 2l/(d·P²)·σ_P ≈ π/(2·P·N^0.5)...
fprintf('  理论上误差 ∝ 1/√N = 1/√%.0e = %.4e\n', N, 1/sqrt(N));
fprintf('\n');

%% ---- 第四步：绘图 ----
figure('Name', '蒲丰投针 Monte Carlo 估计 π', 'NumberTitle', 'off', ...
    'Position', [50, 50, 1200, 800]);

% 子图1：投针示意图（取前1000根针可视化）
subplot(2, 2, 1);
N_show = min(500, N);   % 最多显示500根（太多会卡）
x_show = x_center(1:N_show);   % 截取前 N_show 根针的数据
th_show = theta(1:N_show);

% 绘制平行线（y=0, y=d, y=2d, y=3d, y=4d）
for y_line = 0:d:4*d
    yline(y_line, 'k-', 'LineWidth', 1.5);   % yline：画水平线
end
hold on;
% 绘制每根针：中心 y 坐标均匀分布在 [0, 4d]
y_centers = 4*d * rand(1, N_show);   % 随机中心 y 坐标（仅用于可视化）
cross_show = x_show <= l/2 * sin(th_show);   % 是否相交

% 分别绘制相交（红色）和不相交（蓝色）的针
for i = 1:N_show
    % 针的两端坐标
    % x方向：从 (中心-l/2·cos(θ)) 到 (中心+l/2·cos(θ))
    % y方向：从 (中心-l/2·sin(θ)) 到 (中心+l/2·sin(θ))
    % 这里用 x 代表垂直方向（到直线的距离），y 代表水平方向
    y1 = y_centers(i) - l/2 * sin(th_show(i));
    y2 = y_centers(i) + l/2 * sin(th_show(i));
    x_shift = x_show(i);   % 中心到最近线的距离（不是实际坐标，仅示意）
    if cross_show(i)
        plot([y1, y2], [x_shift, x_shift], 'r-', 'LineWidth', 0.5);
    else
        plot([y1, y2], [x_shift, x_shift], 'b-', 'LineWidth', 0.5);
    end
end
hold off;
title(sprintf('投针示意（前%d根，红=相交，蓝=不相交）', N_show));
xlabel('针的位置（水平）'); ylabel('到最近直线的距离');
xlim([0, 4*d]); ylim([0, d/2]);

% 子图2：π 估计值随 N 收敛
subplot(2, 2, 2);
semilogx(N_checkpoints, pi_running, 'b-', 'LineWidth', 1.5);
hold on;
semilogx(N_checkpoints, pi*ones(size(N_checkpoints)), 'r--', 'LineWidth', 1.5);
hold off;
% semilogx：x轴对数刻度，y轴普通刻度
xlabel('投针次数 N（对数刻度）');
ylabel('π 估计值');
title('π 估计值随实验次数的收敛');
legend({'估计值', '真实值 π'}, 'Location', 'southeast');
grid on;

% 子图3：误差对数图（验证 O(1/√N) 收敛速度）
subplot(2, 2, 3);
loglog(N_checkpoints, err_running, 'b-', 'LineWidth', 1.5, 'DisplayName', '实际误差');
hold on;
loglog(N_checkpoints, 5./sqrt(N_checkpoints), 'r--', 'LineWidth', 1.5, 'DisplayName', 'O(1/√N)');
hold off;
xlabel('投针次数 N（对数刻度）');
ylabel('|π估计值 - π| （对数刻度）');
title('误差随 N 的收敛（斜率=-1/2 验证 O(1/√N)）');
legend;
grid on;
% 在双对数图上，O(1/√N) = O(N^{-1/2}) 是斜率为 -0.5 的直线

% 子图4：重复实验的分布直方图（验证中心极限定理）
subplot(2, 2, 4);
histogram(pi_estimates, 30, 'Normalization', 'pdf', 'FaceColor', [0.3 0.6 0.9]);
% histogram(data, bins, 'Normalization', 'pdf')：归一化为概率密度
hold on;
% 绘制理论正态分布曲线（中心极限定理：大 N 时近似正态）
x_norm = linspace(min(pi_estimates), max(pi_estimates), 200);
y_norm = normpdf(x_norm, pi_mean, pi_std);
% normpdf(x, mu, sigma)：正态分布概率密度函数
plot(x_norm, y_norm, 'r-', 'LineWidth', 2);
xline(pi, 'k--', 'LineWidth', 2);   % 真实 π 值的参考线
hold off;
xlabel('π 估计值');
ylabel('概率密度');
title(sprintf('重复%d次的分布（中心极限定理验证）\n均值=%.6f，标准差=%.4e', ...
    n_repeat, pi_mean, pi_std));
grid on;

sgtitle(sprintf('蒲丰投针 Monte Carlo 估计 π（l=%g，d=%g，N=%.0e）', l, d, N));
