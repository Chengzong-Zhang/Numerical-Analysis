%% =========================================================
%  problem1.m
%  第五章上机第1题
%  ODE初值问题：y' = -1/x^2 - y/x - y^2, 1<=x<=2, y(1)=-1
%
%  第(1)部分：Euler / 改进Euler / Heun / 中点 / 四阶RK 五种方法比较
%  第(2)部分：四阶Adams PECE 与 RK4 比较
% =========================================================
clear; close all; clc;

%% ---- 基本参数 ----
f   = @(x, y) -1./x.^2 - y./x - y.^2;   % ODE 右端函数
x0  = 1;   xend = 2;   y0 = -1;
h   = 0.1;
N   = round((xend - x0) / h);
x   = x0 : h : xend;                     % 1×(N+1) 网格

save_dir = fileparts(mfilename('fullpath'));   % 保存到脚本所在目录

%% =========================================================
%  第(1)部分：五种单步方法
% =========================================================

%% --- 1. 向前 Euler 方法 ---
y_eu = zeros(1, N+1);  y_eu(1) = y0;
for n = 1:N
    y_eu(n+1) = y_eu(n) + h * f(x(n), y_eu(n));
end

%% --- 2. 改进的 Euler 方法（显式梯形 / 2阶RK a2=1）---
% y_{n+1} = y_n + h/2 * [f(x_n,y_n) + f(x_{n+1}, y_n+h*f_n)]
y_me = zeros(1, N+1);  y_me(1) = y0;
for n = 1:N
    k1 = f(x(n),   y_me(n));
    k2 = f(x(n+1), y_me(n) + h*k1);
    y_me(n+1) = y_me(n) + h/2 * (k1 + k2);
end

%% --- 3. Heun 公式（2阶RK a2=2/3）---
% K1 = f(x_n, y_n)
% K2 = f(x_n + 2h/3, y_n + 2h/3 * K1)
% y_{n+1} = y_n + h/4 * (K1 + 3*K2)
y_he = zeros(1, N+1);  y_he(1) = y0;
for n = 1:N
    K1 = f(x(n),         y_he(n));
    K2 = f(x(n)+2*h/3,   y_he(n) + 2*h/3 * K1);
    y_he(n+1) = y_he(n) + h/4 * (K1 + 3*K2);
end

%% --- 4. 中点方法（变形Euler / 2阶RK a2=1/2）---
% y_{n+1} = y_n + h * f(x_n+h/2, y_n + h/2*f(x_n,y_n))
y_mp = zeros(1, N+1);  y_mp(1) = y0;
for n = 1:N
    K1 = f(x(n), y_mp(n));
    y_mp(n+1) = y_mp(n) + h * f(x(n)+h/2, y_mp(n)+h/2*K1);
end

%% --- 5. 经典四阶 Runge-Kutta (RK4) ---
y_rk4 = zeros(1, N+1);  y_rk4(1) = y0;
for n = 1:N
    K1 = f(x(n),     y_rk4(n));
    K2 = f(x(n)+h/2, y_rk4(n)+h/2*K1);
    K3 = f(x(n)+h/2, y_rk4(n)+h/2*K2);
    K4 = f(x(n)+h,   y_rk4(n)+h*K3);
    y_rk4(n+1) = y_rk4(n) + h/6*(K1 + 2*K2 + 2*K3 + K4);
end

%% ---- 输出数值结果表 ----
fprintf('=============================================================\n');
fprintf('  第1题(1): 五种方法数值解对比  h = %.1f\n', h);
fprintf('=============================================================\n');
fprintf('%-6s %-12s %-12s %-12s %-12s %-12s\n', ...
    'x', 'Euler', '改进Euler', 'Heun', '中点', 'RK4');
fprintf('%s\n', repmat('-', 1, 68));
for n = 1:N+1
    fprintf('%-6.2f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f\n', ...
        x(n), y_eu(n), y_me(n), y_he(n), y_mp(n), y_rk4(n));
end

% 以 RK4 作参考，输出各方法与 RK4 的差
fprintf('\n各方法与 RK4 的差（绝对误差）：\n');
fprintf('%-6s %-12s %-12s %-12s %-12s\n', 'x', 'Euler', '改进Euler', 'Heun', '中点');
fprintf('%s\n', repmat('-', 1, 55));
for n = 1:N+1
    fprintf('%-6.2f %-12.2e %-12.2e %-12.2e %-12.2e\n', x(n), ...
        abs(y_eu(n)-y_rk4(n)), abs(y_me(n)-y_rk4(n)), ...
        abs(y_he(n)-y_rk4(n)), abs(y_mp(n)-y_rk4(n)));
end

%% ---- 图1: 五种方法解曲线对比 ----
fig1 = figure('Name','p1_1_solutions','NumberTitle','off','Position',[50 50 800 500]);
plot(x, y_eu, 'b-o',  'LineWidth',1.2, 'MarkerSize',5, 'DisplayName','Euler');
hold on;
plot(x, y_me, 'r-s',  'LineWidth',1.2, 'MarkerSize',5, 'DisplayName','改进Euler');
plot(x, y_he, 'g-d',  'LineWidth',1.2, 'MarkerSize',5, 'DisplayName','Heun');
plot(x, y_mp, 'm-^',  'LineWidth',1.2, 'MarkerSize',5, 'DisplayName','中点');
plot(x, y_rk4,'k-*',  'LineWidth',1.5, 'MarkerSize',6, 'DisplayName','RK4');
xlabel('x','FontSize',12); ylabel('y','FontSize',12);
title('第1题(1): 五种数值方法解曲线对比 (h=0.1)','FontSize',13);
legend('Location','best','FontSize',11);
grid on;
exportgraphics(fig1, fullfile(save_dir, 'p1_part1_solutions.pdf'), 'ContentType', 'vector');
exportgraphics(fig1, fullfile(save_dir, 'p1_part1_solutions.png'), 'Resolution', 150);
fprintf('\n[已保存] p1_part1_solutions.pdf / .png\n');

%% ---- 图2: 各方法与 RK4 的差 ----
fig2 = figure('Name','p1_1_errors','NumberTitle','off','Position',[50 600 800 400]);
semilogy(x(2:end), abs(y_eu(2:end)-y_rk4(2:end)),  'b-o','LineWidth',1.2,'MarkerSize',5,'DisplayName','Euler');
hold on;
semilogy(x(2:end), abs(y_me(2:end)-y_rk4(2:end)),  'r-s','LineWidth',1.2,'MarkerSize',5,'DisplayName','改进Euler');
semilogy(x(2:end), abs(y_he(2:end)-y_rk4(2:end)),  'g-d','LineWidth',1.2,'MarkerSize',5,'DisplayName','Heun');
semilogy(x(2:end), abs(y_mp(2:end)-y_rk4(2:end)),  'm-^','LineWidth',1.2,'MarkerSize',5,'DisplayName','中点');
xlabel('x','FontSize',12); ylabel('|y_{method} - y_{RK4}|','FontSize',12);
title('第1题(1): 各方法与RK4解之差 (h=0.1)','FontSize',13);
legend('Location','best','FontSize',11);
grid on;
exportgraphics(fig2, fullfile(save_dir, 'p1_part1_errors.pdf'), 'ContentType', 'vector');
exportgraphics(fig2, fullfile(save_dir, 'p1_part1_errors.png'), 'Resolution', 150);
fprintf('[已保存] p1_part1_errors.pdf / .png\n');

%% =========================================================
%  第(2)部分：四阶Adams PECE（用RK4提供起步值）与 RK4 比较
% =========================================================
%
%  四步显式Adams-Bashforth (预报):
%    y^P_{n+1} = y_n + h/24*(55f_n - 59f_{n-1} + 37f_{n-2} - 9f_{n-3})
%  三步隐式Adams-Moulton (校正):
%    y_{n+1}   = y_n + h/24*(9f^P_{n+1} + 19f_n - 5f_{n-1} + f_{n-2})
% =========================================================

y_ad = zeros(1, N+1);
fv   = zeros(1, N+1);   % 存储 f(x_n, y_n) 以供多步法使用

% 用 RK4 计算前 4 个起步值 y_0, y_1, y_2, y_3
y_ad(1:4) = y_rk4(1:4);
for i = 1:4
    fv(i) = f(x(i), y_ad(i));
end

% Adams PECE 主循环（从第4步起）
for n = 4:N
    % P: 预报
    y_p = y_ad(n) + h/24 * (55*fv(n) - 59*fv(n-1) + 37*fv(n-2) - 9*fv(n-3));
    % E: 计算预报值的函数值
    f_p = f(x(n+1), y_p);
    % C: 校正
    y_ad(n+1) = y_ad(n) + h/24 * (9*f_p + 19*fv(n) - 5*fv(n-1) + fv(n-2));
    % E: 更新函数值
    fv(n+1) = f(x(n+1), y_ad(n+1));
end

%% ---- 输出对比表 ----
fprintf('\n=============================================================\n');
fprintf('  第1题(2): Adams PECE 与 RK4 比较  h = %.1f\n', h);
fprintf('=============================================================\n');
fprintf('%-6s %-14s %-14s %-14s\n', 'x', 'RK4', 'Adams PECE', '|差值|');
fprintf('%s\n', repmat('-', 1, 52));
for n = 1:N+1
    fprintf('%-6.2f %-14.8f %-14.8f %-14.2e\n', ...
        x(n), y_rk4(n), y_ad(n), abs(y_rk4(n)-y_ad(n)));
end

%% ---- 图3: Adams PECE 与 RK4 对比 ----
fig3 = figure('Name','p1_2_adams','NumberTitle','off','Position',[900 50 800 500]);
plot(x, y_rk4, 'b-o', 'LineWidth',1.5,'MarkerSize',6,'DisplayName','RK4');
hold on;
plot(x, y_ad,  'r-s', 'LineWidth',1.5,'MarkerSize',6,'DisplayName','Adams PECE');
xlabel('x','FontSize',12); ylabel('y','FontSize',12);
title('第1题(2): Adams PECE 与 RK4 对比 (h=0.1)','FontSize',13);
legend('Location','best','FontSize',11);
grid on;
exportgraphics(fig3, fullfile(save_dir, 'p1_part2_adams.pdf'), 'ContentType', 'vector');
exportgraphics(fig3, fullfile(save_dir, 'p1_part2_adams.png'), 'Resolution', 150);
fprintf('[已保存] p1_part2_adams.pdf / .png\n');

%% ---- 图4: Adams PECE 与 RK4 差值（从第4步起） ----
fig4 = figure('Name','p1_2_adams_diff','NumberTitle','off','Position',[900 600 800 400]);
idx_start = 5;   % Adams 从第4步(n=4)才开始独立计算
semilogy(x(idx_start:end), abs(y_ad(idx_start:end)-y_rk4(idx_start:end)), ...
    'k-d','LineWidth',1.5,'MarkerSize',6);
xlabel('x','FontSize',12); ylabel('|y_{Adams} - y_{RK4}|','FontSize',12);
title('第1题(2): Adams PECE 与 RK4 差值 (h=0.1)','FontSize',13);
grid on;
exportgraphics(fig4, fullfile(save_dir, 'p1_part2_diff.pdf'), 'ContentType', 'vector');
exportgraphics(fig4, fullfile(save_dir, 'p1_part2_diff.png'), 'Resolution', 150);
fprintf('[已保存] p1_part2_diff.pdf / .png\n');

fprintf('\n第1题运行完毕。\n');
