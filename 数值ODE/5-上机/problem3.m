%% =========================================================
%  problem3.m
%  第五章上机第3题
%  刚性ODE方程组初值问题
%
%  u' = 32u + 66v + (2/3)(x+1)
%  v' = -66u - 133v - (1/3)(x+1),   x in [0, 0.5]
%  u(0) = v(0) = 1/3
%
%  真解：
%    u = (2/3)x + (2/3)exp(-x) - (1/3)exp(-100x)
%    v = -(1/3)x - (1/3)exp(-x) + (2/3)exp(-100x)
%
%  方法：Euler / 梯形（隐式）/ 经典四阶RK
%  步长：h = 1e-5, 1e-4, 1e-3, 1e-2, 1e-1
%  任务：误差收敛阶 + 计算时间对比 + 刚性分析
% =========================================================
clear; close all; clc;

save_dir = fileparts(mfilename('fullpath'));

%% ---- 系统矩阵与精确解 ----
A = [32, 66; -66, -133];                          % 系数矩阵
g = @(x) [(2/3)*(x+1); -(1/3)*(x+1)];            % 非齐次项
F = @(x, y) A*y + g(x);                           % 右端函数 f(x,y)

u_exact = @(x) (2/3)*x + (2/3)*exp(-x) - (1/3)*exp(-100*x);
v_exact = @(x) -(1/3)*x - (1/3)*exp(-x) + (2/3)*exp(-100*x);

x0 = 0;  xend = 0.5;
y0 = [1/3; 1/3];

%% ---- 说明刚性 ----
% 矩阵 A 的特征值
eigvals = eig(A);
fprintf('=============================================================\n');
fprintf('  第3题: 刚性ODE组 —— Euler / 梯形 / RK4 比较\n');
fprintf('=============================================================\n');
fprintf('系数矩阵 A 的特征值: lambda_1 = %.1f, lambda_2 = %.1f\n', eigvals(1), eigvals(2));
fprintf('刚性比 (stiffness ratio) = %.0f\n', max(abs(eigvals))/min(abs(eigvals)));
fprintf('Euler 绝对稳定要求: h < %.4f\n', 2/max(abs(eigvals)));
fprintf('RK4   绝对稳定要求: h < %.4f\n', 2.785/max(abs(eigvals)));
fprintf('梯形方法 A-stable: 任意步长均稳定\n\n');

%% ---- 梯形法预计算矩阵（线性系统可精确求解）----
% (I - h/2*A) * y_{n+1} = (I + h/2*A) * y_n + h/2*(g(x_n) + g(x_{n+1}))
% 对每个不同的 h 需重新构造
I2 = eye(2);

%% ---- 步长列表 ----
h_list = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1];
nh     = length(h_list);

% 存储每种方法在每个步长下的：误差、运行时间
E_eu   = zeros(1, nh);   t_eu  = zeros(1, nh);
E_tr   = zeros(1, nh);   t_tr  = zeros(1, nh);
E_rk4  = zeros(1, nh);   t_rk4 = zeros(1, nh);

% 收敛阶
p_eu   = nan(1, nh);
p_tr   = nan(1, nh);
p_rk4  = nan(1, nh);

for k = 1:nh
    hk = h_list(k);
    Nk = round((xend - x0) / hk);

    %% ---- Euler 方法 ----
    tic;
    y_eu = zeros(2, Nk+1);
    y_eu(:,1) = y0;
    xn = x0;
    for n = 1:Nk
        y_eu(:,n+1) = y_eu(:,n) + hk * F(xn, y_eu(:,n));
        xn = xn + hk;
    end
    t_eu(k) = toc;
    xend_val = x0 + Nk*hk;
    eu_u = y_eu(1,end);  eu_v = y_eu(2,end);
    E_eu(k) = max(abs(eu_u - u_exact(xend_val)), abs(eu_v - v_exact(xend_val)));

    %% ---- 梯形方法（隐式，线性系统直接求解）----
    tic;
    M_minus = I2 - hk/2 * A;    % (I - h/2*A)，不依赖 x 和 y（线性）
    M_plus  = I2 + hk/2 * A;    % (I + h/2*A)
    y_tr = zeros(2, Nk+1);
    y_tr(:,1) = y0;
    xn = x0;
    for n = 1:Nk
        xn1 = xn + hk;
        rhs = M_plus * y_tr(:,n) + hk/2 * (g(xn) + g(xn1));
        y_tr(:,n+1) = M_minus \ rhs;
        xn = xn1;
    end
    t_tr(k) = toc;
    tr_u = y_tr(1,end);  tr_v = y_tr(2,end);
    E_tr(k) = max(abs(tr_u - u_exact(xend_val)), abs(tr_v - v_exact(xend_val)));

    %% ---- 经典四阶 RK ----
    tic;
    y_rk = zeros(2, Nk+1);
    y_rk(:,1) = y0;
    xn = x0;
    for n = 1:Nk
        w  = y_rk(:,n);
        K1 = F(xn,       w);
        K2 = F(xn+hk/2,  w + hk/2*K1);
        K3 = F(xn+hk/2,  w + hk/2*K2);
        K4 = F(xn+hk,    w + hk*K3);
        y_rk(:,n+1) = w + hk/6*(K1 + 2*K2 + 2*K3 + K4);
        xn = xn + hk;
    end
    t_rk4(k) = toc;
    rk_u = y_rk(1,end);  rk_v = y_rk(2,end);
    E_rk4(k) = max(abs(rk_u - u_exact(xend_val)), abs(rk_v - v_exact(xend_val)));
end

%% ---- 收敛阶计算 ----
for k = 2:nh
    if E_eu(k-1) > 0 && E_eu(k) > 0
        p_eu(k)  = log(E_eu(k-1)/E_eu(k)) / log(h_list(k-1)/h_list(k));
    end
    if E_tr(k-1) > 0 && E_tr(k) > 0
        p_tr(k)  = log(E_tr(k-1)/E_tr(k)) / log(h_list(k-1)/h_list(k));
    end
    if E_rk4(k-1) > 0 && E_rk4(k) > 0
        p_rk4(k) = log(E_rk4(k-1)/E_rk4(k)) / log(h_list(k-1)/h_list(k));
    end
end

%% ---- 输出综合结果表 ----
fprintf('误差与收敛阶（在 x=0.5 处，max(|u_err|, |v_err|)）：\n');
fprintf('%s\n', repmat('=', 1, 95));
fprintf('%-8s | %-14s %-6s %-10s | %-14s %-6s %-10s | %-14s %-6s %-10s\n', ...
    'h', '|Euler误差|', '阶', '时间(s)', '|梯形误差|', '阶', '时间(s)', '|RK4误差|', '阶', '时间(s)');
fprintf('%s\n', repmat('-', 1, 95));
for k = 1:nh
    p_eu_str  = fmt_order(p_eu(k));
    p_tr_str  = fmt_order(p_tr(k));
    p_rk4_str = fmt_order(p_rk4(k));
    fprintf('%-8.0e | %-14.4e %-6s %-10.4f | %-14.4e %-6s %-10.4f | %-14.4e %-6s %-10.4f\n', ...
        h_list(k), ...
        E_eu(k),  p_eu_str,  t_eu(k), ...
        E_tr(k),  p_tr_str,  t_tr(k), ...
        E_rk4(k), p_rk4_str, t_rk4(k));
end
fprintf('%s\n', repmat('=', 1, 95));
fprintf('\n注：h=1e-1 时 Euler 和 RK4 超出绝对稳定区间，误差可能极大或发散。\n');
fprintf('梯形方法 A-stable，对所有步长均稳定收敛。\n\n');

%% ---- 图1: 误差 vs 步长（双对数）----
fig1 = figure('Name','p3_convergence','NumberTitle','off','Position',[50 50 900 500]);
% 只显示收敛范围内的数据（剔除发散结果，仅用 E < 1e6）
eu_valid  = E_eu  < 1e6;
rk4_valid = E_rk4 < 1e6;

loglog(h_list(eu_valid),  E_eu(eu_valid),  'b-o', 'LineWidth',1.5,'MarkerSize',7,'DisplayName','Euler');
hold on;
loglog(h_list,            E_tr,            'r-s', 'LineWidth',1.5,'MarkerSize',7,'DisplayName','梯形');
loglog(h_list(rk4_valid), E_rk4(rk4_valid),'g-d', 'LineWidth',1.5,'MarkerSize',7,'DisplayName','RK4');

% 参考线
h_ref = h_list(2:end);
C1  = E_eu(2)  / h_list(2)^1;
C2  = E_tr(2)  / h_list(2)^2;
C4  = E_rk4(2) / h_list(2)^4;
loglog(h_ref, C1*h_ref.^1, 'b--', 'LineWidth',1.0,'DisplayName','O(h^1)');
loglog(h_ref, C2*h_ref.^2, 'r--', 'LineWidth',1.0,'DisplayName','O(h^2)');
loglog(h_ref, C4*h_ref.^4, 'g--', 'LineWidth',1.0,'DisplayName','O(h^4)');

xlabel('步长 h','FontSize',12); ylabel('最大误差 E(h)','FontSize',12);
title('第3题: 三种方法误差收敛阶（双对数图）','FontSize',13);
legend('Location','best','FontSize',10);
grid on;
exportgraphics(fig1, fullfile(save_dir, 'p3_convergence.pdf'), 'ContentType', 'vector');
exportgraphics(fig1, fullfile(save_dir, 'p3_convergence.png'), 'Resolution', 150);
fprintf('[已保存] p3_convergence.pdf / .png\n');

%% ---- 图2: 计算时间对比（条形图）----
fig2 = figure('Name','p3_timing','NumberTitle','off','Position',[50 600 900 400]);
h_labels = {'1e-5','1e-4','1e-3','1e-2','1e-1'};
bar_data  = [t_eu; t_tr; t_rk4]';
bar(bar_data);
set(gca,'XTickLabel', h_labels,'FontSize',11);
xlabel('步长 h','FontSize',12); ylabel('计算时间 (秒)','FontSize',12);
title('第3题: 三种方法计算时间对比','FontSize',13);
legend({'Euler','梯形','RK4'},'Location','best','FontSize',11);
grid on;
exportgraphics(fig2, fullfile(save_dir, 'p3_timing.pdf'), 'ContentType', 'vector');
exportgraphics(fig2, fullfile(save_dir, 'p3_timing.png'), 'Resolution', 150);
fprintf('[已保存] p3_timing.pdf / .png\n');

%% ---- 图3: 稳定步长下的解曲线（h=1e-3）----
h_demo = 1e-3;
N_demo = round((xend - x0) / h_demo);
x_demo = x0 : h_demo : xend;

% 三种方法
y_eu_demo  = zeros(2, N_demo+1);  y_eu_demo(:,1)  = y0;
y_tr_demo  = zeros(2, N_demo+1);  y_tr_demo(:,1)  = y0;
y_rk_demo  = zeros(2, N_demo+1);  y_rk_demo(:,1)  = y0;
M_minus_d  = I2 - h_demo/2 * A;
M_plus_d   = I2 + h_demo/2 * A;

for n = 1:N_demo
    xn  = x_demo(n);   xn1 = x_demo(n+1);
    % Euler
    y_eu_demo(:,n+1)  = y_eu_demo(:,n) + h_demo * F(xn, y_eu_demo(:,n));
    % 梯形
    rhs = M_plus_d * y_tr_demo(:,n) + h_demo/2*(g(xn)+g(xn1));
    y_tr_demo(:,n+1)  = M_minus_d \ rhs;
    % RK4
    w  = y_rk_demo(:,n);
    K1 = F(xn,          w);
    K2 = F(xn+h_demo/2, w+h_demo/2*K1);
    K3 = F(xn+h_demo/2, w+h_demo/2*K2);
    K4 = F(xn+h_demo,   w+h_demo*K3);
    y_rk_demo(:,n+1)  = w + h_demo/6*(K1+2*K2+2*K3+K4);
end

% 真解
u_ex = u_exact(x_demo);
v_ex = v_exact(x_demo);

fig3 = figure('Name','p3_solution_h1e3','NumberTitle','off','Position',[900 50 900 600]);
subplot(2,1,1);
plot(x_demo, u_ex,              'k-',  'LineWidth',2.0, 'DisplayName','u 真解');
hold on;
plot(x_demo, y_eu_demo(1,:),    'b--', 'LineWidth',1.2, 'DisplayName','u Euler');
plot(x_demo, y_tr_demo(1,:),    'r--', 'LineWidth',1.2, 'DisplayName','u 梯形');
plot(x_demo, y_rk_demo(1,:),    'g--', 'LineWidth',1.2, 'DisplayName','u RK4');
ylabel('u','FontSize',12); title(['第3题: 解曲线对比 (h=',num2str(h_demo),')'],'FontSize',13);
legend('Location','best','FontSize',10); grid on;

subplot(2,1,2);
plot(x_demo, v_ex,              'k-',  'LineWidth',2.0, 'DisplayName','v 真解');
hold on;
plot(x_demo, y_eu_demo(2,:),    'b--', 'LineWidth',1.2, 'DisplayName','v Euler');
plot(x_demo, y_tr_demo(2,:),    'r--', 'LineWidth',1.2, 'DisplayName','v 梯形');
plot(x_demo, y_rk_demo(2,:),    'g--', 'LineWidth',1.2, 'DisplayName','v RK4');
xlabel('x','FontSize',12); ylabel('v','FontSize',12);
legend('Location','best','FontSize',10); grid on;
exportgraphics(fig3, fullfile(save_dir, 'p3_solution.pdf'), 'ContentType', 'vector');
exportgraphics(fig3, fullfile(save_dir, 'p3_solution.png'), 'Resolution', 150);
fprintf('[已保存] p3_solution.pdf / .png\n');

fprintf('\n第3题运行完毕。\n');

%% ---- 辅助函数：格式化收敛阶输出 ----
function s = fmt_order(p)
    if isnan(p)
        s = '---';
    elseif p > 100 || p < -100
        s = 'N/A';    % 发散或精度到头
    else
        s = sprintf('%.2f', p);
    end
end
