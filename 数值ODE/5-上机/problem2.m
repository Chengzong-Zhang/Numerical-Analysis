%% =========================================================
%  problem2.m
%  第五章上机第2题
%  二阶ODE初值问题（化为一阶系统后用经典四阶RK求解）
%
%  方程：y'' + (2/x)*y' - (6/x^2)*y = 5 - 6x + 7x^2,  1 < x < 2
%  初值：y(1) = 1/2,  y'(1) = 2
%  真解：y = x^2 - x^3 + (1/2)*x^4 + x^2*ln(x)
%
%  步长：h = 1/200, 1/100, 1/50, 1/25, 1/10, 1/5, 1/4
%  任务：四阶RK计算各步长下的最大误差，分析收敛阶
% =========================================================
clear; close all; clc;

save_dir = fileparts(mfilename('fullpath'));

%% ---- 将二阶 ODE 化为一阶系统 ----
%  令 w(1) = y,  w(2) = y'
%  则 w'(1) = w(2)
%      w'(2) = -(2/x)*w(2) + (6/x^2)*w(1) + 5 - 6x + 7x^2
F = @(x, w) [w(2); -(2/x)*w(2) + (6/x^2)*w(1) + 5 - 6*x + 7*x^2];

%% ---- 精确解与参数 ----
y_exact = @(x) x.^2 - x.^3 + 0.5*x.^4 + x.^2 .* log(x);

x0 = 1;   xend = 2;
w0 = [0.5; 2];   % [y(1); y'(1)]

%% ---- 各步长下的误差计算 ----
h_list = [1/200, 1/100, 1/50, 1/25, 1/10, 1/5, 1/4];
nh     = length(h_list);
E_max  = zeros(1, nh);   % 最大绝对误差

for k = 1:nh
    hk    = h_list(k);
    Nk    = round((xend - x0) / hk);
    xk    = x0 : hk : xend;

    % ---- 经典四阶 RK（内联实现）----
    w = zeros(2, Nk+1);
    w(:,1) = w0;
    for n = 1:Nk
        xn = x0 + (n-1)*hk;
        K1 = F(xn,       w(:,n));
        K2 = F(xn+hk/2,  w(:,n) + hk/2*K1);
        K3 = F(xn+hk/2,  w(:,n) + hk/2*K2);
        K4 = F(xn+hk,    w(:,n) + hk*K3);
        w(:,n+1) = w(:,n) + hk/6*(K1 + 2*K2 + 2*K3 + K4);
    end

    y_num    = w(1,:);               % y 分量数值解
    y_ex     = y_exact(xk);
    E_max(k) = max(abs(y_num - y_ex));
end

%% ---- 收敛阶计算 ----
%  p_k ≈ log(E(h_{k-1})/E(h_k)) / log(h_{k-1}/h_k)
p_order    = nan(1, nh);
for k = 2:nh
    if E_max(k-1) > 0 && E_max(k) > 0
        p_order(k) = log(E_max(k-1)/E_max(k)) / log(h_list(k-1)/h_list(k));
    end
end

%% ---- 输出结果表 ----
fprintf('=============================================================\n');
fprintf('  第2题: 四阶RK求解二阶ODE —— 误差与收敛阶\n');
fprintf('  真解: y = x^2 - x^3 + (1/2)*x^4 + x^2*ln(x)\n');
fprintf('=============================================================\n');
fprintf('%-12s %-20s %-12s\n', 'h', '最大绝对误差 E(h)', '收敛阶 p');
fprintf('%s\n', repmat('-', 1, 48));
for k = 1:nh
    if isnan(p_order(k))
        fprintf('%-12s %-20.6e   ---\n', num2str(h_list(k),'%.5f'), E_max(k));
    else
        fprintf('%-12s %-20.6e   %-8.4f\n', num2str(h_list(k),'%.5f'), E_max(k), p_order(k));
    end
end
fprintf('\n理论收敛阶为 4（经典RK4 方法）。\n');

%% ---- 图1: 最小步长解曲线与真解对比 ----
h_best = h_list(1);   % h = 1/200
Nb     = round((xend - x0) / h_best);
x_best = x0 : h_best : xend;
wb     = zeros(2, Nb+1);
wb(:,1)= w0;
for n = 1:Nb
    xn = x0 + (n-1)*h_best;
    K1 = F(xn,           wb(:,n));
    K2 = F(xn+h_best/2,  wb(:,n) + h_best/2*K1);
    K3 = F(xn+h_best/2,  wb(:,n) + h_best/2*K2);
    K4 = F(xn+h_best,    wb(:,n) + h_best*K3);
    wb(:,n+1) = wb(:,n) + h_best/6*(K1 + 2*K2 + 2*K3 + K4);
end
y_exact_best = y_exact(x_best);

fig1 = figure('Name','p2_solution','NumberTitle','off','Position',[50 50 800 500]);
plot(x_best, y_exact_best,  'k-',  'LineWidth',2.0, 'DisplayName','真解');
hold on;
plot(x_best, wb(1,:),       'r--', 'LineWidth',1.5, 'DisplayName','RK4 (h=1/200)');
xlabel('x','FontSize',12); ylabel('y','FontSize',12);
title('第2题: 二阶ODE数值解与真解对比 (h=1/200)','FontSize',13);
legend('Location','best','FontSize',11); grid on;
exportgraphics(fig1, fullfile(save_dir, 'p2_solution.pdf'), 'ContentType', 'vector');
exportgraphics(fig1, fullfile(save_dir, 'p2_solution.png'), 'Resolution', 150);
fprintf('\n[已保存] p2_solution.pdf / .png\n');

%% ---- 图2: 误差 vs 步长（双对数图）----
fig2 = figure('Name','p2_convergence','NumberTitle','off','Position',[50 600 800 450]);
loglog(h_list, E_max, 'b-o', 'LineWidth',1.5,'MarkerSize',8,'DisplayName','RK4误差');
hold on;
h_ref = linspace(h_list(1), h_list(end), 50);
C_ref = E_max(1) / h_list(1)^4;
loglog(h_ref, C_ref * h_ref.^4, 'k--','LineWidth',1.2,'DisplayName','O(h^4)');
xlabel('步长 h','FontSize',12); ylabel('最大绝对误差 E(h)','FontSize',12);
title('第2题: 误差收敛阶（双对数图）','FontSize',13);
legend('Location','best','FontSize',11); grid on;
exportgraphics(fig2, fullfile(save_dir, 'p2_convergence.pdf'), 'ContentType', 'vector');
exportgraphics(fig2, fullfile(save_dir, 'p2_convergence.png'), 'Resolution', 150);
fprintf('[已保存] p2_convergence.pdf / .png\n');

%% ---- 图3: 不同步长的空间误差分布 ----
fig3 = figure('Name','p2_spatial_error','NumberTitle','off','Position',[900 50 800 500]);
h_show  = [1/200, 1/100, 1/50, 1/10, 1/4];
clr_arr = {'b','r','g','m','k'};
for ki = 1:length(h_show)
    hk = h_show(ki);
    Nk = round((xend - x0) / hk);
    xk = x0 : hk : xend;
    wk = zeros(2, Nk+1);  wk(:,1) = w0;
    for n = 1:Nk
        xn = x0 + (n-1)*hk;
        K1 = F(xn,       wk(:,n));
        K2 = F(xn+hk/2,  wk(:,n)+hk/2*K1);
        K3 = F(xn+hk/2,  wk(:,n)+hk/2*K2);
        K4 = F(xn+hk,    wk(:,n)+hk*K3);
        wk(:,n+1) = wk(:,n) + hk/6*(K1+2*K2+2*K3+K4);
    end
    semilogy(xk, abs(wk(1,:)-y_exact(xk)), [clr_arr{ki},'-o'], ...
        'LineWidth',1.2,'MarkerSize',4,'DisplayName',['h=',num2str(hk,'%.4f')]);
    hold on;
end
xlabel('x','FontSize',12); ylabel('|y_{num} - y_{exact}|','FontSize',12);
title('第2题: 不同步长的误差空间分布','FontSize',13);
legend('Location','best','FontSize',10); grid on;
exportgraphics(fig3, fullfile(save_dir, 'p2_spatial_error.pdf'), 'ContentType', 'vector');
exportgraphics(fig3, fullfile(save_dir, 'p2_spatial_error.png'), 'Resolution', 150);
fprintf('[已保存] p2_spatial_error.pdf / .png\n');

fprintf('\n第2题运行完毕。\n');
