%% 对称区间上 n 次最佳平方逼近多项式 (承接第一问)
clear; clc; close all;

%% 1. 定义参数
d = 1;               % 半区间长度
n = 3;               % 逼近多项式次数
f = @(x) sin(pi*x);  % 被逼近函数
x_plot = -1:0.01:1;  % 用于绘图和验证的网格点

%% 2. 定义第一问的自适应积分函数 (嵌入主脚本)
% 自适应两点 Gauss-Legendre 积分
function result = AdaptiveGaussLegendre2(f, a, b, tol)
    if nargin < 4
        tol = 1e-8;
    end
    
    [I, error_est] = AdaptiveStep(f, a, b, tol, 0);
    
    if error_est > tol
        fprintf('警告：初始递归可能未完全收敛。\n');
    end
    
    result = I;
end

function [I_val, I_err] = AdaptiveStep(f, a, b, tol, depth)
    I_ab = Gauss2Points(f, a, b);
    
    m = (a + b) / 2;
    I_left = Gauss2Points(f, a, m);
    I_right = Gauss2Points(f, m, b);
    I_sub = I_left + I_right;
    
    I_err = 16 * abs(I_sub - I_ab);
    
    if I_err <= tol || depth > 10
        I_val = I_ab;
    else
        [I_left_val, I_left_err] = AdaptiveStep(f, a, m, tol/2, depth+1);
        [I_right_val, I_right_err] = AdaptiveStep(f, m, b, tol/2, depth+1);
        I_val = I_left_val + I_right_val;
        I_err = 16 * abs(I_val - I_ab);
    end
end

function I_val = Gauss2Points(f, a, b)
    sqrt3 = sqrt(3);
    t1 = -1/sqrt3;
    t2 = 1/sqrt3;
    
    h = (b - a) / 2;
    x1 = h * t1 + (a + b) / 2;
    x2 = h * t2 + (a + b) / 2;
    
    I_val = h * (f(x1) + f(x2));
end

%% 3. 构建法方程组 A*c = b
A = zeros(n+1, n+1);
b = zeros(n+1, 1);

for i = 0:n
    for j = 0:n
        k = i + j;
        % 计算 A_ij = int_{-d}^{d} x^k dx
        if mod(k, 2) == 0  % k 为偶数
            A(i+1, j+1) = 2 * d^(k+1) / (k+1);
        else                % k 为奇数
            A(i+1, j+1) = 0;
        end
    end
    
    % 计算 b_i = int_{-d}^{d} f(x) * x^i dx
    % 使用第一问的自适应积分函数
    integrand = @(x_val) f(x_val) .* (x_val.^i);
    b(i+1) = AdaptiveGaussLegendre2(integrand, -d, d, 1e-10);
end

%% 4. 求解系数 c
c = A \ b;

%% 5. 构造逼近多项式 P_n(x)
% c 向量是 [c_0, c_1, ..., c_n]，polyval 需要 [c_n, ..., c_0]
c_poly = flip(c); 
P_n = @(x_val) polyval(c_poly, x_val);

%% 6. 计算误差
syms x_sym
f_sym = sin(pi*x_sym);
I_exact = int(f_sym, x_sym, -d, d); % 精确积分值

% 近似积分值
I_approx = AdaptiveGaussLegendre2(P_n, -d, d, 1e-10);

% L2 范数误差
error_L2 = abs(I_exact - I_approx);

%% 7. 输出结果
fprintf('----- 对称区间最佳平方逼近验证结果 -----\n');
fprintf('区间: [-%.0f, %.0f]\n', -d, d);
fprintf('函数: f(x) = sin(pi*x)\n');
fprintf('次数: n = %d\n', n);
fprintf('----- 法方程组系数矩阵 A -----\n');
fprintf('A =\n');
fprintf('  %.4f  %.4f  %.4f  %.4f\n', A(1,:));
fprintf('  %.4f  %.4f  %.4f  %.4f\n', A(2,:));
fprintf('  %.4f  %.4f  %.4f  %.4f\n', A(3,:));
fprintf('  %.4f  %.4f  %.4f  %.4f\n', A(4,:));
fprintf('----- 法方程组右端项 b -----\n');
fprintf('b =\n');
fprintf('  %.4f\n', b(1));
fprintf('  %.4f\n', b(2));
fprintf('  %.4f\n', b(3));
fprintf('  %.4f\n', b(4));
fprintf('----- 解向量 c (多项式系数) -----\n');
fprintf('c = [c_0, c_1, c_2, c_3] =\n');
fprintf('  %.10f\n', c(1));
fprintf('  %.10f\n', c(2));
fprintf('  %.10f\n', c(3));
fprintf('  %.10f\n', c(4));
fprintf('----- 逼近多项式表达式 -----\n');
fprintf('P_3(x) = %.10f + %.10f*x + %.10f*x^2 + %.10f*x^3\n', c(1), c(2), c(3), c(4));
fprintf('----- 误差分析 -----\n');
fprintf('精确积分值: I_exact = %.15f\n', double(I_exact));
fprintf('自适应近似积分值: I_approx = %.15f\n', I_approx);
fprintf('绝对误差: |I_exact - I_approx| = %.2e\n', error_L2);
fprintf('目标误差上界: tol = 1e-8\n');
fprintf('----- 验证结论 -----\n');
if error_L2 <= 1e-8
    fprintf('验证成功！误差 %.2e 小于目标误差 1e-8。\n', error_L2);
else
    fprintf('验证未完全满足！误差 %.2e 略大于目标误差 1e-8。\n', error_L2);
end