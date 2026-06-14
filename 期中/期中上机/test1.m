function result = AdaptiveGaussLegendre2(f, a, b, tol)
% AdaptiveGaussLegendre2: 自适应两点 Gauss-Legendre 积分
%   result = AdaptiveGaussLegendre2(f, a, b, tol) 计算积分 I = int_a^b f(x) dx
%   要求总误差上界为 tol (默认 1e-8)
%
%   输入:
%     f  - 函数句柄
%     a, b - 积分下限和上限
%     tol - 目标误差上界 (默认 1e-8)
%
%   输出:
%     result - 积分近似值

    if nargin < 4
        tol = 1e-8; % 默认误差上界
    end
    
    % 初始化全局变量用于存储精确参考值 (可选，这里采用纯自适应策略)
    % 自适应策略：比较 [a,b] 的整体估计与 [a,m] + [m,b] 的估计
    
    [I, error_est] = AdaptiveStep(f, a, b, tol, 0);
    
    if error_est > tol
        fprintf('警告：初始递归可能未完全收敛，建议增加 tol 精度或检查函数光滑性。\n');
    end
    
    result = I;
end

function [I_val, I_err] = AdaptiveStep(f, a, b, tol, depth)
    % 递归计算积分
    % 如果深度过深或区间足够小，直接计算
    
    % 1. 计算整体区间 [a, b] 的两点 Gauss 近似值
    I_ab = Gauss2Points(f, a, b);
    
    % 2. 计算两个子区间 [a, m] 和 [m, b] 的两点 Gauss 近似值
    m = (a + b) / 2;
    I_left = Gauss2Points(f, a, m);
    I_right = Gauss2Points(f, m, b);
    I_sub = I_left + I_right;
    
    % 3. 误差估计 (Richardson 外推思想)
    % 对于两点 Gauss，整体误差 E ~ C * h^4 (h 为区间长度)
    % 子区间误差 E_sub ~ C * (h/2)^4 = C * h^4 / 16
    % 所以 E ~ 16 * (I_sub - I_ab)
    % 这是一个常用的误差估计公式，适用于高斯求积
    I_err = 16 * abs(I_sub - I_ab);
    
    % 4. 递归判断
    if I_err <= tol || depth > 10 % 防止无限递归
        I_val = I_ab;
    else
        % 递归细分
        [I_left_val, I_left_err] = AdaptiveStep(f, a, m, tol/2, depth+1); % 放宽子区间精度
        [I_right_val, I_right_err] = AdaptiveStep(f, m, b, tol/2, depth+1);
        I_val = I_left_val + I_right_val;
        % 重新估计总误差
        I_err = 16 * abs(I_val - I_ab);
    end
end

function I_val = Gauss2Points(f, a, b)
    % 两点 Gauss-Legendre 积分公式
    % 节点在 [-1, 1] 上: -1/sqrt(3), 1/sqrt(3)
    % 权重均为 1
    sqrt3 = sqrt(3);
    t1 = -1/sqrt3;
    t2 = 1/sqrt3;
    
    % 变换到 [a, b]
    h = (b - a) / 2;
    x1 = h * t1 + (a + b) / 2;
    x2 = h * t2 + (a + b) / 2;
    
    % 积分值
    I_val = h * (f(x1) + f(x2));
end

% ==================== 主程序：验证算法 ====================
% 定义被积函数
f = @(x) 0.5 * (5*x.^3 - 3*x) .* sin(pi*x);

% 积分区间 [-1, 1]
a = -1;
b = 1;

% 目标误差
tol = 1e-8;

% 调用自适应积分
I_approx = AdaptiveGaussLegendre2(f, a, b, tol);

% 计算精确解 (符号积分)
syms x_sym
f_sym = 0.5 * (5*x_sym^3 - 3*x_sym) * sin(pi*x_sym);
I_exact_sym = int(f_sym, x_sym, a, b);
I_exact = double(I_exact_sym);

% 输出结果
fprintf('----- 自适应两点 Gauss-Legendre 积分验证结果 -----\n');
fprintf('积分区间: [%.0f, %.0f]\n', a, b);
fprintf('被积函数: f(x) = 0.5*(5*x^3 - 3*x)*sin(pi*x)\n');
fprintf('精确解 (符号计算): I_exact = %.15f\n', I_exact);
fprintf('自适应近似解:       I_approx = %.15f\n', I_approx);
fprintf('----- 误差分析 -----\n');
fprintf('绝对误差: |I_exact - I_approx| = %.2e\n', abs(I_exact - I_approx));
fprintf('相对误差: |I_exact - I_approx| / |I_exact| = %.2e\n', abs(I_exact - I_approx) /abs(I_exact));
fprintf('目标误差上界:       tol = %.1e\n', tol);
fprintf('----- 验证结论 -----\n');
if abs(I_exact - I_approx) <= tol
    fprintf('验证成功！误差 %.2e 小于目标误差 %.1e。\n', abs(I_exact - I_approx), tol);
else
    fprintf('验证未完全满足！误差 %.2e 略大于目标误差 %.1e。\n', abs(I_exact - I_approx), tol);
end