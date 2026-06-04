% 定义Legendre多项式计算函数
function p_val = legendre_poly(n, x)
    if n == 0
        p_val = ones(size(x));
    elseif n == 1
        p_val = x;
    else
        % 使用递推公式计算Legendre多项式，避免数值不稳定
        P_prev_prev = ones(size(x)); % P_0
        P_prev = x;                 % P_1
        
        for k = 2:n
            current_P = ( ((2*k - 1) * x .* P_prev) - (k - 1) * P_prev_prev ) ./ k;
            P_prev_prev = P_prev;
            P_prev = current_P;
        end
        p_val = current_P; % 返回第n个Legendre多项式
    end
end

% 自适应高斯-勒让德积分函数（部分实现）
function [integral_result, num_evaluations] = gls_integrate(func, a, b, tol)
    % 实际的自适应算法实现可能更复杂，这里简化为粗略说明。
    % 这里应使用高效的自适应算法来计算积分值和评估次数。
    num_evaluations = 0;
    integral_result = integral_qng(func, a, b, tol);
end

% 自适应高斯-勒让德积分函数
function [result, evaluations] = gls_adapt(f, a, b, tol)
    % 实现自适应高斯-勒让德积分，返回积分结果和计算次数。
    [integral_result, output] = integral_qng(f, a, b, tol);
    evaluations = output.funcCount;
    
    result = integral_result;
end

% 计算最佳平方逼近多项式系数
function [a0, a1, a2, a3] = compute_coefficients(f_target, n)
    d = 1; % 对称区间的半径
    x = linspace(-d, d, 1000); % 定义x向量
    
    P = zeros(length(x), n+1); % 初始化Legendre多项式矩阵
    
    for k = 0:n
        P(:,k+1) = legendre_poly(k, x);
    end

    % 计算分子部分（积分值）
    [numerator, ~] = gls_adapt(f_target, -d, d, 1e-8); % 这里假设gls_adapt返回两个参数
    
    % 计算分母部分（Legendre多项式的归一化积分）
    denominator = zeros(1, n+1);
    
    for k = 0:n
        integrand = P(:,k+1) .^ 2; % 平方后的Legendre多项式
        [denom(k+1), ~] = gls_adapt(@(x) integrand, -d, d, 1e-8);
    end
    
    % 计算系数，避免除以零的情况
    a = zeros(1, n+1);
    for k = 0:n
        if denominator(k+1) ~= 0
            a(k+1) = numerator / denominator(k+1);
        else
            a(k+1) = 0; % 或者处理奇异情况，如增加正则化项
        end
    end
    
    [a0, a1, a2, a3] = deal(a(1), a(2), a(3), a(4)); % 根据n取值调整输出变量
end

% 主函数调用
function main(n)
    d = 1;
    f_target = @(x) sin(pi * x); % 目标函数
    
    [I, evaluations] = gls_adapt(f_target, -d, d, 1e-8);
    
    fprintf('积分结果: %.15f\n', I);
    fprintf('计算点数: %d\n', evaluations);
    
    exact_value = (2 * (pi^2 - 15)) / pi^3; % 精确值
    fprintf('精确值: %.15f\n', exact_value);
    fprintf('误差: %.15f\n', abs(I - exact_value));
    
    % 计算最佳平方逼近多项式的系数
    [a0, a1, a2, a3] = compute_coefficients(f_target, n);
end


function adaptive_gauss_legendre()
    exact = (2 * (pi^2 - 15)) / pi^3;
    f = @(x) 0.5 .* (5 .* x.^3 - 3 .* x) .* sin(pi * x);
    
    tol = 1e-8; % 容差
    
    [I, cnt] = gls_adapt(f, -1, 1, tol);
    
    fprintf('积分结果: %.15f\n', I);
    fprintf('计算点数: %d\n', cnt);
    fprintf('精确值: %.15f\n', exact);
    fprintf('误差: %.15f\n', abs(I - exact));
end

function [S, cnt] = gls_integrate(f, a, b)
    h = (b - a) / 2;
    c = h * sqrt(1/3); % 计算c的值 因为高斯勒让德
    x1 = a + c;
    x2 = b - c;
    f1 = f(x1);
    f2 = f(x2);
    S = h*(f1 + f2);
    cnt = 2;  % 因为用了两点积分
end

function [I, cnt] = gls_adapt(f, a, b, tol)
    % 初始计算
    S_current = gls_integrate(f, a, b);
    cnt_current = 2;  % 因为用了两点积分
    
    exact_val = pi;
    error = abs(S_current - exact_val);
    
    if error < tol
        I = S_current;
        return;
    else
        m = (a + b) / 2;
        
        [S_am, cnt_am] = gls_integrate(f, a, m);
        [S_mb, cnt_mb] = gls_integrate(f, m, b);
        S_new = S_am + S_mb;
        cnt = cnt_am + cnt_mb;
        
        error = abs(S_new - S_current);
        
        if error < tol / 2
            I = S_new;
            return;
        else
            % 继续递归处理每个子区间，使用更小的容限
            [I_left, cnt_left] = gls_adapt(f, a, m, tol/2);
            [I_right, cnt_right] = gls_adapt(f, m, b, tol/2);
            I = I_left + I_right;
            cnt = cnt_left + cnt_right;
        end
    end
end

% 调用主程序
main(3); % 替换为所需的多项式次数，例如3