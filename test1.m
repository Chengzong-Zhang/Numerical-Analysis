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
    c = h * sqrt(1/3); % 计算c的值
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