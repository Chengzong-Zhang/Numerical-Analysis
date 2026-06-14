%% 隐式 ODE 方法所需：不动点迭代与 Newton-Raphson 迭代
clear; clc;                                                % 清理环境。
G = @(z) z-cos(z);                                         % 【必须替换】把隐式一步公式整理成 G(z)=0。
dG = @(z) 1+sin(z);                                        % 【必须替换】G 对 z 的导数；向量问题应换成 Jacobian 矩阵。
g = @(z) cos(z);                                           % 【必须替换】若能整理成 z=g(z)，供不动点迭代使用。
z0 = 0.5; tol = 1e-12; maxIter = 100;                     % 【必须替换】初始猜测、容差、迭代上限。
zFixed = z0;                                               % 不动点迭代当前值。
for k = 1:maxIter                                          % 不动点迭代要求局部 |g'(z*)|<1 才易收敛。
    zNew = g(zFixed);                                      % 计算 z^{k+1}=g(z^k)。
    if abs(zNew-zFixed)<tol, break; end                    % 相邻迭代差作为实用停止准则。
    zFixed = zNew;                                         % 更新迭代值。
end                                                        % 不动点迭代结束。
zNewton = z0;                                              % Newton 迭代当前值。
for k = 1:maxIter                                          % Newton 在根附近通常二次收敛。
    delta = G(zNewton)/dG(zNewton);                         % 标量 Newton 修正量；向量问题写成 dG(z)\G(z)，不要用 inv。
    zNewtonNew = zNewton-delta;                            % Newton 更新 z^{k+1}=z^k-G/G'。
    if abs(zNewtonNew-zNewton)<tol, break; end              % 修正量足够小时停止。
    zNewton = zNewtonNew;                                  % 更新迭代值。
end                                                        % Newton 迭代结束。
fprintf('fixed-point = %.15g, residual = %.3e\n',zNew,abs(G(zNew))); % residual 检查是否真正满足方程。
fprintf('Newton      = %.15g, residual = %.3e\n',zNewtonNew,abs(G(zNewtonNew))); % 比较 Newton 结果。
