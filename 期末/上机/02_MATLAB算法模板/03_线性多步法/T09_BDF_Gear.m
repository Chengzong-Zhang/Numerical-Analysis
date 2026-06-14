%% BDF/Gear 1 至 6 阶隐式模板
clear; clc; close all;                                      % 清理环境。
f = @(t,y) -10*(y-cos(t))-sin(t);                          % 【必须替换】ODE 右端；默认例题精确解为 cos(t)，并带快速衰减项。
fy = @(t,y) -10;                                           % 【必须替换】偏导 df/dy；Newton 解隐式方程时必须提供。
t0 = 0; T = 1; y0 = 1; h = 0.05; k = 2;                   % 【必须替换】区间、初值、步长、BDF 阶数 1...6。
tol = 1e-12; maxIter = 200;                                % Newton 迭代的停止阈值与最大迭代次数。
alpha = {[1 -1],[1 -4/3 1/3],[1 -18/11 9/11 -2/11],...    % 归一化 BDF1 至 BDF3 左端系数。
         [1 -48/25 36/25 -16/25 3/25],...                 % BDF4 左端系数。
         [1 -300/137 300/137 -200/137 75/137 -12/137],... % BDF5 左端系数。
         [1 -360/147 450/147 -400/147 225/147 -72/147 10/147]}; % BDF6 左端系数。
beta = [1,2/3,6/11,12/25,60/137,60/147];                  % 右端 h*beta_k*f(t_{n+1},y_{n+1}) 系数。
assert(ismember(k,1:6),'k 必须是 1 到 6。');               % 检查阶数合法。
a = alpha{k}; b = beta(k);                                 % 取所选 BDF 系数。
N = round((T-t0)/h); t = t0+(0:N)*h;                      % 构造网格。
y = zeros(1,N+1); y(1)=y0;                                % 预分配数值解。
for n = 1:k-1                                              % 用 RK4 提供 k-1 个起步值。
    q1=f(t(n),y(n)); q2=f(t(n)+h/2,y(n)+h*q1/2);           % RK4 阶段 1、2。
    q3=f(t(n)+h/2,y(n)+h*q2/2); q4=f(t(n)+h,y(n)+h*q3);   % RK4 阶段 3、4。
    y(n+1)=y(n)+h*(q1+2*q2+2*q3+q4)/6;                    % 保存起步值。
end                                                        % 起步结束。
for n = k:N                                                % 从已有 k 个历史值开始推进。
    known = 0;                                             % known 表示移项后的已知历史组合。
    for j = 1:k                                            % a(2)...a(k+1) 分别乘 y_n,...,y_{n-k+1}。
        known = known-a(j+1)*y(n-j+1);                     % 因 a(1)=1，移项后历史项取负号。
    end                                                    % 历史组合完成。
    z = y(n);                                              % 用最近解作为 y_{n+1} 的 Newton 初猜。
    for iter = 1:maxIter                                   % 解 G(z)=z-known-h*b*f(t_{n+1},z)=0。
        G = z-known-h*b*f(t(n+1),z);                       % 隐式 BDF 方程的残差。
        dG = 1-h*b*fy(t(n+1),z);                           % 残差导数；向量题应换成 I-h*b*Jacobian。
        zNew = z-G/dG;                                     % 标量 Newton 更新；向量题写成 z-dG\G。
        if abs(zNew-z)<tol, break; end                     % Newton 修正足够小时停止。
        z = zNew;                                          % 更新猜测并继续 Newton 迭代。
    end                                                    % 隐式迭代结束。
    assert(iter<maxIter,'BDF 隐式迭代未收敛，请改用 Newton 或减小 h。'); % 显式提示失败。
    y(n+1)=zNew;                                           % 保存本步 BDF 解。
end                                                        % BDF 推进结束。
disp(table(t.',y.','VariableNames',{'t',sprintf('BDF%d_y',k)})); % 输出结果。
plot(t,y,'o-'); hold on; plot(t,cos(t),'k--'); grid on; legend('BDF','exact'); title(sprintf('BDF%d',k)); % 对比默认真解。
