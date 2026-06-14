%% Adams-Moulton 隐式内插与 Adams 预测-校正 PECE
clear; clc; close all;                                      % 清理环境。
f = @(t,y) -y+exp(-t);                                     % 【必须替换】ODE 右端。
t0 = 0; T = 1; y0 = 0; h = 0.1; order = 4;                % 【必须替换】区间、初值、步长、AM 阶数 1...6。
mode = "PECE"; tol = 1e-12; maxIter = 100;                 % 【必须替换】PECE 或 implicit；隐式迭代参数。
AB = {[1],[3 -1]/2,[23 -16 5]/12,[55 -59 37 -9]/24,[1901 -2774 2616 -1274 251]/720,[4277 -7923 9982 -7298 2877 -475]/1440}; % 预测权重。
AM = {[1],[1 1]/2,[5 8 -1]/12,[9 19 -5 1]/24,[251 646 -264 106 -19]/720,[475 1427 -798 482 -173 27]/1440}; % [f_{n+1},f_n,...]。
assert(ismember(order,1:6),'order 必须是 1 到 6。');         % 检查阶数。
ab = AB{order}; am = AM{order};                            % 取相同阶数的预测与校正系数。
N = round((T-t0)/h); t = t0+(0:N)*h;                      % 构造固定网格。
y = zeros(1,N+1); y(1) = y0;                              % 预分配数值解。
for n = 1:order-1                                          % RK4 提供多步法所需起步值。
    k1=f(t(n),y(n)); k2=f(t(n)+h/2,y(n)+h*k1/2);           % RK4 阶段 1、2。
    k3=f(t(n)+h/2,y(n)+h*k2/2); k4=f(t(n)+h,y(n)+h*k3);   % RK4 阶段 3、4。
    y(n+1)=y(n)+h*(k1+2*k2+2*k3+k4)/6;                    % 四阶起步值。
end                                                        % 起步结束。
for n = order:N                                            % 主多步循环。
    fHist = zeros(order,1);                                % [f_n,f_{n-1},...]。
    for j = 1:order                                        % 收集历史函数值。
        fHist(j)=f(t(n-j+1),y(n-j+1));                     % 对应 Adams 系数排列。
    end                                                    % 收集结束。
    zPred = y(n)+h*ab*fHist;                               % P：Adams-Bashforth 显式预报。
    if mode=="PECE"                                        % PECE 只校正一次，计算量小。
        z = y(n)+h*(am(1)*f(t(n+1),zPred)+am(2:end)*fHist(1:order-1)); % E+C：用预报函数值代入 AM。
    else                                                   % 完全隐式 AM：反复求解封闭公式。
        z = zPred;                                         % 用预测值作为隐式迭代初猜。
        for k = 1:maxIter                                  % 不动点迭代。
            zNew = y(n)+h*(am(1)*f(t(n+1),z)+am(2:end)*fHist(1:order-1)); % 更新隐式未知量。
            if abs(zNew-z)<tol, break; end                 % 达到容差则退出。
            z = zNew;                                      % 更新猜测。
        end                                                % 迭代结束。
        z = zNew;                                          % 保存最终迭代值。
    end                                                    % 模式选择结束。
    y(n+1)=z;                                              % E：保存校正值，供后续步骤使用。
end                                                        % 主循环结束。
disp(table(t.',y.','VariableNames',{'t',sprintf('AM%d_y',order)})); % 输出表格。
plot(t,y,'o-'); grid on; xlabel('t'); ylabel('y'); title(sprintf('AM%d %s',order,mode)); % 画图。
