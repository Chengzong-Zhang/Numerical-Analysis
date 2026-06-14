%% Adams-Bashforth 显式外插公式：统一覆盖 1 至 6 步
clear; clc; close all;                                      % 清理环境。
f = @(t,y) -y+exp(-t);                                     % 【必须替换】ODE 右端。
t0 = 0; T = 1; y0 = 0; h = 0.1; steps = 4;                % 【必须替换】区间、初值、步长、AB 步数 1...6。
AB = {[1],[3 -1]/2,[23 -16 5]/12,[55 -59 37 -9]/24,...    % cell 中依次存 AB1 到 AB4 的函数值权重。
      [1901 -2774 2616 -1274 251]/720,...                  % AB5 权重；续行符 ... 表示下一行仍属同一语句。
      [4277 -7923 9982 -7298 2877 -475]/1440};             % AB6 权重。
assert(ismember(steps,1:6),'steps 必须是 1 到 6。');         % ismember 检查选择是否合法。
coef = AB{steps};                                          % 用花括号取出所选阶数的数值向量。
N = round((T-t0)/h); assert(N>=steps-1,'区间太短，无法提供多步起步值。'); % 多步法需要足够节点。
t = t0+(0:N)*h; y = zeros(1,N+1); y(1) = y0;              % 构造网格并预分配。
for n = 1:steps-1                                          % 用 RK4 提供 y_1,...,y_{steps-1}。
    k1=f(t(n),y(n)); k2=f(t(n)+h/2,y(n)+h*k1/2);           % RK4 前两个阶段。
    k3=f(t(n)+h/2,y(n)+h*k2/2); k4=f(t(n)+h,y(n)+h*k3);   % RK4 后两个阶段。
    y(n+1)=y(n)+h*(k1+2*k2+2*k3+k4)/6;                    % 高精度起步，避免起步误差破坏多步法阶数。
end                                                        % 起步结束。
for n = steps:N                                            % 从已有 steps 个函数值的位置开始 AB 推进。
    fHist = zeros(steps,1);                                % 保存 [f_n,f_{n-1},...]。
    for j = 1:steps                                        % 收集历史函数值。
        fHist(j) = f(t(n-j+1),y(n-j+1));                   % MATLAB 下标换算：j=1 对应当前 f_n。
    end                                                    % 历史函数值收集结束。
    y(n+1) = y(n)+h*coef*fHist;                            % Adams-Bashforth 外插公式，完全显式。
end                                                        % 多步推进结束。
disp(table(t.',y.','VariableNames',{'t',sprintf('AB%d_y',steps)})); % 显示结果。
plot(t,y,'o-'); grid on; xlabel('t'); ylabel('y'); title(sprintf('AB%d',steps)); % 画图。
