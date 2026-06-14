%% 方程组与高阶 ODE：先化为一阶系统，再统一使用 RK4
clear; clc; close all;                                      % 清理环境。
F = @(t,w) [w(2); exp(2*t)*sin(t)-2*w(1)+2*w(2)];          % 【必须替换】二阶 y''=g(t,y,y') 写成 [y';g]；一般方程组也直接写列向量。
t0 = 0; T = 1; w0 = [-0.4;-0.6]; h = 0.1;                 % 【必须替换】w0=[y(t0);y'(t0);...]，区间与步长。
N = round((T-t0)/h); assert(abs(t0+N*h-T)<1e-12,'步长必须整除区间长度。'); % 检查固定网格。
t = t0+(0:N)*h;                                            % 生成时间节点。
m = numel(w0);                                             % 系统一阶方程个数，也等于状态向量维数。
w = zeros(m,N+1); w(:,1)=w0(:);                           % 每列保存一个时刻的全部分量；(:) 强制列向量。
for n = 1:N                                                % 对向量系统逐步推进。
    K1 = F(t(n),w(:,n));                                   % 第一阶段斜率向量。
    K2 = F(t(n)+h/2,w(:,n)+h*K1/2);                        % 第二阶段在预测中点计算。
    K3 = F(t(n)+h/2,w(:,n)+h*K2/2);                        % 第三阶段再次修正中点。
    K4 = F(t(n)+h,w(:,n)+h*K3);                            % 第四阶段在右端点计算。
    w(:,n+1)=w(:,n)+h*(K1+2*K2+2*K3+K4)/6;                % 向量加权公式与标量 RK4 完全相同。
end                                                        % 推进结束。
disp(array2table([t.',w.'],'VariableNames',["t","w"+string(1:m)])); % array2table 把数值矩阵转成带列名表格。
plot(t,w.','LineWidth',1.3); grid on; xlabel('t'); ylabel('components'); legend("w"+string(1:m)); % 一次画全部分量。
