%% Milne、Hamming、Nyström 多步法模板
clear; clc; close all;                                      % 清理环境。
f = @(t,y) -y+exp(-t);                                     % 【必须替换】ODE 右端。
t0 = 0; T = 1; y0 = 0; h = 0.1;                           % 【必须替换】区间、初值、步长。
method = "Milne";                                          % 【必须替换】Milne/Hamming/Nystrom。
N = round((T-t0)/h); t = t0+(0:N)*h;                      % 固定步长网格。
y = zeros(1,N+1); y(1)=y0;                                % 预分配数值解。
startCount = 3*(method~="Nystrom")+1*(method=="Nystrom");  % Milne/Hamming 需 4 个点，Nyström 需 2 个点。
for n = 1:startCount                                       % 用 RK4 计算所需起步值。
    k1=f(t(n),y(n)); k2=f(t(n)+h/2,y(n)+h*k1/2);           % RK4 前两阶段。
    k3=f(t(n)+h/2,y(n)+h*k2/2); k4=f(t(n)+h,y(n)+h*k3);   % RK4 后两阶段。
    y(n+1)=y(n)+h*(k1+2*k2+2*k3+k4)/6;                    % 保存起步值。
end                                                        % 起步结束。
if method=="Nystrom"                                       % Nyström/Leapfrog：y_{n+1}=y_{n-1}+2h f_n。
    for n = 2:N                                            % 从已有 y_0,y_1 开始。
        y(n+1)=y(n-1)+2*h*f(t(n),y(n));                    % 二步显式中心格式；寄生根 -1 会导致振荡。
    end                                                    % Nyström 推进结束。
else                                                       % Milne 与 Hamming 共用 Milne 预测式。
    for n = 4:N                                            % 已有 y_0,...,y_3 后开始推进。
        fp=f(t(n),y(n)); fm1=f(t(n-1),y(n-1)); fm2=f(t(n-2),y(n-2)); % 历史函数值。
        yp=y(n-3)+4*h*(2*fp-fm1+2*fm2)/3;                  % Milne 四步显式预测公式。
        fPred=f(t(n+1),yp);                                % 在预测点计算 f_{n+1}。
        if method=="Milne"                                 % Milne 校正公式。
            y(n+1)=y(n-1)+h*(fPred+4*fp+fm1)/3;            % Simpson 型隐式校正，只代入预测值一次。
        else                                               % Hamming 四阶校正公式。
            y(n+1)=(9*y(n)-y(n-2)+3*h*(fPred+2*fp-fm1))/8; % Hamming 通过待定系数提高误差性质。
        end                                                % 校正方法选择结束。
    end                                                    % Milne/Hamming 推进结束。
end                                                        % 方法大分支结束。
disp(table(t.',y.','VariableNames',{'t',char(method)}));    % 输出结果。
plot(t,y,'o-'); grid on; xlabel('t'); ylabel('y'); title(method); % 画图。

