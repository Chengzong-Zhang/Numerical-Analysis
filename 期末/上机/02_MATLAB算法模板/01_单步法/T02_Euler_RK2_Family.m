%% Euler 与教材二阶单步法统一模板
clear; clc; close all;                                      % 清理运行环境，确保结果只由本脚本产生。
f = @(t,y) -y+t.^2+1;                                      % 【必须替换】ODE 右端 f(t,y)；@ 创建匿名函数，.^ 支持逐元素幂。
t0 = 0; T = 1; y0 = 1; h = 0.1;                           % 【必须替换】区间、初值和固定步长。
method = "improvedEuler";                                  % 【必须替换】可选 forwardEuler/backwardEuler/trapezoid/improvedEuler/midpoint/heun2。
tol = 1e-12; maxIter = 100;                                % 隐式法不动点迭代的停止阈值与最大迭代次数。
N = round((T-t0)/h);                                       % 固定步长方法的总步数。
assert(abs(t0+N*h-T)<1e-12,'步长必须整除区间长度。');       % 保证最后一个节点恰好等于 T。
t = t0+(0:N)*h;                                            % 构造等距网格 t_0,...,t_N。
y = zeros(1,N+1); y(1) = y0;                              % 预分配数值解并写入初值。
for n = 1:N                                                % 从数学下标 n=0 推进到 n=N-1。
    tn = t(n); yn = y(n);                                  % 保存当前节点，便于公式与教材一一对应。
    switch method                                          % switch 根据字符串选择一种算法分支。
        case "forwardEuler"                                % 向前 Euler：显式、一阶、稳定区间 (-2,0)。
            y(n+1) = yn+h*f(tn,yn);                        % y_{n+1}=y_n+h f_n。
        case "backwardEuler"                               % 向后 Euler：隐式、一阶、A 稳定，适合刚性问题。
            z = yn+h*f(tn,yn);                             % 用向前 Euler 值作为隐式方程的初始猜测。
            for k = 1:maxIter                              % 不动点迭代求 z=yn+h*f(t_{n+1},z)。
                zNew = yn+h*f(t(n+1),z);                   % 一次不动点更新；非线性 f 时不能直接解出。
                if abs(zNew-z)<tol, break; end             % 单行 if：相邻迭代足够接近就停止。
                z = zNew;                                  % 更新猜测值，继续迭代。
            end                                            % 结束隐式迭代。
            y(n+1) = zNew;                                 % 保存向后 Euler 的隐式解。
        case "trapezoid"                                   % 隐式梯形：二阶、A 稳定。
            z = yn+h*f(tn,yn);                             % 用 Euler 预测值启动迭代。
            for k = 1:maxIter                              % 解 z=yn+h/2*(f_n+f_{n+1})。
                zNew = yn+h/2*(f(tn,yn)+f(t(n+1),z));      % 梯形公式来自对积分形式使用梯形求积。
                if abs(zNew-z)<tol, break; end             % 达到容差后停止。
                z = zNew;                                  % 更新当前迭代值。
            end                                            % 结束隐式迭代。
            y(n+1) = zNew;                                 % 保存隐式梯形结果。
        case "improvedEuler"                               % 改进 Euler，也叫显式梯形或 Euler 预测校正。
            k1 = f(tn,yn);                                 % 左端斜率，预测阶段 P。
            k2 = f(tn+h,yn+h*k1);                          % 在 Euler 预测点计算右端斜率，计算阶段 E。
            y(n+1) = yn+h*(k1+k2)/2;                       % 两端斜率平均，校正阶段 C；整体二阶。
        case "midpoint"                                    % 变形 Euler/中点 RK2，整体二阶。
            k1 = f(tn,yn);                                 % 当前点斜率。
            k2 = f(tn+h/2,yn+h*k1/2);                      % 用 k1 预测中点，再取中点斜率。
            y(n+1) = yn+h*k2;                              % 用中点斜率推进整步。
        case "heun2"                                       % 教材二阶 Heun：a2=2/3，权重 1/4 与 3/4。
            k1 = f(tn,yn);                                 % 第一阶段斜率。
            k2 = f(tn+2*h/3,yn+2*h*k1/3);                  % 在 2/3 步位置取第二阶段斜率。
            y(n+1) = yn+h*(k1+3*k2)/4;                     % 加权组合满足二阶 RK 阶条件。
        otherwise                                          % 捕获拼写错误或未支持的方法名。
            error('未知 method：%s',method);                % error 主动终止并显示错误信息。
    end                                                    % 结束算法选择。
end                                                        % 结束时间推进。
disp(table(t.',y.','VariableNames',{'t',char(method)}));    % char 把 string 转成 table 可用的变量名。
plot(t,y,'o-','LineWidth',1.3); grid on; xlabel('t'); ylabel('y'); title(method); % 画数值解。
