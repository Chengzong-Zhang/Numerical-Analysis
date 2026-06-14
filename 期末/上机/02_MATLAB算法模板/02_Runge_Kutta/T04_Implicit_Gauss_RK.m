%% 隐式 Gauss-Runge-Kutta：一级二阶、二级四阶、三级六阶
clear; clc; close all;                                      % 清理环境。
f = @(t,y) -10*(y-cos(t))-sin(t);                          % 【必须替换】ODE 右端；默认例题含快速衰减项，体现刚性。
t0 = 0; T = 1; y0 = 1; h = 0.1;                           % 【必须替换】区间、初值、步长。
method = "Gauss4"; tol = 1e-12; maxIter = 100;             % 【必须替换】Gauss2/Gauss4/Gauss6；迭代容差与上限。
switch method                                              % 装入隐式 RK 的 Butcher 表。
    case "Gauss2"                                          % 一级二阶隐式中点法。
        A = 1/2; b = 1; c = 1/2;                           % 单阶段也要解隐式阶段方程。
    case "Gauss4"                                          % 二级四阶 Gauss-RK。
        q = sqrt(3);                                       % sqrt(3) 是 Gauss 两节点系数的一部分。
        A = [1/4,1/4-q/6;1/4+q/6,1/4];                    % A 含对角及上三角项，因此各阶段耦合。
        b = [1/2;1/2]; c = [1/2-q/6;1/2+q/6];             % Gauss 节点与求积权重。
    case "Gauss6"                                          % 三级六阶 Gauss-RK。
        q = sqrt(15);                                      % 三节点 Gauss 系数。
        A = [5/36,2/9-q/15,5/36-q/30;5/36+q/24,2/9,5/36-q/24;5/36+q/30,2/9+q/15,5/36]; % 阶段矩阵。
        b = [5/18;4/9;5/18]; c = [1/2-q/10;1/2;1/2+q/10]; % 六阶方法的权重和节点。
    otherwise                                              % 方法名错误时终止。
        error('未知 method：%s',method);                    % 报告错误。
end                                                        % 系数选择结束。
N = round((T-t0)/h); assert(abs(t0+N*h-T)<1e-12,'步长必须整除区间长度。'); % 固定网格检查。
t = t0+(0:N)*h; y = zeros(1,N+1); y(1) = y0;              % 生成网格并预分配标量解。
s = numel(b);                                              % 隐式 RK 的阶段数。
for n = 1:N                                                % 逐步推进。
    K = repmat(f(t(n),y(n)),s,1);                          % repmat 复制初始斜率，作为全部阶段的初猜。
    for iter = 1:maxIter                                   % 对耦合阶段方程做不动点迭代。
        Knew = zeros(s,1);                                 % 预分配新阶段向量。
        for r = 1:s                                        % 更新每个阶段。
            stageY = y(n)+h*A(r,:)*K;                      % 当前全部阶段共同决定第 r 个阶段状态。
            Knew(r) = f(t(n)+c(r)*h,stageY);               % 计算新的第 r 阶段斜率。
        end                                                % 本轮阶段更新结束。
        if norm(Knew-K,inf)<tol, break; end                 % 无穷范数控制所有阶段的最大变化。
        K = Knew;                                          % 接受新阶段值并继续迭代。
    end                                                    % 隐式阶段迭代结束。
    assert(iter<maxIter,'隐式阶段未收敛：请减小 h 或改用 Newton 迭代。'); % 防止静默返回错误解。
    y(n+1) = y(n)+h*b.'*Knew;                              % 使用隐式阶段的加权和推进一步。
end                                                        % 时间推进结束。
disp(table(t.',y.','VariableNames',{'t',char(method)}));    % 显示结果。
plot(t,y,'o-'); grid on; xlabel('t'); ylabel('y'); title(method); % 画图。
