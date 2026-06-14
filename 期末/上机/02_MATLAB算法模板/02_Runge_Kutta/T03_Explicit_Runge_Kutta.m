%% 一般显式 Runge-Kutta：用 Butcher 表覆盖三阶 Heun、三阶 RK、经典 RK4
clear; clc; close all;                                      % 清理环境。
f = @(t,y) -y+t.^2+1;                                      % 【必须替换】标量或列向量 ODE 右端。
t0 = 0; T = 1; y0 = 1; h = 0.1;                           % 【必须替换】区间、初值、步长；向量初值写成 [a;b]。
method = "RK4";                                            % 【必须替换】可选 Heun3/RK3/RK4。
switch method                                              % 根据方法装入 Butcher 系数 A,b,c。
    case "Heun3"                                           % 教材三级三阶 Heun 方法。
        A = [0 0 0; 1/3 0 0; 0 2/3 0];                    % A 严格下三角，所以阶段可依次显式计算。
        b = [1/4;0;3/4]; c = [0;1/3;2/3];                 % b 是最终权重，c 是阶段时间位置。
    case "RK3"                                             % 教材公式 (5.3.13) 的三级三阶 RK。
        A = [0 0 0; 1/2 0 0; 0 3/4 0];                    % 每行给出该阶段对前面阶段的组合。
        b = [2/9;3/9;4/9]; c = [0;1/2;3/4];               % 权重和为 1，保证至少一阶相容。
    case "RK4"                                             % 最常用的经典四级四阶 RK。
        A = [0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0];    % 经典 RK4 的阶段矩阵。
        b = [1/6;1/3;1/3;1/6]; c = [0;1/2;1/2;1];         % 经典 RK4 权重与节点。
    otherwise                                              % 未知名称时停止。
        error('未知 method：%s',method);                    % 显示错误方法名。
end                                                        % Butcher 表选择结束。
N = round((T-t0)/h);                                       % 计算总步数。
assert(abs(t0+N*h-T)<1e-12,'步长必须整除区间长度。');       % 检查终点。
t = t0+(0:N)*h;                                            % 生成等距时间网格。
y0 = y0(:);                                                % 冒号把标量/行向量统一成列向量。
y = zeros(numel(y0),N+1); y(:,1) = y0;                    % 每一列保存一个时间点的状态向量。
s = numel(b);                                              % RK 级数等于权重个数。
for n = 1:N                                                % 逐步推进。
    K = zeros(numel(y0),s);                                % K(:,r) 保存第 r 个阶段斜率。
    for r = 1:s                                            % 依次计算所有显式阶段。
        stageY = y(:,n)+h*K*A(r,:).';                      % 矩阵乘法组合已算阶段；.' 是非共轭转置。
        K(:,r) = f(t(n)+c(r)*h,stageY);                    % 在阶段点计算斜率 f。
    end                                                    % 阶段计算结束。
    y(:,n+1) = y(:,n)+h*K*b;                               % 用权重 b 加权所有阶段，完成一步。
end                                                        % 时间推进结束。
disp(table(t.',y(1,:).','VariableNames',{'t','第一分量'})); % 展示第一分量；方程组可自行增加其他列。
plot(t,y(1,:),'o-','LineWidth',1.3); grid on; xlabel('t'); ylabel('y_1'); title(method); % 画第一分量。
