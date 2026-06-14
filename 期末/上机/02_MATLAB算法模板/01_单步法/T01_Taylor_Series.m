%% Taylor 级数法模板：y_{n+1}=y_n+sum_{j=1}^p h^j/j! * y^(j)(t_n)
clear; clc; close all;                                      % MATLAB 语法：清变量、清命令窗、关图窗，避免旧数据干扰。
fDeriv = {@(t,y) -y+t.^2+1, @(t,y) y-t.^2+2*t-1};         % 【必须替换】依次填写 y'、y''、...；花括号创建 cell 数组。
t0 = 0; T = 1; y0 = 1; h = 0.1;                           % 【必须替换】左端点、右端点、初值、固定步长。
p = numel(fDeriv);                                         % numel 返回导数函数个数，也就是 Taylor 方法阶数 p。
N = round((T-t0)/h);                                       % round 将步数取整；本模板要求 (T-t0)/h 为整数。
assert(abs(t0+N*h-T)<1e-12,'步长必须整除区间长度。');       % assert 在条件不成立时停止，防止终点与题目不一致。
t = t0+(0:N)*h;                                            % 冒号语法 0:N 生成整数网格，向量化得到所有节点。
y = zeros(1,N+1);                                          % zeros 预分配内存；N 步包含 N+1 个节点。
y(1) = y0;                                                 % MATLAB 下标从 1 开始，所以 y(1) 对应数学上的 y_0。
for n = 1:N                                                % for 循环逐步推进；每步仅使用当前点信息，属于单步法。
    increment = 0;                                         % 累积本步 Taylor 增量，先置零。
    for j = 1:p                                            % 依次加入 h*y'、h^2/2!*y'' 等项。
        increment = increment + h^j/factorial(j)*fDeriv{j}(t(n),y(n)); % cell 用 {j} 取函数；factorial(j)=j!。
    end                                                    % 结束内层 Taylor 求和。
    y(n+1) = y(n)+increment;                               % 截断到 p 阶导数，局部截断误差通常为 O(h^(p+1))。
end                                                        % 结束时间推进。
disp(table(t.',y.','VariableNames',{'t','Taylor_y'}));      % 点转置 .' 把行向量变列向量；table 生成结果表。
plot(t,y,'o-','LineWidth',1.3); grid on;                    % plot 画数值解；grid on 打开网格便于观察。
xlabel('t'); ylabel('y'); title(sprintf('%d 阶 Taylor 方法',p)); % sprintf 把阶数写入标题。
