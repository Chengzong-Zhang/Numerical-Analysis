%% 误差、收敛阶、Richardson 外推与绝对稳定性模板
clear; clc; close all;                                      % 清理环境。
exact = @(t) exp(-t);                                      % 【必须替换】精确解；若未知可用更小步长高阶解作参考。
approx = @(h) (1-h)^(1/h);                                 % 【必须替换】返回终点数值解的函数；默认是 Euler 解 y'=-y 到 t=1。
hs = [0.2,0.1,0.05,0.025];                                % 【必须替换】逐次减半最方便估计阶数。
E = zeros(size(hs));                                       % size(hs) 创建同形状误差数组。
for i = 1:numel(hs)                                        % 遍历全部步长。
    E(i)=abs(approx(hs(i))-exact(1));                       % 终点绝对误差；若题目要求最大误差，应对所有网格点取 max。
end                                                        % 误差计算结束。
p = nan(size(hs));                                         % nan 表示第一项没有前一组误差可比较。
for i = 2:numel(hs)                                        % 从第二个步长开始估计收敛阶。
    p(i)=log(E(i-1)/E(i))/log(hs(i-1)/hs(i));              % 一般步长比公式；减半时分母就是 log(2)。
end                                                        % 阶数估计结束。
assumedOrder = 1;                                          % 【必须替换】Richardson 外推所依据的方法理论阶数。
yCoarse = approx(hs(end-1)); yFine = approx(hs(end));      % 取最后两组粗细网格结果。
yRich = yFine+(yFine-yCoarse)/((hs(end-1)/hs(end))^assumedOrder-1); % 消去主误差项，提高一阶或更多。
disp(table(hs.',E.',p.','VariableNames',{'h','error','estimated_order'})); % 输出误差与阶数表。
fprintf('Richardson 外推终值 = %.15g\n',yRich);              % 打印外推结果。
z = linspace(-5,1,1200);                                   % 【必须替换】z=lambda*h；沿实轴检查稳定区间。
R_Euler = 1+z;                                             % 向前 Euler 稳定函数 R(z)=1+z。
R_ImprovedEuler = 1+z+z.^2/2;                              % 改进 Euler/RK2 稳定函数。
R_RK4 = 1+z+z.^2/2+z.^3/6+z.^4/24;                        % 经典 RK4 稳定函数。
plot(z,abs(R_Euler),z,abs(R_ImprovedEuler),z,abs(R_RK4),'LineWidth',1.2); % 画 |R(z)|。
yline(1,'k--'); grid on; ylim([0,3]);                       % |R(z)|<1 的区域对应绝对稳定。
xlabel('z=\lambda h'); ylabel('|R(z)|'); legend('Euler','RK2','RK4','|R|=1'); title('实轴绝对稳定性'); % 图示说明。
