%% Runge-Kutta-Fehlberg 4(5) 自适应步长模板
clear; clc; close all;                                      % 清理环境。
f = @(t,y) -y+t.^2+1;                                      % 【必须替换】ODE 右端。
t0 = 0; T = 1; y0 = 1;                                    % 【必须替换】区间与初值。
h = 0.1; hMin = 1e-6; hMax = 0.2; tol = 1e-6;             % 【必须替换】初始/最小/最大步长与每步局部误差限。
t = t0; y = y0; hs = [];                                   % 自适应网格长度未知，先用动态数组保存。
while t(end)<T                                             % while 在尚未到终点时继续。
    hn = min([h,hMax,T-t(end)]);                           % min 保证不超过最大步长，并让最后一步恰到 T。
    tn = t(end); yn = y(end);                              % end 关键字表示数组最后一个元素。
    K1 = f(tn,yn);                                         % Fehlberg 第 1 阶段斜率。
    K2 = f(tn+hn/4,yn+hn*K1/4);                            % 第 2 阶段。
    K3 = f(tn+3*hn/8,yn+hn*(3*K1+9*K2)/32);                % 第 3 阶段。
    K4 = f(tn+12*hn/13,yn+hn*(1932*K1-7200*K2+7296*K3)/2197); % 第 4 阶段。
    K5 = f(tn+hn,yn+hn*(439*K1/216-8*K2+3680*K3/513-845*K4/4104)); % 第 5 阶段。
    K6 = f(tn+hn/2,yn+hn*(-8*K1/27+2*K2-3544*K3/2565+1859*K4/4104-11*K5/40)); % 第 6 阶段。
    y4 = yn+hn*(25*K1/216+1408*K3/2565+2197*K4/4104-K5/5); % 嵌套四阶近似。
    y5 = yn+hn*(16*K1/135+6656*K3/12825+28561*K4/56430-9*K5/50+2*K6/55); % 嵌套五阶近似。
    err = abs(y5-y4);                                      % 两个嵌套近似之差作为可计算的局部误差估计。
    if err<=tol                                            % 只有误差合格时才接受本步。
        t(end+1) = tn+hn;                                  % end+1 在数组末尾追加新时间点。
        y(end+1) = y5;                                     % 接受精度更高的五阶值。
        hs(end+1) = hn;                                    % 记录被接受的实际步长。
    end                                                    % 拒绝步时不更新 t 和 y，只调整 h 后重算。
    factor = 4*(err==0)+0.9*(tol/max(err,eps))^(1/5)*(err~=0); % 安全系数 0.9；eps 防止除零；五阶误差用 1/5 次幂。
    factor = min(4,max(0.1,factor));                       % 限制单次步长变化倍数，避免剧烈震荡。
    h = min(hMax,max(hMin,hn*factor));                     % 把新步长限制在 [hMin,hMax]。
    assert(err<=tol || h>hMin,'达到 hMin 仍不能满足误差限。'); % 无法继续时明确报错。
end                                                        % 自适应循环结束。
disp(table(t.',y.','VariableNames',{'t','RKF45_y'}));       % 自适应节点不等距，直接打印实际节点。
stairs(t(1:end-1),hs,'LineWidth',1.3); grid on; xlabel('t'); ylabel('accepted h'); title('RKF45 实际步长'); % 画步长变化。
