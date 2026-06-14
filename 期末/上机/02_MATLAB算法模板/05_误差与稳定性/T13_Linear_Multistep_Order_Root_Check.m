%% 一般线性多步法：阶数条件、特征根条件与绝对稳定检查
clear; clc;                                                % 清理环境。
alpha = [0,-1,1];                                          % 【必须替换】按 [alpha_0,...,alpha_k] 填写；默认是二步 Adams-Bashforth 左端。
beta = [-1/2,3/2,0];                                       % 【必须替换】按 [beta_0,...,beta_k] 填写；显式法最后一项通常为 0。
maxOrderToCheck = 8;                                       % 【必须替换】最多检查到多少阶。
assert(numel(alpha)==numel(beta),'alpha 与 beta 长度必须相同。'); % 同一个线性多步公式必须有同样多的系数。
k = numel(alpha)-1;                                        % k 是多步法步数；系数下标从数学上的 0 到 k。
j = 0:k;                                                   % j 保存公式中的节点编号。
c = zeros(1,maxOrderToCheck+1);                            % c(q+1) 对应教材中的阶条件系数 c_q。
c(1) = sum(alpha);                                         % c_0=sum(alpha_j)；相容方法必须满足 c_0=0。
for q = 1:maxOrderToCheck                                  % 依次检查 c_1,c_2,...。
    c(q+1)=sum(alpha.*j.^q)/factorial(q)-sum(beta.*j.^(q-1))/factorial(q-1); % Taylor 展开后的通用阶条件。
end                                                        % 阶条件计算结束。
firstNonzero = find(abs(c)>1e-12,1,'first');               % find(...,1,'first') 找第一个不满足零条件的位置。
if isempty(firstNonzero)                                   % 若检查范围内所有 c_q 都近似为零。
    fprintf('至少为 %d 阶；请增大 maxOrderToCheck 继续检查。\n',maxOrderToCheck); % 当前只能给出阶数下界。
else                                                       % 找到了首个非零 c_q。
    order = firstNonzero-2;                                % c_0 到 c_p 为零、c_{p+1} 非零时方法为 p 阶。
    fprintf('判定阶数 p = %d，首个非零系数 c_%d = %.6g。\n',order,order+1,c(firstNonzero)); % 输出阶数。
end                                                        % 阶数判断结束。
rhoRoots = roots(fliplr(alpha));                           % roots 要求按最高次到常数项排列，所以用 fliplr 反转。
rootCondition = all(abs(rhoRoots)<=1+1e-10) && all(arrayfun(@(r) sum(abs(rhoRoots-r)<1e-8)==1,rhoRoots(abs(abs(rhoRoots)-1)<1e-8))); % 单位圆内且单位圆根单重。
disp(table(rhoRoots,abs(rhoRoots),'VariableNames',{'rho_root','modulus'})); % 显示特征多项式 rho 的根。
fprintf('特征根条件是否满足：%d（1=满足，0=不满足）\n',rootCondition); % 相容加根条件可推出收敛。
z = -1;                                                    % 【必须替换】测试 z=lambda*h；衰减问题通常 z<0。
stabilityRoots = roots(fliplr(alpha-z*beta));              % 绝对稳定特征方程 rho(xi)-z*sigma(xi)=0。
fprintf('z=%.6g 时最大根模 = %.6g；小于 1 表示绝对稳定。\n',z,max(abs(stabilityRoots))); % 输出固定 z 的稳定性判断依据。
