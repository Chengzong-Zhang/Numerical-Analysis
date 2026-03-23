%% =========================================================
%  machine_precision.m
%  功能：计算并验证计算机的机器精度、下溢值、上溢值
%  =========================================================
%
%  【数值分析背景知识】
%  计算机使用"浮点数"来表示实数，遵循 IEEE 754 双精度标准。
%  双精度浮点数的结构：
%    符号位(1位) + 指数位(11位) + 尾数位(52位) = 共64位
%
%  由此产生三个重要极限：
%  1. 机器精度 (eps)：1.0 与比1.0大的最小浮点数之间的差
%                     理论值 = 2^(-52) ≈ 2.2204e-16
%  2. 下溢值 (realmin)：最小正规格化浮点数
%                       理论值 = 2^(-1022) ≈ 2.2251e-308
%  3. 上溢值 (realmax)：最大有限浮点数
%                       理论值 = (2-2^(-52)) * 2^1023 ≈ 1.7977e+308
%
%  =========================================================

% 清除工作区变量、关闭所有图窗、清空命令行窗口
% clear：删除工作区中所有变量，避免上次运行的残留数据干扰
% close all：关闭所有打开的图形窗口
% clc：清空命令行窗口的显示内容（不删除变量）
clear; close all; clc;

fprintf('==============================================\n');
fprintf('   计算机浮点数特性：机器精度、下溢、上溢\n');
fprintf('==============================================\n\n');
% fprintf：格式化输出函数，类似C语言的printf
% \n：换行符
% \t：制表符（后面会用到）


%% =========================================================
%  第一部分：计算机器精度 (Machine Epsilon)
%  =========================================================
%
%  【数值分析】
%  机器精度 eps 定义为：满足 1.0 + eps > 1.0 的最小正数 eps
%  它衡量了浮点运算的相对误差上界。
%  算法思路：
%    从 eps = 1 开始，不断将其除以2，
%    直到 1.0 + eps 在计算机中与 1.0 无法区分（即相等）为止，
%    此时上一步的 eps 就是机器精度。
%  =========================================================

fprintf('--- 第一部分：机器精度 (Machine Epsilon) ---\n\n');

% 初始化 eps_calc 为 1.0
% 变量名后不加分号";"会在命令行显示结果；加";"则不显示
% 这里加";"是为了保持输出整洁
eps_calc = 1.0;

% while 循环：当括号内条件为真时，反复执行循环体
% 条件：(1.0 + eps_calc/2) > 1.0
%   即：将 eps_calc 再缩小一半后，加到1.0上，结果还是否大于1.0
%   如果大于，说明还没到精度极限，继续缩小
while (1.0 + eps_calc / 2) > 1.0
    eps_calc = eps_calc / 2;  % 将机器精度候选值折半
end
% 循环结束时，eps_calc 就是机器精度
% 此时 eps_calc/2 已经小到计算机无法区分 1.0 与 1.0+eps_calc/2

% 输出我们算法计算的结果
fprintf('【算法计算】机器精度 eps = %.16e\n', eps_calc);
% %.16e：科学计数法格式，保留16位有效小数

% 输出 MATLAB 内置常量进行对比
% eps 是 MATLAB 预定义的机器精度常量
fprintf('【MATLAB内置】eps       = %.16e\n', eps);
fprintf('【理论值】    2^(-52)   = %.16e\n\n', 2^(-52));

% 验证：计算相对误差
% abs()：取绝对值的函数
relative_error_eps = abs(eps_calc - eps) / eps;
fprintf('算法结果与内置eps的相对误差 = %.2e\n\n', relative_error_eps);
% %.2e：科学计数法，保留2位小数


%% =========================================================
%  第二部分：计算下溢值 (Underflow) ——正规化边界检测算法
%  =========================================================
%
%  【数值分析】
%  下溢值 = 最小正规格化浮点数 realmin = 2^(-1022)
%  正规化浮点数形如 1.f × 2^e，尾数前导位隐含为1；
%  当指数 e 降到最小值 e_min = -1022 后无法继续减小，此点即为下溢边界。
%
%  【正确算法思路：基于正规化检测】
%  关键观察：正规化数 x 的相对精度为 eps（即 x + x*eps > x）。
%  当 x 变成非规格化数时，x*eps 会小于最小可表示正数，舍入为零，
%  导致 x + x*eps == x，相对精度"消失"。
%  利用这一性质：若 x/2 + (x/2)*eps == x/2，说明 x/2 已是非规格化数，
%  则 x 就是最小正规格化数（下溢值）。
%  =========================================================

fprintf('--- 第二部分：下溢值 (Underflow) ---\n\n');

underflow_calc = 1.0;

% 正规化检测：只要 x/2 的相对精度 eps 仍可观测，说明 x/2 还是规格化数
% 条件：x/2 + (x/2)*eps > x/2  ← x/2 是正规化数
% 反之：x/2 + (x/2)*eps == x/2 ← x/2 已是非规格化数，停止
while underflow_calc / 2 + (underflow_calc / 2) * eps > underflow_calc / 2
    underflow_calc = underflow_calc / 2;
end
% 循环结束：underflow_calc/2 刚好变成非规格化数
% 所以 underflow_calc 是最小的正规格化浮点数，即下溢值

fprintf('【算法计算】下溢值（最小正规格化数）= %.16e\n', underflow_calc);
fprintf('【MATLAB内置】realmin             = %.16e\n', realmin);
fprintf('【理论值】    2^(-1022)           = %.16e\n\n', 2^(-1022));

relative_error_underflow = abs(underflow_calc - realmin) / realmin;
fprintf('算法结果与realmin的相对误差 = %.2e\n\n', relative_error_underflow);


%% =========================================================
%  第三部分：计算上溢值 (Overflow)
%  =========================================================
%
%  【数值分析】
%  上溢 (Overflow)：当浮点数超过最大可表示值时，结果变为 Inf（无穷大）。
%  最大有限浮点数 realmax = (2 - 2^(-52)) * 2^1023 ≈ 1.7977e+308
%
%  精确算法（二进制逼近法）分两阶段：
%  第一阶段：翻倍找到最大的2的幂次（即 2^1023）
%  第二阶段：从该幂次出发，用二分法逐步加上更小的量向 realmax 逼近
%    每次尝试加上当前步长 step，若不溢出则接受，否则将 step 减半继续
%    这本质上是在二进制位上从高位到低位逐一填充尾数的1
%  =========================================================

fprintf('--- 第三部分：上溢值 (Overflow) ---\n\n');

% ---- 第一阶段：翻倍，找到最大的2的整数次幂 ----
overflow_calc = 1.0;

% ~isinf()：~是逻辑非，"不是无穷大"时继续循环
while ~isinf(overflow_calc * 2)
    overflow_calc = overflow_calc * 2;  % 翻倍
end
% 循环结束时 overflow_calc = 2^1023
fprintf('第一阶段：找到最大2次幂 = %.16e\n', overflow_calc);

% ---- 第二阶段：二分逼近，逐步填充尾数位 ----
%
%  realmax 的二进制表示（规格化形式）：
%    1.1111...1 × 2^1023   （尾数52位全为1）
%  即：2^1023 + 2^1022 + 2^1021 + ... + 2^971
%
%  第二阶段从 2^1023 开始，step 初始为 2^1022（即当前值的一半），
%  每次尝试将 step 加入结果：
%    - 如果加上后不溢出 → 接受，overflow_calc += step
%    - 如果加上后溢出   → 拒绝，step 减半后继续尝试
%  直到 step 小于1（所有尾数位已填满）
%
step = overflow_calc / 2;  % step 初始为 2^1022

% step >= 1：当步长大于等于1时继续（步长缩小到1以下就没意义了）
while step >= 1
    if ~isinf(overflow_calc + step)
        % 加上 step 后没有溢出 → 接受这一位为1
        overflow_calc = overflow_calc + step;
    end
    % 无论接受还是拒绝，步长都减半，尝试填充下一个更低的尾数位
    step = step / 2;
end
% 循环结束后，overflow_calc 已经精确逼近到真实的 realmax

fprintf('第二阶段：二分逼近后   = %.16e\n', overflow_calc);
fprintf('【MATLAB内置】realmax  = %.16e\n', realmax);
fprintf('【理论值】(2-2^-52)*2^1023 = %.16e\n\n', (2 - 2^(-52)) * 2^1023);

% 验证精度：计算与 realmax 的相对误差
relative_error_overflow = abs(overflow_calc - realmax) / realmax;
fprintf('算法结果与realmax的相对误差 = %.2e\n\n', relative_error_overflow);


%% =========================================================
%  第四部分：单精度浮点数 (Single Precision, float32)
%  =========================================================
%
%  【数值分析】
%  IEEE 754 单精度：1位符号 + 8位指数 + 23位尾数 = 共32位（比双精度少32位）
%  理论极限值：
%    机器精度    = 2^(-23)         ≈ 1.1921e-07
%    最小正规格化 = 2^(-126)        ≈ 1.1755e-38
%    最小非规格化 = 2^(-149)        ≈ 1.4013e-45
%    最大有限数  = (2-2^(-23))*2^127 ≈ 3.4028e+38
%
%  MATLAB关键语法：single(x) 将 x 转换为单精度类型
%  单精度变量参与运算时，MATLAB保持单精度精度
%  注意：必须保证两个操作数都是 single，否则会升格为双精度！
%    正确：single(1.0) + eps_s   （两者都是single）
%    错误：1.0 + eps_s           （1.0是双精度，会把eps_s提升为双精度）
%  =========================================================

fprintf('--- 第四部分：单精度浮点数 ---\n\n');

% ---- 4-1. 单精度机器精度 ----
eps_s = single(1.0);   % single()：将字面量转为单精度，class(eps_s)返回'single'

while single(1.0) + eps_s / single(2.0) > single(1.0)
    eps_s = eps_s / single(2.0);
end

fprintf('【单精度机器精度】\n');
fprintf('  算法值             = %.8e\n',  eps_s);
fprintf('  eps(''single'')     = %.8e\n',  eps('single'));
% eps('single')：内置单精度机器精度，以双精度返回便于fprintf显示
fprintf('  理论值 2^(-23)     = %.8e\n\n', 2^(-23));

% ---- 4-2. 单精度下溢值（同样用正规化检测） ----
underflow_s = single(1.0);

% 注意：eps_s 是前面算出的单精度机器精度（single 类型）
% 条件与双精度完全一致：检测 x/2 是否仍维持相对精度 eps_s
while underflow_s / single(2) + (underflow_s / single(2)) * eps_s > underflow_s / single(2)
    underflow_s = underflow_s / single(2);
end
% 结果是最小正单精度规格化数 = 2^(-126)

fprintf('【单精度下溢值】\n');
fprintf('  算法值（最小正规格化）= %.4e\n', underflow_s);
fprintf('  realmin(''single'')   = %.4e\n', realmin('single'));
fprintf('  理论 2^(-126)        = %.4e\n\n', 2^(-126));

% ---- 4-3. 单精度上溢值（同样的两阶段算法）----
overflow_s = single(1.0);

% 第一阶段：找最大2次幂
while ~isinf(overflow_s * single(2.0))
    overflow_s = overflow_s * single(2.0);
end
% 此时 overflow_s = 2^127（单精度）

% 第二阶段：填充23位尾数
step_s = overflow_s / single(2.0);   % step_s = 2^126
while step_s >= single(1.0)
    if ~isinf(overflow_s + step_s)
        overflow_s = overflow_s + step_s;
    end
    step_s = step_s / single(2.0);
end

fprintf('【单精度上溢值】\n');
fprintf('  算法值               = %.8e\n', overflow_s);
fprintf('  realmax(''single'')   = %.8e\n', realmax('single'));
fprintf('  理论 (2-2^-23)*2^127 = %.8e\n\n', (2 - 2^(-23)) * 2^127);


%% =========================================================
%  第五部分：综合对比汇总（双精度 vs 单精度）
%  =========================================================

fprintf('==============================================\n');
fprintf('       综合对比汇总：双精度 vs 单精度\n');
fprintf('==============================================\n\n');

% 表头
fprintf('%-12s %-10s %-22s %-22s\n', '量', '精度', '算法计算值', 'MATLAB内置值');
fprintf('%s\n', repmat('-', 1, 70));
% repmat('-', 1, 70)：将字符'-'重复70次，生成分隔线

% 双精度行
fprintf('%-12s %-10s %-22.6e %-22.6e\n', '机器精度', '双精度', eps_calc,       eps);
fprintf('%-12s %-10s %-22.6e %-22.6e\n', '',         '单精度', eps_s,           eps('single'));
fprintf('%-12s %-10s %-22.6e %-22.6e\n', '下溢值',   '双精度', underflow_calc,  realmin('double'));
fprintf('%-12s %-10s %-22.6e %-22.6e\n', '',         '单精度', underflow_s,     realmin('single'));
fprintf('%-12s %-10s %-22.6e %-22.6e\n', '上溢值',   '双精度', overflow_calc,   realmax);
fprintf('%-12s %-10s %-22.6e %-22.6e\n', '',         '单精度', overflow_s,      realmax('single'));

fprintf('\n');

%% =========================================================
%  第五部分：验证这些极限的实际效果
%  =========================================================

fprintf('--- 第五部分：验证极限的实际效果 ---\n\n');

% 验证1：机器精度的效果
% 在机器精度范围内的差异，计算机无法区分
a = 1.0;
b = 1.0 + eps/2;    % 加上 eps/2（比机器精度小一半）
c = 1.0 + eps;      % 加上恰好一个 eps

fprintf('验证机器精度：\n');
fprintf('  1.0 == 1.0 + eps/2 ?  结果：%d（1=相等，0=不等）\n', a == b);
% ==：比较运算符，在MATLAB中返回1(真)或0(假)
fprintf('  1.0 == 1.0 + eps   ?  结果：%d\n\n', a == c);

% 验证2：下溢效果
fprintf('验证下溢：\n');
very_small = underflow_calc;            % 最小正浮点数
even_smaller = underflow_calc / 2;     % 再缩小一半，应该变为0
fprintf('  最小正浮点数           = %.4e\n', very_small);
fprintf('  最小正浮点数 / 2       = %.4e（下溢归零）\n\n', even_smaller);

% 验证3：上溢效果
fprintf('验证上溢：\n');
big_num = realmax;                      % 最大有限浮点数
overflow_result = realmax * 2;          % 超出范围，产生Inf
fprintf('  realmax                = %.4e\n', big_num);
fprintf('  realmax * 2            = %f（上溢变Inf）\n\n', overflow_result);
% 上溢后结果为 Inf（无穷大），%f格式会显示为 Inf

% isinf()验证
fprintf('  isinf(realmax * 2) 返回：%d（1表示确实是Inf）\n\n', isinf(overflow_result));


%% =========================================================
%  第六部分：图形可视化（可选）
%  =========================================================

fprintf('--- 第六部分：图形可视化 ---\n\n');

% figure：创建一个新的图形窗口
% 'Name'：设置图窗标题，'NumberTitle','off'：不显示默认的图号
figure('Name', '浮点数特性可视化', 'NumberTitle', 'off', 'Position', [100 100 800 500]);
% 'Position'：[左边距 下边距 宽度 高度]，单位为像素

% subplot(m, n, k)：将图窗分成 m行×n列 的网格，当前在第k个位置绘图
subplot(1, 2, 1);

% 绘制机器精度随迭代次数的变化
% 双精度需要约53次迭代（2^53次折半）；预分配60个元素避免循环内动态扩容
% zeros(1, n)：生成 1行n列 的全零行向量，预先分配内存
eps_history = zeros(1, 60);   % 预分配：避免 SAGROW 警告，提升循环性能
eps_temp = 1.0;
iter = 0;

while (1.0 + eps_temp / 2) > 1.0
    iter = iter + 1;
    eps_history(iter) = eps_temp;
    eps_temp = eps_temp / 2;
end
eps_history = eps_history(1:iter);  % 截取实际使用的部分，去掉多余的零

% semilogy：绘制 y 轴为对数刻度的折线图
% 1:length(eps_history)：从1到数组长度的整数序列，作为x轴
% 'b-o'：蓝色(b)实线(-)圆点标记(o)
semilogy(1:length(eps_history), eps_history, 'b-o', 'MarkerSize', 3);

% 添加图表元素
xlabel('迭代次数');         % x轴标签
ylabel('eps候选值（对数刻度）');  % y轴标签
title('机器精度搜索过程');   % 图标题
grid on;                    % 显示网格线

% 在图上标注最终结果
% hold on：保持当前图形，在同一坐标系上继续绘图（不清除原有内容）
hold on;
% plot：在指定位置画一个点
% length(eps_history)：最后一次迭代的x坐标
% eps_history(end)：数组最后一个元素（end是MATLAB中表示最后一个索引的关键字）
% 'r*'：红色(r)星号(*)标记，'MarkerSize',10：标记大小
plot(length(eps_history), eps_history(end), 'r*', 'MarkerSize', 10);
% 添加文字标注
% 'HorizontalAlignment','right'：文字右对齐
text(length(eps_history)-2, eps_history(end)*100, ...
    sprintf('eps=%.2e', eps_history(end)), 'Color', 'red', 'FontSize', 9);
% ...：MATLAB中的续行符，表示本行代码在下一行继续
% sprintf：将格式化字符串返回为字符串变量（不直接输出，而是返回字符串）

subplot(1, 2, 2);

% 可视化三个关键值在数轴上的位置（用对数坐标表示）
key_values = [underflow_calc, realmin, eps, realmax];
% 行向量：[a, b, c, d] 表示包含4个元素的行向量
key_labels = {'算法下溢', 'realmin', 'eps', 'realmax'};
% cell数组：用{}定义，可以存储不同类型的数据（这里存字符串）
colors = {'b', 'c', 'g', 'r'};
% 颜色代码：b=蓝, c=青, g=绿, r=红

% 用茎状图(stem)展示这些值
for i = 1:length(key_values)
    % for循环：i 从 1 递增到 key_values 的元素数量
    semilogy(i, key_values(i), 'o', 'Color', colors{i}, ...
        'MarkerSize', 12, 'MarkerFaceColor', colors{i});
    % colors{i}：用{}访问cell数组的第i个元素
    hold on;
end

% 设置x轴刻度标签
set(gca, 'XTick', 1:4, 'XTickLabel', key_labels);
% gca：get current axes，获取当前坐标轴对象
% 'XTick'：指定刻度位置，'XTickLabel'：指定刻度标签文字
% 'XTickLabel' 接受 cell 数组作为标签

ylabel('数值（对数刻度）');
title('浮点数关键极限值对比');
grid on;

% 调整图形整体标题
sgtitle('IEEE 754 双精度浮点数特性');
% sgtitle：为整个 figure 添加总标题（subplot group title）

fprintf('图形已生成，请查看弹出的图形窗口。\n\n');
fprintf('程序运行完毕。\n');
