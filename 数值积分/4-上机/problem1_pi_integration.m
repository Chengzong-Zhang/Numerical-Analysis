%% =========================================================
%  problem1_pi_integration.m
%  功能：利用数值积分公式计算圆周率 π 的近似值
%
%  【数学背景】
%  数学上已证明：
%    ∫₀¹ 4/(1+x²) dx = π
%  因为 ∫₀¹ 4/(1+x²)dx = 4·arctan(x)|₀¹ = 4·(π/4 - 0) = π
%
%  本程序分三部分：
%  (1) 复合梯形公式 + 复合 Simpson 公式，改变步长 h，分析误差随 h 的变化
%  (2) Romberg 积分法（Richardson 外推加速）
%  (3) 自适应 Simpson 积分法
%  =========================================================

clear; close all; clc;

fprintf('==============================================\n');
fprintf('   第一题：数值积分计算圆周率 π 的近似值\n');
fprintf('==============================================\n\n');

% 定义被积函数 f(x) = 4/(1+x²)
% @(x) 是 MATLAB 匿名函数语法：定义临时函数，x 是输入参数
% 4./(1+x.^2)：./和.^是逐元素运算，使函数能接受向量输入（向量化）
f = @(x) 4 ./ (1 + x.^2);

% 积分区间 [a, b]
a = 0;   % 积分下限
b = 1;   % 积分上限

% 精确值：MATLAB 内置 pi 就是双精度最精确的 π 值
pi_exact = pi;

fprintf('圆周率精确参考值：π = %.15f\n\n', pi_exact);


%% =========================================================
%  第一部分：复合梯形公式 + 复合 Simpson 公式
%  分析误差随步长 h 的变化规律
%  =========================================================
%
%  【复合梯形公式】（教材公式 4.3.1）
%    将 [a,b] 分成 n 等份，步长 h = (b-a)/n，节点 xᵢ = a + ih
%    T_n = h/2·[f(a) + f(b) + 2·∑_{i=1}^{n-1} f(xᵢ)]
%    离散误差：E_n = -(b-a)h²/12·f''(ξ)，精度阶数 O(h²)
%
%  【复合 Simpson 公式】（教材公式 4.3.5）
%    将 [a,b] 分成 2m 等份（n=2m，偶数个子区间），步长 h = (b-a)/(2m)
%    S_m = h/3·[f(a)+f(b)+4·∑奇数项+2·∑偶数项]
%    离散误差：E_m = -h⁴(b-a)/180·f⁽⁴⁾(ξ)，精度阶数 O(h⁴)
%
%  【为什么存在最优 h？】
%  总误差 = 截断误差 + 舍入误差
%  - 截断误差（离散误差）：随 h 减小而减小（梯形为 O(h²)，Simpson 为 O(h⁴)）
%  - 舍入误差：每次浮点运算引入机器精度 eps≈2.2e-16，
%    计算 n=1/h 个函数值并累加，舍入误差约为 O(eps/h)
%  - 总误差 ≈ C₁h^p + C₂·eps/h，对 h 求导令其为零，得最优 h
%  - 梯形：h_opt ~ eps^(1/3) ≈ 6e-6；Simpson：h_opt ~ eps^(1/5) ≈ 1e-3
%  =========================================================

fprintf('--- 第一部分：复合梯形 vs 复合 Simpson（误差随 h 变化）---\n\n');

% 生成不同步长 h 的测试序列（从粗到细，指数间隔）
% logspace(a, b, n)：生成从10^a到10^b的n个对数均匀间隔的数
% 这里 h 从 10^0=1 到 10^(-9)，共100个测试点
h_vals = logspace(0, -9, 100);

% 对应的子区间数 n（梯形）和 m（Simpson 需偶数个子区间）
% ceil(x)：向上取整，确保 n 不小于 1
n_vals = ceil((b - a) ./ h_vals);

% 预分配误差数组（避免循环中动态扩容）
% zeros(1, N)：创建 1×N 的零向量，数值分析中常用于预分配
err_trap = zeros(1, length(h_vals));    % 梯形误差
err_simp = zeros(1, length(h_vals));    % Simpson 误差
actual_h = zeros(1, length(h_vals));    % 实际步长

fprintf('计算中（共 %d 个步长）...\n', length(h_vals));

for k = 1:length(h_vals)
    n = n_vals(k);           % 梯形：n 等份
    m = max(1, ceil(n/2));   % Simpson：m 对，共 2m 等份（保证偶数个子区间）

    h_trap = (b - a) / n;    % 梯形实际步长
    h_simp = (b - a) / (2*m); % Simpson 实际步长
    actual_h(k) = h_trap;

    % ---- 复合梯形公式 ----
    % 生成等间距节点：linspace(a, b, n+1) 生成从 a 到 b 的 n+1 个等间距点
    x_trap = linspace(a, b, n + 1);
    % f(x_trap)：向量化调用，一次计算所有节点的函数值（利用了 ./ 和 .^）
    fx_trap = f(x_trap);
    % 复合梯形公式：T = h/2*(端点之和 + 2*内点之和)
    % sum(fx_trap(2:end-1))：对中间节点（不含首尾）求和
    T_n = h_trap / 2 * (fx_trap(1) + fx_trap(end) + 2 * sum(fx_trap(2:end-1)));
    err_trap(k) = abs(T_n - pi_exact);   % abs()：取绝对值，计算误差

    % ---- 复合 Simpson 公式 ----
    x_simp = linspace(a, b, 2*m + 1);   % 共 2m+1 个节点
    fx_simp = f(x_simp);
    % Simpson 系数：端点×1，奇数位×4（中点），偶数位（非端点）×2
    % fx_simp(2:2:end-1)：下标从2开始，步长2到倒数第二，即奇数节点（1-indexed中的偶数下标）
    % x_simp(1), x_simp(3), x_simp(5), ... 对应 x₀,x₂,x₄,... （偶数下标，内部偶数点）
    % 注意：MATLAB 下标从 1 起，节点 x₀,x₁,...,x_{2m}：
    %   端点：fx(1) 和 fx(end)
    %   中点（奇数节点，1-indexed偶数下标 2,4,...,2m）×4：fx_simp(2:2:end-1)
    %   内偶数节点（1-indexed奇数下标 3,5,...,2m-1）×2：fx_simp(3:2:end-2)
    S_m = h_simp / 3 * (fx_simp(1) + fx_simp(end) + ...
          4 * sum(fx_simp(2:2:end-1)) + ...
          2 * sum(fx_simp(3:2:end-2)));
    err_simp(k) = abs(S_m - pi_exact);
end

fprintf('计算完成。\n\n');

% ---- 分析最优步长 ----
% min(v)：返回向量 v 的最小值
% [val, idx] = min(v)：同时返回最小值 val 和其对应的下标 idx
[min_err_trap, idx_trap] = min(err_trap);
[min_err_simp, idx_simp] = min(err_simp);

fprintf('【梯形公式】最小误差 = %.4e，对应 h = %.4e，n = %d\n', ...
    min_err_trap, actual_h(idx_trap), n_vals(idx_trap));
fprintf('            理论最优 h ≈ eps^(1/3) = %.4e\n', eps^(1/3));
fprintf('【Simpson】 最小误差 = %.4e，对应 h = %.4e，n = %d\n', ...
    min_err_simp, actual_h(idx_simp), 2*ceil(n_vals(idx_simp)/2));
fprintf('            理论最优 h ≈ eps^(1/5) = %.4e\n\n', eps^(1/5));

fprintf('【结论】h 低于最优值后误差不再减小，原因：\n');
fprintf('  舍入误差 O(eps/h) 开始主导，抵消截断误差减小带来的收益。\n');
fprintf('  梯形最优 h ~ eps^(1/3)，Simpson 最优 h ~ eps^(1/5)。\n\n');

% ---- 绘图：误差随 h 的变化 ----
figure('Name', '题目1(1)：误差随步长h的变化', 'NumberTitle', 'off', ...
       'Position', [50 50 900 500]);
% loglog：双对数坐标绘图（x轴和y轴都是对数刻度）
% 在双对数坐标下，误差 O(h^p) 表现为斜率为 p 的直线
loglog(actual_h, err_trap, 'b-', 'LineWidth', 1.5, 'DisplayName', '复合梯形误差');
hold on;
% hold on：保持当前图形，继续在同一坐标系绘图（不覆盖已有内容）
loglog(actual_h, err_simp, 'r-', 'LineWidth', 1.5, 'DisplayName', '复合Simpson误差');

% 标注理论斜率参考线：梯形 O(h²)，Simpson O(h⁴)
h_ref = logspace(0, -5, 50);   % 参考 h 范围
% 调整参考线纵截距以便与实际误差曲线对齐
loglog(h_ref, 0.3*h_ref.^2, 'b--', 'LineWidth', 1.0, 'DisplayName', 'O(h²) 理论斜率');
loglog(h_ref, 0.05*h_ref.^4, 'r--', 'LineWidth', 1.0, 'DisplayName', 'O(h⁴) 理论斜率');

% 标注最优步长位置
% xline：在指定 x 值处画竖直线（参数：x值，线型，标签，属性）
xline(actual_h(idx_trap), 'b:', '梯形最优h', 'LineWidth', 1.5, 'FontSize', 9);
xline(actual_h(idx_simp), 'r:', 'Simpson最优h', 'LineWidth', 1.5, 'FontSize', 9);

% 设置坐标轴标签、标题、图例
xlabel('步长 h（对数刻度）');
ylabel('|π近似值 - π精确值|（对数刻度）');
title('复合梯形 vs 复合Simpson：误差随步长h变化的双对数图');
% legend：显示图例，'Location'：指定图例位置，'NorthEast' 表示右上角
legend('Location', 'NorthEast', 'FontSize', 9);
grid on;   % 显示网格线，方便读取数据
% axis：设置坐标轴范围，[xmin xmax ymin ymax]
axis([1e-9 1 1e-16 10]);


%% =========================================================
%  第二部分：Romberg 积分法
%  =========================================================
%
%  【数值分析背景——Romberg 积分】
%  Romberg 积分是对复合梯形公式进行 Richardson 外推加速的结果。
%
%  核心思想：
%  由 Euler-Maclaurin 公式知，n 等份复合梯形值 T_{n} 满足
%    I = T_n + c₂h² + c₄h⁴ + ...  (h = (b-a)/n)
%  当区间逐次对分时（h → h/2），误差主项按 h² 规律变化，可用
%  Richardson 外推消去 h² 项，提升收敛阶数。
%
%  【Romberg 表（教材公式 4.4.4）】
%  记 T_{m,1} 为第 m 次对分（2^{m-1} 等份）后的梯形值：
%    T_{1,1} = (b-a)/2·[f(a)+f(b)]
%    T_{m,1} = 1/2·[T_{m-1,1} + h_{m-1}·∑ f(a+(k-1/2)h_{m-1})]
%  其中 h_{m-1} = (b-a)/2^{m-2}，新增的 2^{m-2} 个中点。
%
%  然后逐列 Richardson 外推：
%    T_{m,j} = (4^{j-1}·T_{m,j-1} - T_{m-1,j-1}) / (4^{j-1} - 1)
%  对角线元素 T_{m,m} 的精度阶数为 O(h^{2m})，收敛极快。
%  =========================================================

fprintf('--- 第二部分：Romberg 积分法 ---\n\n');

% Romberg 参数
max_level = 10;   % 最大对分次数（T 表大小为 max_level × max_level）
tol_romberg = 1e-12;   % 收敛判据：相邻对角元素之差小于此值即停止

% 初始化 Romberg 表（全零矩阵）
% zeros(m, n)：创建 m 行 n 列的全零矩阵
T = zeros(max_level, max_level);

% 第一行第一列：最粗的梯形值（只用两端点，n=1，h=b-a）
% T_{1,1} = (b-a)/2·[f(a)+f(b)]
T(1, 1) = (b - a) / 2 * (f(a) + f(b));

fprintf('Romberg 表（逐行填充，显示对角元素收敛情况）：\n');
fprintf('%-5s  %-25s  %-20s\n', '对分次m', 'T_{m,m}（对角元素）', '|T_{m,m}-π|');
fprintf('%s\n', repmat('-', 1, 55));
fprintf('%5d  %-25.15f  %-20.4e\n', 1, T(1,1), abs(T(1,1) - pi_exact));

romberg_result = T(1, 1);    % 保存最终结果
converged_level = 1;          % 记录收敛时的层数

for m = 2:max_level
    % ---- 第一列：逐次对分的梯形值 ----
    % 第 m 次对分：共 2^{m-1} 个等份，步长 h_{m-1} = (b-a)/2^{m-2}
    % 新增节点为上一层节点的中点
    h_prev = (b - a) / 2^(m - 2);   % 上一层步长（第 m-1 次对分的步长）

    % 新增中点的位置：a + (k-1/2)·h_prev，k = 1, 2, ..., 2^{m-2}
    % (1:2^(m-2)) 生成从 1 到 2^{m-2} 的整数行向量
    k_vals = 1 : 2^(m - 2);
    % 新中点 x 坐标（向量）：x = a + (k-0.5)*h_prev
    x_new = a + (k_vals - 0.5) * h_prev;

    % 变步长梯形递推公式（教材 4.4.4）：
    % T_{m,1} = 1/2·[T_{m-1,1} + h_{m-1}·∑ f(新增中点)]
    h_cur = (b - a) / 2^(m - 1);   % 当前层步长
    T(m, 1) = 0.5 * T(m-1, 1) + h_cur * sum(f(x_new));
    % 解释：h_cur = h_prev/2，上式等价于：
    % 新梯形值 = 旧梯形值/2 + 新增中点的函数值之和×当前步长
    % 这避免了重复计算已有节点的函数值，节省约一半计算量

    % ---- 第 j 列：Richardson 外推 ----
    for j = 2 : m
        % Richardson 外推公式（教材 Romberg 算法）：
        % T_{m,j} = (4^{j-1}·T_{m,j-1} - T_{m-1,j-1}) / (4^{j-1} - 1)
        % 物理含义：消去 h^{2(j-1)} 误差项，使精度从 O(h^{2(j-2)}) 提升到 O(h^{2(j-1)})
        % 4^{j-1} 来自相邻梯形值步长比为 2：误差主项比为 4^{j-1}
        coeff = 4^(j - 1);   % 外推系数
        T(m, j) = (coeff * T(m, j-1) - T(m-1, j-1)) / (coeff - 1);
    end

    % 打印当前对角元素
    fprintf('%5d  %-25.15f  %-20.4e\n', m, T(m, m), abs(T(m, m) - pi_exact));

    % 收敛判断：相邻两个对角元素之差足够小
    if m >= 2 && abs(T(m, m) - T(m-1, m-1)) < tol_romberg
        romberg_result = T(m, m);
        converged_level = m;
        fprintf('\n  收敛！在第 %d 层（共 %d 个函数求值点）达到精度 %.0e\n', ...
            m, 2^(m-1) + 1, tol_romberg);
        break;
    end
    romberg_result = T(m, m);
    converged_level = m;
end

fprintf('\nRomberg 积分结果：%.15f\n', romberg_result);
fprintf('与 π 的误差：       %.4e\n', abs(romberg_result - pi_exact));
fprintf('函数求值次数：       %d 次\n\n', 2^(converged_level-1) + 1);

% 打印完整 Romberg 表（截取到收敛层）
fprintf('完整 Romberg 表（T_{m,j}，对角线收敛最快）：\n');
for row = 1:converged_level
    fprintf('m=%d: ', row);
    for col = 1:row
        fprintf('%12.8f  ', T(row, col));
    end
    fprintf('\n');
end
fprintf('\n');


%% =========================================================
%  第三部分：自适应 Simpson 积分法
%  =========================================================
%
%  【数值分析背景——自适应积分】
%  自适应积分的核心思想：
%    - 在函数变化剧烈的区域自动细化步长（加密节点）
%    - 在函数变化平缓的区域保持较大步长（减少计算）
%    - 从而在保证精度的前提下最小化函数求值次数
%
%  【自适应 Simpson 误差估计】
%  对子区间 [a,b]，比较：
%    S(a,b)     = (b-a)/6·[f(a)+4f((a+b)/2)+f(b)]   （整区间 Simpson）
%    S(a,m)+S(m,b)  = 对分后两个子区间各自的 Simpson 之和
%    m = (a+b)/2 为中点
%  由 Simpson 误差阶数 O(h⁵) 可知：
%    S(a,m)+S(m,b) 的误差 ≈ [S(a,b) 的误差]/16
%  因此局部误差估计为：
%    |误差| ≈ |S(a,m)+S(m,b) - S(a,b)| / 15
%  若此估计值 ≤ tol，则接受 S(a,m)+S(m,b)+校正项 作为最终结果；
%  否则递归对 [a,m] 和 [m,b] 各自继续细化。
%  =========================================================

fprintf('--- 第三部分：自适应 Simpson 积分法 ---\n\n');

tol_adaptive = 1e-10;    % 绝对误差容限（允许的最大误差）

% 初始的整区间 Simpson 值（递归起点）
S_initial = simp1(f, a, b);
% 调用自适应递归函数（定义在文件末尾）
[pi_adaptive, func_count] = adapt_simp(f, a, b, tol_adaptive, S_initial);

fprintf('自适应 Simpson 结果：%.15f\n', pi_adaptive);
fprintf('与 π 的误差：         %.4e\n', abs(pi_adaptive - pi_exact));
fprintf('函数求值次数：         %d 次\n', func_count);
fprintf('\n');


%% =========================================================
%  第四部分：三种方法结果汇总比较
%  =========================================================

fprintf('==============================================\n');
fprintf('         三种方法计算 π 的结果汇总\n');
fprintf('==============================================\n\n');

% 选几个有代表性的 h 值展示梯形和 Simpson 的结果
h_demo = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6];
fprintf('%-10s  %-22s  %-22s\n', 'h', '复合梯形（π近似值）', '复合Simpson（π近似值）');
fprintf('%s\n', repmat('-', 1, 58));
for h = h_demo
    n = round((b-a)/h);
    m = max(1, ceil(n/2));
    h_t = (b-a)/n;
    h_s = (b-a)/(2*m);
    x_t = linspace(a, b, n+1);
    fx_t = f(x_t);
    T_val = h_t/2*(fx_t(1)+fx_t(end)+2*sum(fx_t(2:end-1)));
    x_s = linspace(a, b, 2*m+1);
    fx_s = f(x_s);
    S_val = h_s/3*(fx_s(1)+fx_s(end)+4*sum(fx_s(2:2:end-1))+2*sum(fx_s(3:2:end-2)));
    fprintf('h=%.0e    %.15f        %.15f\n', h, T_val, S_val);
end

fprintf('\n%-25s  %-20s  %-15s\n', '方法', 'π近似值', '误差');
fprintf('%s\n', repmat('-', 1, 63));
fprintf('%-25s  %-20.15f  %-15.4e\n', 'Romberg积分法', romberg_result, abs(romberg_result-pi_exact));
fprintf('%-25s  %-20.15f  %-15.4e\n', '自适应Simpson', pi_adaptive, abs(pi_adaptive-pi_exact));
fprintf('%-25s  %-20.15f  %-15s\n', 'π精确参考值', pi_exact, '0');
fprintf('\n');


%% =========================================================
%  绘图：Romberg 收敛过程对比
%  =========================================================

figure('Name', '题目1(2-3)：Romberg与自适应Simpson', 'NumberTitle', 'off', ...
       'Position', [100 100 900 400]);

subplot(1, 2, 1);
% 绘制 Romberg 表对角线的收敛过程
romberg_diag = zeros(1, converged_level);
for m = 1:converged_level
    romberg_diag(m) = abs(T(m, m) - pi_exact);
end
% semilogy：y轴用对数刻度，x轴普通刻度
semilogy(1:converged_level, romberg_diag, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 7);
xlabel('对分次数 m');
ylabel('|T_{m,m} - π|（对数刻度）');
title('Romberg 对角线收敛过程');
grid on;

subplot(1, 2, 2);
% 绘制梯形、Simpson、Romberg、自适应的误差对比（最终精度）
methods = {'梯形(h=1e-4)', 'Simpson(h=1e-4)', 'Romberg', '自适应Simpson'};
% 计算梯形和Simpson在h=1e-4时的误差
h0 = 1e-4;
n0 = round((b-a)/h0); m0 = ceil(n0/2);
x0_t = linspace(a,b,n0+1); fx0_t = f(x0_t);
T0 = h0/2*(fx0_t(1)+fx0_t(end)+2*sum(fx0_t(2:end-1)));
x0_s = linspace(a,b,2*m0+1); fx0_s = f(x0_s);
S0 = (h0/2)/3*(fx0_s(1)+fx0_s(end)+4*sum(fx0_s(2:2:end-1))+2*sum(fx0_s(3:2:end-2)));
errs = [abs(T0-pi_exact), abs(S0-pi_exact), abs(romberg_result-pi_exact), abs(pi_adaptive-pi_exact)];

% bar：绘制条形图
bar(1:4, errs);
% set(gca, ...)：设置当前坐标轴（get current axes）的属性
set(gca, 'YScale', 'log');   % y轴对数刻度
set(gca, 'XTickLabel', methods);   % 设置 x 轴刻度标签
% xtickangle：设置 x 轴刻度标签的旋转角度（避免标签重叠）
xtickangle(20);
ylabel('误差绝对值（对数刻度）');
title('各方法精度比较');
grid on;

sgtitle('第一题：数值积分计算圆周率 π');
% sgtitle：为整个 figure 添加总标题（所有子图共用的大标题）

fprintf('图形已生成，请查看弹出的图形窗口。\n\n');
fprintf('程序运行完毕。\n');


%% =========================================================
%  辅助函数定义（MATLAB 要求辅助函数定义在文件末尾，主脚本之后）
%  =========================================================

function S = simp1(f, a, b)
    %SIMP1  对单个区间 [a,b] 应用 Simpson 公式
    %  S = (b-a)/6·[f(a) + 4f(m) + f(b)]，m=(a+b)/2
    %
    %  【数值含义】
    %  Simpson 公式用通过 a, m=(a+b)/2, b 三点的抛物线来近似被积函数，
    %  代数精度为 3（能精确积分 1, x, x², x³ 这四种函数）
    m = (a + b) / 2;          % 中点（Simpson 公式的中间节点）
    S = (b - a) / 6 * (f(a) + 4*f(m) + f(b));
    % (b-a)/6：h/3 的变形，h=(b-a)/2，故 h/3=(b-a)/6
end


function [I, cnt] = adapt_simp(f, a, b, tol, S_ab)
    %ADAPT_SIMP  自适应 Simpson 积分（递归实现）
    %
    %  输入：
    %    f    - 被积函数（匿名函数）
    %    a, b - 积分区间端点
    %    tol  - 本子区间允许的最大误差（整体误差的局部分配）
    %    S_ab - 已计算好的整区间 [a,b] 的 Simpson 值（避免重复计算）
    %
    %  输出：
    %    I   - 积分近似值
    %    cnt - 本次递归消耗的函数求值次数
    %
    %  【算法原理】
    %  将 [a,b] 从中点 m 对分，分别在 [a,m] 和 [m,b] 上算 Simpson：
    %    S(a,m) + S(m,b)  vs  S(a,b)
    %  误差估计：|误差| ≈ |S(a,m)+S(m,b) - S(a,b)| / 15
    %  （15 = 2⁴-1，来自 Simpson 误差阶 O(h⁵)，步长减半后误差减为 1/16）
    %  若误差满足容限，接受带校正的结果；否则递归继续细化。

    m = (a + b) / 2;         % 对分中点

    % 计算左半段 [a,m] 和右半段 [m,b] 的 Simpson 值
    S_am = simp1(f, a, m);   % 此处各调用 f 两次（除已知端点外），共约 3 次
    S_mb = simp1(f, m, b);   % 类似

    % 两个子区间 Simpson 值之和
    S2 = S_am + S_mb;

    % 局部误差估计（Richardson 外推估计真实误差的下界）
    err_est = abs(S2 - S_ab) / 15;
    % 分母 15 = 4^2 - 1：由于步长折半，Simpson 误差从 O(h⁵) 变为 O((h/2)⁵) = O(h⁵/32)
    % 两个子区间总误差 ≈ S_ab_err/16，故 S2 误差 ≈ (S_ab - S2)/15

    % 本次计算新增的函数求值次数（中点 m + 两个 1/4 点）
    % 注：S_am 用到 f(a), f((a+m)/2), f(m)；S_mb 用到 f(m), f((m+b)/2), f(b)
    % 其中 f(a), f(b) 由上层传入（S_ab已算），新增 f((a+m)/2), f(m), f((m+b)/2) 共3个
    cnt = 3;

    if err_est <= tol
        % 误差满足容限：接受结果，并加上 Richardson 校正项
        % 校正：S2 + (S2 - S_ab)/15 ≈ 真实积分值（消去主误差项）
        I = S2 + (S2 - S_ab) / 15;
    else
        % 误差不满足容限：递归细化左右两半区间
        % 误差平均分配给两个子区间，各允许 tol/2 的误差
        [I_left,  cnt_left]  = adapt_simp(f, a, m, tol/2, S_am);
        [I_right, cnt_right] = adapt_simp(f, m, b, tol/2, S_mb);
        I   = I_left + I_right;
        cnt = cnt + cnt_left + cnt_right;
    end
end
