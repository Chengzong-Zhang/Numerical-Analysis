%% =========================================================
%  problem2_planck_integral.m
%  功能：用各种数值积分方法计算 Planck 黑体辐射积分
%
%  【物理背景】
%  Planck 黑体辐射理论（1900年）推导出单位频率的辐射能量密度：
%    B(ν) ∝ ν³ / (e^{hν/kT} - 1)
%
%  对所有频率积分后，总辐射能量正比于 Stefan-Boltzmann 常数，
%  其数学核心是以下无穷积分：
%    I = ∫₀^{+∞} x³/(eˣ+1) dx
%
%  注意：物理上 Stefan-Boltzmann 定律用的是 eˣ-1（玻色子），
%  这里题目给的是 eˣ+1（费米子型），精确值为：
%    I = 7π⁴/120 ≈ 5.68219...
%
%  (Fermi-Dirac 分布：∫₀^∞ x^{n-1}/(eˣ+1)dx = (1-2^{1-n})Γ(n)ζ(n))
%  n=4：I = (1-2^{-3})·6·π⁴/90 = (7/8)·(6π⁴/90) = 7π⁴/120
%
%  【数值挑战】
%  1. 积分上限为 +∞，需截断处理（当 x 很大时 x³·e^{-x}→0）
%  2. 被积函数在 x→0 时有 0/0 型不定式（需要 L'Hôpital，极限为 0）
%  3. 函数峰值在 x≈2.8 附近，衰减较慢，需选取足够大的上限
%  =========================================================

clear; close all; clc;

fprintf('==============================================\n');
fprintf('   第二题：Planck 黑体辐射积分\n');
fprintf('   I = ∫₀^{+∞} x³/(eˣ+1) dx\n');
fprintf('==============================================\n\n');

% 被积函数：f(x) = x³/(eˣ+1)
% 注意 x=0 时：x³=0 且 eˣ+1=2，故极限值为 0/2=0，无奇点
% 数值上当 x 很大时 eˣ 极大，浮点会溢出，需判断：
%   当 eˣ > 1/eps（大约 x > 37），直接用 x³·e^{-x} 近似（避免 Inf）
f_planck = @(x) planck_func(x);   % 调用稳健的 Planck 函数（定义在文件末尾）

% 精确值：7π⁴/120（由 Fermi-Dirac 积分解析公式）
I_exact = 7 * pi^4 / 120;
fprintf('精确值：I = 7π⁴/120 = %.10f\n\n', I_exact);

% 积分截断上限：选取足够大使尾部可忽略
% 当 x=50 时，f(50) = 50³·e^{-50} ≈ 1.25e5·1.93e-22 ≈ 2.4e-17 ≈ 0（极小）
% 积分尾部 ∫_{X}^∞ ≤ ∫_{X}^∞ x³e^{-x}dx → 0 极快
x_min = 1e-10;   % 避免 x=0（实际上 f(0)=0，数值上用极小正数）
x_max = 50;      % 截断上限

fprintf('积分区间截断为 [%.0e, %d]（实际误差可忽略不计）\n\n', x_min, x_max);

% 显示函数在关键点的值（验证代码正确性）
fprintf('函数值验证：\n');
check_x = [0.01, 1, 2.8, 5, 10, 20, 50];
for x = check_x
    fprintf('  f(%.2f) = %.6e\n', x, f_planck(x));
end
fprintf('\n');


%% =========================================================
%  方法一：复合梯形公式
%  =========================================================
%
%  【复合梯形公式回顾】（教材公式 4.3.1）
%  T_n = h/2·[f(a)+f(b)+2∑f(xᵢ)]，h=(b-a)/n，xᵢ = a+ih
%  离散误差 O(h²)
%  =========================================================

fprintf('--- 方法一：复合梯形公式 ---\n\n');

n_trap = 10000;   % 子区间数（步长越小精度越高，但受舍入误差限制）
h_trap = (x_max - x_min) / n_trap;   % 步长
x_trap = linspace(x_min, x_max, n_trap + 1);   % 等间距节点
fx_trap = arrayfun(f_planck, x_trap);
% arrayfun(func, array)：对数组每个元素依次调用 func，等价于循环
% 这里用 arrayfun 是因为 planck_func 内部有 if 语句，不支持向量整体输入

% 复合梯形求和
I_trap = h_trap / 2 * (fx_trap(1) + fx_trap(end) + 2*sum(fx_trap(2:end-1)));

fprintf('复合梯形（n=%d）：I = %.10f，误差 = %.4e\n\n', ...
    n_trap, I_trap, abs(I_trap - I_exact));


%% =========================================================
%  方法二：复合 Simpson 公式
%  =========================================================
%
%  【复合 Simpson 公式回顾】（教材公式 4.3.5）
%  S_m = h/3·[f(a)+f(b)+4∑奇数位+2∑偶数内点]
%  离散误差 O(h⁴)，比梯形精度高两阶
%  =========================================================

fprintf('--- 方法二：复合 Simpson 公式 ---\n\n');

m_simp = 5000;   % m 对子区间（总共 2m 个等份）
h_simp = (x_max - x_min) / (2 * m_simp);   % 步长 h = (b-a)/(2m)
x_simp = linspace(x_min, x_max, 2*m_simp + 1);   % 2m+1 个节点
fx_simp = arrayfun(f_planck, x_simp);

% 复合 Simpson 公式
% 奇数下标（1-indexed 2,4,...,2m）对应数学上的"中点"×4
% 偶数下标（1-indexed 3,5,...,2m-1）对应数学上的"内偶数点"×2
I_simp = h_simp / 3 * (fx_simp(1) + fx_simp(end) + ...
         4*sum(fx_simp(2:2:end-1)) + ...
         2*sum(fx_simp(3:2:end-2)));

fprintf('复合Simpson（m=%d,共%d点）：I = %.10f，误差 = %.4e\n\n', ...
    m_simp, 2*m_simp+1, I_simp, abs(I_simp - I_exact));


%% =========================================================
%  方法三：Romberg 积分法
%  =========================================================
%
%  【Romberg 对非有界积分的处理】
%  对截断区间 [x_min, x_max] 应用 Romberg 积分。
%  Romberg 的优势：无需预先指定步长，自动达到高精度，
%  且使用的函数求值次数远少于同精度的简单复合规则。
%  =========================================================

fprintf('--- 方法三：Romberg 积分法 ---\n\n');

max_lv = 12;       % 最大对分层数
tol_rom = 1e-8;    % 收敛阈值
T_rom = zeros(max_lv, max_lv);

% 第一层：最粗的梯形（仅用两端点）
T_rom(1, 1) = (x_max - x_min) / 2 * (f_planck(x_min) + f_planck(x_max));

fprintf('Romberg 表对角线收敛过程：\n');
fprintf('%-5s  %-22s  %-15s\n', 'm', 'T_{m,m}', '|误差|');
fprintf('%s\n', repmat('-', 1, 45));
fprintf('%5d  %-22.12f  %-15.4e\n', 1, T_rom(1,1), abs(T_rom(1,1)-I_exact));

I_romberg = T_rom(1,1);
rom_level = 1;

for m = 2:max_lv
    % 变步长梯形递推：新增中点
    h_rom_prev = (x_max - x_min) / 2^(m-2);   % 上一层步长
    h_rom_cur  = (x_max - x_min) / 2^(m-1);   % 当前层步长
    k_new = 1:2^(m-2);                          % 新增点编号
    x_new = x_min + (k_new - 0.5) * h_rom_prev;% 新增中点坐标

    % T_{m,1}：变步长梯形递推
    T_rom(m, 1) = 0.5 * T_rom(m-1, 1) + h_rom_cur * sum(arrayfun(f_planck, x_new));

    % Richardson 外推填充各列
    for j = 2:m
        c = 4^(j-1);
        T_rom(m, j) = (c * T_rom(m,j-1) - T_rom(m-1,j-1)) / (c - 1);
    end

    fprintf('%5d  %-22.12f  %-15.4e\n', m, T_rom(m,m), abs(T_rom(m,m)-I_exact));

    if abs(T_rom(m,m) - T_rom(m-1,m-1)) < tol_rom
        I_romberg = T_rom(m, m);
        rom_level = m;
        fprintf('  收敛于第 %d 层\n', m);
        break;
    end
    I_romberg = T_rom(m, m);
    rom_level = m;
end

fprintf('\nRomberg 结果：I = %.10f，误差 = %.4e\n\n', ...
    I_romberg, abs(I_romberg - I_exact));


%% =========================================================
%  方法四：自适应 Simpson 积分法
%  =========================================================
%
%  自适应积分最适合处理此类"函数在某区域集中"的情形：
%  Planck 函数在 x≈2.8 附近有峰值，在两侧迅速衰减，
%  自适应算法会在峰值附近自动加密节点，在尾部保持稀疏，
%  从而以最少的函数求值次数达到目标精度。
%  =========================================================

fprintf('--- 方法四：自适应 Simpson 积分法 ---\n\n');

tol_adapt = 1e-8;
S_init = simp_one(f_planck, x_min, x_max);
[I_adapt, fcnt_adapt] = adaptive_simp(f_planck, x_min, x_max, tol_adapt, S_init);

fprintf('自适应Simpson结果：I = %.10f，误差 = %.4e，函数求值 %d 次\n\n', ...
    I_adapt, abs(I_adapt - I_exact), fcnt_adapt + 3);


%% =========================================================
%  方法五：MATLAB 内置自适应积分（作为参考基准）
%  =========================================================
%
%  integral(fun, a, b)：MATLAB 内置的高精度自适应积分函数
%  使用 Gauss-Kronrod 规则，能处理广义积分（可用 Inf 作为上限）
%  这里用 Inf 作上限，直接计算原始广义积分，作为精度验证基准
%  =========================================================

fprintf('--- 方法五：MATLAB 内置积分（基准）---\n\n');

% integral(fun, 0, Inf)：直接计算广义积分，MATLAB 内部自动处理 Inf
% 'RelTol'：相对误差容限；'AbsTol'：绝对误差容限
I_matlab = integral(f_planck, 0, Inf, 'RelTol', 1e-12, 'AbsTol', 1e-14);
fprintf('MATLAB integral(0,Inf)：I = %.10f，误差 = %.4e\n\n', ...
    I_matlab, abs(I_matlab - I_exact));


%% =========================================================
%  方法六：变量替换 + Gauss-Laguerre 型积分（选做加分）
%  =========================================================
%
%  【思路】
%  对无穷积分，令 t = e^{-x}，则 x = -ln(t), dx = -1/t dt
%  当 x: 0→∞ 时，t: 1→0
%  I = ∫₀^1 (-ln t)³/(1/t + 1) · (1/t) dt
%    = ∫₀^1 (-ln t)³/(1 + t) dt
%  这样将无穷积分转化为 [0,1] 上的有界积分，再用复合 Simpson 求解
%  =========================================================

fprintf('--- 方法六：变量替换法（∫转化为[0,1]上的积分）---\n\n');

% 变换后的被积函数：g(t) = (-ln t)³ / (1 + t)，t∈(0,1)
% t→0+ 时：(-ln t)³→∞，但 (1+t)→1，g(t)∼(-ln t)³→∞ 很慢
% t→1- 时：(-ln t)→0，g(t)→0
% 注意 t=0 是奇点（ln0=-∞），需从极小正数开始
g_transform = @(t) (-log(t)).^3 ./ (1 + t);

% 在 [eps^0.1, 1-eps^0.1] 上用复合 Simpson（避免端点奇点）
a_g = 1e-12;   % 极小正数避开 t=0
b_g = 1 - 1e-12;
m_g = 50000;   % 需要较多节点（因为 t→0 时函数变化剧烈）
h_g = (b_g - a_g) / (2*m_g);
x_g = linspace(a_g, b_g, 2*m_g + 1);
fx_g = g_transform(x_g);
I_transform = h_g / 3 * (fx_g(1) + fx_g(end) + ...
              4*sum(fx_g(2:2:end-1)) + 2*sum(fx_g(3:2:end-2)));

fprintf('变量替换法（m=%d）：I = %.10f，误差 = %.4e\n\n', ...
    m_g, I_transform, abs(I_transform - I_exact));


%% =========================================================
%  综合比较与可视化
%  =========================================================

fprintf('==============================================\n');
fprintf('          六种方法计算结果综合比较\n');
fprintf('==============================================\n\n');

methods_name = {'复合梯形(n=10000)', '复合Simpson(m=5000)', 'Romberg', ...
                '自适应Simpson', 'MATLAB integral', '变量替换+Simpson'};
results = [I_trap, I_simp, I_romberg, I_adapt, I_matlab, I_transform];

fprintf('%-28s  %-18s  %-15s\n', '方法', '计算值', '误差');
fprintf('%s\n', repmat('-', 1, 65));
for k = 1:length(methods_name)
    fprintf('%-28s  %-18.10f  %-15.4e\n', methods_name{k}, results(k), abs(results(k)-I_exact));
end
fprintf('%-28s  %-18.10f  %-15s\n', '精确值 7π⁴/120', I_exact, '0');

% ---- 绘图 ----
figure('Name', '第二题：Planck 黑体辐射积分', 'NumberTitle', 'off', ...
       'Position', [50 50 1000 600]);

% 子图1：被积函数图像
subplot(2, 2, 1);
x_plot = linspace(0.01, 15, 500);   % 绘图用的 x 点序列
y_plot = arrayfun(f_planck, x_plot);
% fill：用于绘制填充多边形，这里填充曲线下面积以直观显示积分区域
fill([x_plot, fliplr(x_plot)], [y_plot, zeros(1,length(y_plot))], ...
     [0.7 0.85 1], 'EdgeColor', 'b', 'LineWidth', 1.2);
% fliplr(v)：将向量 v 左右翻转，用于构造填充多边形的闭合路径
xlabel('x');
ylabel('f(x) = x³/(eˣ+1)');
title('Planck 被积函数（阴影为积分区域）');
grid on;

% 子图2：各方法误差对比（条形图）
subplot(2, 2, 2);
errs_all = abs(results - I_exact);
bar(errs_all);
set(gca, 'YScale', 'log');
set(gca, 'XTickLabel', {'梯形', 'Simpson', 'Romberg', '自适应', '内置', '变换'});
xtickangle(30);
ylabel('误差（对数刻度）');
title('六种方法精度比较');
grid on;

% 子图3：Romberg 表对角线收敛
subplot(2, 2, 3);
romberg_err_diag = zeros(1, rom_level);
for m = 1:rom_level
    romberg_err_diag(m) = abs(T_rom(m,m) - I_exact);
end
semilogy(1:rom_level, romberg_err_diag, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 7);
xlabel('对分层数 m');
ylabel('|T_{m,m} - I_{exact}|');
title('Romberg 积分收敛过程');
grid on;

% 子图4：复合梯形不同步长下的误差
subplot(2, 2, 4);
n_range = round(logspace(1, 4, 30));   % 子区间数从 10 到 10000
err_range = zeros(1, length(n_range));
for k = 1:length(n_range)
    n_k = n_range(k);
    h_k = (x_max - x_min) / n_k;
    x_k = linspace(x_min, x_max, n_k+1);
    fx_k = arrayfun(f_planck, x_k);
    T_k = h_k/2*(fx_k(1)+fx_k(end)+2*sum(fx_k(2:end-1)));
    err_range(k) = abs(T_k - I_exact);
end
loglog(n_range, err_range, 'b-o', 'MarkerSize', 4, 'LineWidth', 1.2);
hold on;
% 理论 O(1/n²) 参考线
loglog(n_range, 10./n_range.^2, 'r--', 'LineWidth', 1.2, 'DisplayName', 'O(n^{-2})');
xlabel('子区间数 n');
ylabel('误差');
title('梯形公式收敛阶数验证');
legend({'实际误差', 'O(n^{-2})'}, 'Location', 'NorthEast');
grid on;

sgtitle('第二题：Planck 黑体辐射积分 ∫₀^∞ x³/(eˣ+1)dx = 7π⁴/120');

fprintf('\n图形已生成，请查看弹出的图形窗口。\n');
fprintf('程序运行完毕。\n');


%% =========================================================
%  辅助函数定义
%  =========================================================

function y = planck_func(x)
    %PLANCK_FUNC  稳健计算 f(x) = x³/(eˣ+1)
    %
    %  数值挑战：
    %  - x=0：0/0 型，极限为 0，直接赋值
    %  - x 很大：exp(x) 溢出为 Inf，用 x³·e^{-x} 近似（误差 ≈ e^{-x} 极小）
    %  - 正常情况：直接计算

    if x <= 0
        y = 0;             % f(0) = 0（极限值），x<0 物理上不用（但防御性编程）
    elseif x > 500
        y = 0;             % f(x)≈0（极快衰减，可忽略）
    elseif x > 37          % exp(x) > 1/eps，精度丧失，改用近似
        % 当 x 很大时：x³/(eˣ+1) ≈ x³·e^{-x}·1/(1+e^{-x}) ≈ x³·e^{-x}
        y = x^3 * exp(-x);
    else
        y = x^3 / (exp(x) + 1);   % 正常计算
    end
end


function S = simp_one(f, a, b)
    %SIMP_ONE  单区间 Simpson 公式：(b-a)/6·[f(a)+4f(m)+f(b)]
    m = (a + b) / 2;
    S = (b - a) / 6 * (f(a) + 4*f(m) + f(b));
end


function [I, cnt] = adaptive_simp(f, a, b, tol, S_ab)
    %ADAPTIVE_SIMP  自适应 Simpson 积分（递归）
    %  原理同第一题，此处处理更一般的 [a,b] 区间
    m = (a + b) / 2;
    S_am = simp_one(f, a, m);
    S_mb = simp_one(f, m, b);
    S2   = S_am + S_mb;
    cnt  = 3;   % 新增：f((a+m)/2), f(m), f((m+b)/2)

    if abs(S2 - S_ab) / 15 <= tol
        I = S2 + (S2 - S_ab) / 15;
    else
        [I_l, c_l] = adaptive_simp(f, a, m, tol/2, S_am);
        [I_r, c_r] = adaptive_simp(f, m, b, tol/2, S_mb);
        I   = I_l + I_r;
        cnt = cnt + c_l + c_r;
    end
end
