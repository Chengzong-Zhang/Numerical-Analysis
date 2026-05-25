%% =========================================================
%  模板T09：广义积分（无穷上限）的数值处理
%
%  【问题描述】
%  计算 ∫₀^{+∞} f(x) dx（无穷积分）
%  例：∫₀^∞ x³/(eˣ+1) dx = 7π⁴/120 ≈ 5.68219
%
%  【主要策略】
%  方法1：截断上限——选取足够大的 X，使尾部 ∫_X^∞ 可忽略
%  方法2：变量替换——令 t=e^{-x}，把 [0,∞) 映射到 [0,1]
%  方法3：MATLAB 内置 integral（直接传 Inf 作上限）
%
%  【数值陷阱处理】
%  (a) x→0 时 0/0 型不定式：极限赋值
%  (b) x 很大时 exp(x) 溢出（超过约 x=709）：改用 x³·e^{-x} 近似
%  =========================================================

clear; close all; clc;

%% ===== ① 修改这里：定义被积函数和精确值 =====
% 被积函数：f(x) = x³/(eˣ+1)，x∈[0,+∞)
% 使用稳健版本（处理数值奇点）
f_planck = @(x) planck_func(x);   % 稳健函数见文件末尾

I_exact = 7 * pi^4 / 120;   % 精确值：7π⁴/120（Fermi-Dirac 积分解析值）
fprintf('Planck 积分：∫₀^∞ x³/(eˣ+1) dx\n');
fprintf('精确值 = 7π⁴/120 = %.10f\n\n', I_exact);
%% =====================================================

%% ===================================================
%  方法1：截断上限 + 复合梯形
%  选择截断点 X 使得 ∫_X^∞ f dx ≈ 0
%  ===================================================
fprintf('--- 方法1：截断上限 + 复合梯形 ---\n');

% 选择截断点：当 x=50 时，f(50) = 50³·e⁻⁵⁰ ≈ 10⁻¹⁷ ≈ 0（可忽略）
x_min = 1e-10;   % 避开 x=0（实际上 f(0)=0，用极小正数安全）
x_max = 50;      % 截断上限（根据函数衰减速度选择）
n_trap = 10000;  % 子区间数

h = (x_max - x_min) / n_trap;
x_trap = linspace(x_min, x_max, n_trap + 1);
% 注意：planck_func 内部有 if 语句，不能整体向量化，需用 arrayfun
fx_trap = arrayfun(f_planck, x_trap);
% arrayfun(func, array)：对数组每个元素依次调用 func（处理非向量化函数）

I_trap = h/2 * (fx_trap(1) + fx_trap(end) + 2*sum(fx_trap(2:end-1)));
fprintf('复合梯形（n=%d，区间[%.0e,%d]）：\n', n_trap, x_min, x_max);
fprintf('  结果 = %.10f，误差 = %.4e\n\n', I_trap, abs(I_trap - I_exact));

%% ===================================================
%  方法2：截断上限 + 复合 Simpson
%  ===================================================
fprintf('--- 方法2：截断上限 + 复合 Simpson ---\n');
m_simp = 5000;   % m 对子区间（总共 2m 等份）
h_s = (x_max - x_min) / (2*m_simp);
x_simp = linspace(x_min, x_max, 2*m_simp + 1);
fx_simp = arrayfun(f_planck, x_simp);
I_simp = h_s/3 * (fx_simp(1) + fx_simp(end) + ...
         4*sum(fx_simp(2:2:end-1)) + 2*sum(fx_simp(3:2:end-2)));
fprintf('复合Simpson（m=%d）：\n', m_simp);
fprintf('  结果 = %.10f，误差 = %.4e\n\n', I_simp, abs(I_simp - I_exact));

%% ===================================================
%  方法3：截断上限 + Romberg
%  ===================================================
fprintf('--- 方法3：截断上限 + Romberg ---\n');
max_lv = 12; tol_rom = 1e-8;
T_rom = zeros(max_lv, max_lv);
T_rom(1,1) = (x_max-x_min)/2 * (f_planck(x_min) + f_planck(x_max));

I_romberg = T_rom(1,1);
for m = 2:max_lv
    h_p = (x_max-x_min)/2^(m-2);
    h_c = (x_max-x_min)/2^(m-1);
    k_new = 1:2^(m-2);
    x_new = x_min + (k_new - 0.5)*h_p;
    T_rom(m,1) = 0.5*T_rom(m-1,1) + h_c*sum(arrayfun(f_planck, x_new));
    for j = 2:m
        c4 = 4^(j-1);
        T_rom(m,j) = (c4*T_rom(m,j-1) - T_rom(m-1,j-1))/(c4-1);
    end
    I_romberg = T_rom(m,m);
    if abs(T_rom(m,m)-T_rom(m-1,m-1)) < tol_rom
        fprintf('  Romberg 收敛于第 %d 层\n', m);
        break;
    end
end
fprintf('Romberg（截断区间[%.0e,%d]）：\n', x_min, x_max);
fprintf('  结果 = %.10f，误差 = %.4e\n\n', I_romberg, abs(I_romberg - I_exact));

%% ===================================================
%  方法4：MATLAB 内置 integral（直接算到 Inf）
%  这是最简便的参考方法，考试中可用来验证其他方法
%  ===================================================
fprintf('--- 方法4：MATLAB integral(0, Inf) ---\n');
% integral(fun, 0, Inf)：MATLAB 内置自适应积分，可直接用 Inf 作上限
% 'RelTol'：相对误差容限；'AbsTol'：绝对误差容限
I_matlab = integral(f_planck, 0, Inf, 'RelTol', 1e-12, 'AbsTol', 1e-14);
fprintf('MATLAB integral：\n');
fprintf('  结果 = %.10f，误差 = %.4e\n\n', I_matlab, abs(I_matlab - I_exact));

%% ===================================================
%  方法5：变量替换法（将无穷积分变为有界积分）
%
%  令 t = e^{-x}，则 x = -ln(t), dx = -1/t·dt
%  当 x: 0→∞ 时，t: 1→0
%  ∫₀^∞ x³/(eˣ+1) dx = ∫₀¹ (-ln t)³/(1/t + 1) · (1/t) dt
%                      = ∫₀¹ (-ln t)³/(1 + t) dt
%
%  这把无穷积分变成了 [0,1] 上的有界积分！
%  ===================================================
fprintf('--- 方法5：变量替换法（t=e^{-x}，映射到[0,1]）---\n');

% 变换后的被积函数：g(t) = (-ln t)³ / (1+t)，t∈(0,1)
% t→0+ 时有奇点（-ln t → ∞），但可积，需从极小正数开始
g_transform = @(t) (-log(t)).^3 ./ (1 + t);
% -log(t)：MATLAB 的 log 是自然对数 ln；t→0 时 -log(t)→∞

a_g = 1e-12;   % 避开 t=0 的奇点（从极小正数开始）
b_g = 1 - 1e-12;
m_g = 50000;   % 需要较多节点（t→0 时函数变化剧烈）
h_g = (b_g - a_g) / (2*m_g);
x_g = linspace(a_g, b_g, 2*m_g + 1);
fx_g = g_transform(x_g);   % 此函数支持向量化（.^ 和 ./ 都是逐元素）
I_transform = h_g/3 * (fx_g(1)+fx_g(end) + ...
              4*sum(fx_g(2:2:end-1)) + 2*sum(fx_g(3:2:end-2)));
fprintf('变量替换 + 复合Simpson（m=%d）：\n', m_g);
fprintf('  结果 = %.10f，误差 = %.4e\n\n', I_transform, abs(I_transform - I_exact));

%% ---- 汇总对比 ----
fprintf('==============================================\n');
fprintf('         五种方法结果汇总\n');
fprintf('==============================================\n');
methods = {'截断+梯形', '截断+Simpson', '截断+Romberg', 'MATLAB integral', '变量替换'};
results_val = [I_trap, I_simp, I_romberg, I_matlab, I_transform];
fprintf('%-20s  %-18s  %-15s\n', '方法', '计算结果', '误差');
fprintf('%s\n', repmat('-', 1, 57));
for k = 1:length(methods)
    fprintf('%-20s  %-18.10f  %-15.4e\n', methods{k}, results_val(k), abs(results_val(k)-I_exact));
end
fprintf('%-20s  %-18.10f  %-15s\n', '精确值', I_exact, '0');

%% =========================================================
%  局部函数
%  =========================================================

function y = planck_func(x)
    % 稳健计算 f(x) = x³/(eˣ+1)，处理数值奇点
    %
    % 三种情况：
    % (a) x ≤ 0：极限值 f(0) = 0（0/0 型，分子→0 更快）
    % (b) x > 500：函数值 ≈ 0（极快衰减，可忽略）
    % (c) x > 37：exp(x) 将溢出（超过双精度最大值约 1.8e308）
    %             改用近似：x³/(eˣ+1) ≈ x³·e^{-x}（误差 ≈ e^{-x} 极小）
    % (d) 正常情况：直接计算
    if x <= 0
        y = 0;             % f(0)=0（极限），x<0 防御性处理
    elseif x > 500
        y = 0;             % 完全可忽略
    elseif x > 37         % exp(x) 将超出双精度范围
        y = x^3 * exp(-x);  % x³/(eˣ+1) ≈ x³/eˣ = x³·e^{-x}（相对误差≈e^{-x}）
    else
        y = x^3 / (exp(x) + 1);   % 正常计算（注意是 eˣ+1，费米子型）
    end
end
