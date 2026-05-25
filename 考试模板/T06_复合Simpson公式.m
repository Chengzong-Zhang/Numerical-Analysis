%% =========================================================
%  模板T06：复合 Simpson 公式
%
%  【公式】
%  将 [a,b] 分成 2m 等份（必须是偶数！），步长 h = (b-a)/(2m)
%  共 2m+1 个节点：x₀, x₁, x₂, ..., x_{2m}
%
%  系数规则：
%    端点 x₀ 和 x_{2m}：系数 1
%    奇数下标节点 x₁, x₃, ..., x_{2m-1}（中点类型）：系数 4
%    偶数内点 x₂, x₄, ..., x_{2m-2}（非端点）：系数 2
%
%    S_m = h/3 · [f(x₀) + f(x_{2m}) + 4·Σ奇 + 2·Σ偶内]
%
%  【MATLAB 下标说明】（MATLAB 下标从 1 开始！）
%  节点 x₀ → fx(1)，x₁ → fx(2)，...，x_{2m} → fx(2m+1) = fx(end)
%    奇数下标节点（x₁,x₃,...,x_{2m-1}）→ fx(2), fx(4), ..., fx(2m)   → fx(2:2:end-1)
%    偶数内点（x₂,x₄,...,x_{2m-2}）→ fx(3), fx(5), ..., fx(2m-1) → fx(3:2:end-2)
%
%  【误差】
%  截断误差 = -h⁴(b-a)/180 · f⁽⁴⁾(ξ)，精度阶数 O(h⁴)
%  h 减半，误差缩小为 1/16（比梯形快得多！）
%
%  【最优步长】
%  Simpson 最优步长：h_opt ≈ eps^(1/5) ≈ 1.4×10⁻³
%  =========================================================

clear; close all; clc;

%% ===== ① 修改这里：定义你的积分 =====
f = @(x) 4 ./ (1 + x.^2);   % 被积函数（逐元素运算！）
a = 0;                        % 积分下限
b = 1;                        % 积分上限
I_exact = pi;                 % 精确值
%% =====================================

fprintf('复合 Simpson 公式\n');
fprintf('积分区间：[%g, %g]，精确值：%.10f\n\n', a, b, I_exact);

%% ---- 用法一：给定子区间对数 m，直接计算 ----
m = 500;                          % m 对子区间（总共 2m 等份，2m 必须为偶数，自动满足）
h = (b - a) / (2 * m);           % 步长 h = (b-a)/(2m)
x = linspace(a, b, 2*m + 1);     % 生成 2m+1 个节点
% linspace(a, b, 2m+1)：共 2m+1 个均匀点，含两端点
fx = f(x);                        % 在所有节点处计算函数值

% 套公式（下标规则见上面的说明）
S = h/3 * (fx(1) + fx(end) + ...
    4 * sum(fx(2:2:end-1)) + ...   % 奇数节点（MATLAB偶数下标 2,4,...,2m）
    2 * sum(fx(3:2:end-2)));       % 偶数内点（MATLAB奇数下标 3,5,...,2m-1）
% 如果 m=1（只有3个节点），fx(3:2:end-2) 为空，sum 返回0，公式仍正确

fprintf('m=%d（h=%.2e，共%d个节点）：S_m = %.10f，误差 = %.4e\n\n', ...
    m, h, 2*m+1, S, abs(S - I_exact));

%% ---- 辅助：单区间 Simpson 公式（3点）----
% S(a,b) = (b-a)/6 · [f(a) + 4f(m) + f(b)]，m=(a+b)/2
% 注意：(b-a)/6 = h/3，其中 h=(b-a)/2（单区间步长）
m_pt = (a + b) / 2;
S_single = (b-a)/6 * (f(a) + 4*f(m_pt) + f(b));
fprintf('单区间 Simpson（仅3点）：S = %.10f，误差 = %.4e\n\n', ...
    S_single, abs(S_single - I_exact));

%% ---- 用法二：测试不同 m（观察收敛阶数）----
fprintf('--- 改变 m，观察误差随 h 的变化 ---\n');
fprintf('%-10s  %-10s  %-20s  %-15s\n', 'm', 'h', 'S_m 结果', '误差 |S_m - I|');
fprintf('%s\n', repmat('-', 1, 58));

m_test = [5, 50, 500, 5000];
for m_k = m_test
    h_k = (b - a) / (2*m_k);
    x_k = linspace(a, b, 2*m_k + 1);
    fx_k = f(x_k);
    S_k = h_k/3 * (fx_k(1) + fx_k(end) + ...
                   4*sum(fx_k(2:2:end-1)) + 2*sum(fx_k(3:2:end-2)));
    err_k = abs(S_k - I_exact);
    fprintf('%-10d  %-10.2e  %-20.15f  %-15.4e\n', m_k, h_k, S_k, err_k);
end
fprintf('\n');
fprintf('观察：m 增大10倍（h缩小10倍），误差应缩小约10000倍（因为阶数是 O(h⁴)）\n\n');

%% ---- 用法三：与梯形公式对比（误差随 h 变化的双对数图）----
h_vals = logspace(0, -7, 70);
n_trap_vals = ceil((b-a)./h_vals);
m_simp_vals = max(1, ceil(n_trap_vals/2));

err_trap = zeros(size(h_vals));
err_simp = zeros(size(h_vals));
actual_h  = zeros(size(h_vals));

for k = 1:length(h_vals)
    % 梯形
    n_k = n_trap_vals(k);
    h_t = (b-a)/n_k;
    x_t = linspace(a, b, n_k+1);
    fx_t = f(x_t);
    T_k = h_t/2*(fx_t(1)+fx_t(end)+2*sum(fx_t(2:end-1)));
    err_trap(k) = abs(T_k - I_exact);

    % Simpson
    m_k = m_simp_vals(k);
    h_s = (b-a)/(2*m_k);
    x_s = linspace(a, b, 2*m_k+1);
    fx_s = f(x_s);
    S_k = h_s/3*(fx_s(1)+fx_s(end)+4*sum(fx_s(2:2:end-1))+2*sum(fx_s(3:2:end-2)));
    err_simp(k) = abs(S_k - I_exact);

    actual_h(k) = h_t;
end

figure('Name', '梯形 vs Simpson：误差随步长变化', 'NumberTitle', 'off');
loglog(actual_h, err_trap, 'b-', 'LineWidth', 1.5, 'DisplayName', '复合梯形 O(h²)');
hold on;
loglog(actual_h, err_simp, 'r-', 'LineWidth', 1.5, 'DisplayName', '复合Simpson O(h⁴)');
% 参考斜率线
h_ref = logspace(0, -5, 50);
loglog(h_ref, 0.1*h_ref.^2, 'b--', 'LineWidth', 1, 'DisplayName', 'O(h²) 参考');
loglog(h_ref, 0.002*h_ref.^4, 'r--', 'LineWidth', 1, 'DisplayName', 'O(h⁴) 参考');
hold off;
xlabel('步长 h（对数刻度）');
ylabel('|误差|（对数刻度）');
title('复合梯形 vs 复合 Simpson：误差随步长变化（双对数图）');
legend('Location', 'NorthEast');
grid on;
axis([1e-7, 1, 1e-16, 10]);

fprintf('【关键观察】\n');
fprintf('  在双对数图上：梯形误差曲线斜率≈2，Simpson误差曲线斜率≈4\n');
fprintf('  斜率就是收敛阶数！这是验证数值方法精度阶数的标准方法\n\n');
fprintf('  Simpson 精度阶数比梯形高2阶，在同样步长下误差小得多\n');
fprintf('  但 h 超过最优步长后两者误差都会因舍入误差而增大\n');
