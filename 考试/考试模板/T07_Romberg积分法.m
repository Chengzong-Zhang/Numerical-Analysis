%% =========================================================
%  模板T07：Romberg 积分法（Richardson 外推加速）
%
%  【核心思想】
%  对复合梯形公式进行逐次对分，然后用 Richardson 外推消去误差主项，
%  使精度从 O(h²) 迅速提升到 O(h^{2m})（第 m 层）。
%
%  【Romberg 表构造】
%  第一列 T_{m,1}：第 m 次对分后的梯形值（子区间数 = 2^{m-1}）
%  第 j 列（j≥2）：Richardson 外推
%
%  T_{m,j} = [4^{j-1} · T_{m,j-1} - T_{m-1,j-1}] / (4^{j-1} - 1)
%
%  对角线 T_{m,m} 的精度为 O(h^{2m})，收敛极快。
%
%  【变步长梯形递推（节省计算）】
%  T_{m,1} = (1/2)·T_{m-1,1} + h_m·Σ f(新增中点)
%  每次对分只需计算新增的 2^{m-2} 个中点，不重复计算已有节点。
%
%  【Romberg 表示意图（以4层为例）】
%  T_{1,1}
%  T_{2,1}  T_{2,2}
%  T_{3,1}  T_{3,2}  T_{3,3}
%  T_{4,1}  T_{4,2}  T_{4,3}  T_{4,4}  ← 对角线最精确
%
%  阶数：第1列 O(h²)，第2列 O(h⁴)，第3列 O(h⁶)，...
%  =========================================================

clear; close all; clc;

%% ===== ① 修改这里：定义你的积分 =====
f = @(x) 4 ./ (1 + x.^2);   % 被积函数（逐元素运算！）
a = 0;                        % 积分下限
b = 1;                        % 积分上限
I_exact = pi;                 % 精确值
%% =====================================

max_level = 10;     % 最大对分次数（Romberg 表最大层数）
tol = 1e-12;        % 收敛判据：相邻对角元素之差小于此值即停止

fprintf('Romberg 积分法\n');
fprintf('积分区间：[%g, %g]，精确值：%.15f\n\n', a, b, I_exact);

%% ---- 第一步：初始化 Romberg 表 ----
T = zeros(max_level, max_level);   % 预分配 Romberg 表（全零矩阵）
% T(m, j) 存储 T_{m,j}

% 第 1 层第 1 列：最粗的梯形（只用两端点，n=1，h=b-a）
% T_{1,1} = (b-a)/2 · [f(a) + f(b)]
T(1, 1) = (b - a) / 2 * (f(a) + f(b));

fprintf('Romberg 表（逐层填充，只显示对角线）：\n');
fprintf('%-6s  %-25s  %-20s\n', '层数 m', 'T_{m,m}（对角元素）', '|T_{m,m} - 精确值|');
fprintf('%s\n', repmat('-', 1, 55));
fprintf('%6d  %-25.15f  %-20.4e\n', 1, T(1,1), abs(T(1,1) - I_exact));

%% ---- 第二步：逐层填充 Romberg 表 ----
converged_level = 1;
for m = 2 : max_level

    %% -- 2a：计算第一列 T_{m,1}（变步长梯形递推）--
    % 上一层步长 h_{m-1} = (b-a) / 2^{m-2}
    h_prev = (b - a) / 2^(m - 2);
    % 当前层步长 h_m = (b-a) / 2^{m-1}  （= h_{m-1}/2）
    h_cur  = (b - a) / 2^(m - 1);

    % 新增中点的编号：k = 1, 2, ..., 2^{m-2}
    % 新增中点坐标：x = a + (k - 0.5) · h_{m-1}，即每个旧区间的中点
    k_new = 1 : 2^(m-2);
    x_new = a + (k_new - 0.5) * h_prev;
    % (k_new - 0.5)：将 k=1,2,...,2^{m-2} 映射到 0.5, 1.5, 2.5,...

    % 变步长递推公式：
    % T_{m,1} = (1/2)·T_{m-1,1} + h_m · Σ f(新增中点)
    % 解释：旧梯形值/2（端点贡献折半）+ 新增中点的贡献
    T(m, 1) = 0.5 * T(m-1, 1) + h_cur * sum(f(x_new));
    % sum(f(x_new))：对向量 x_new 的所有元素求 f 值并求和（向量化）

    %% -- 2b：Richardson 外推填充第 2 到第 m 列 --
    for j = 2 : m
        coeff = 4^(j - 1);   % 外推系数 4^{j-1}
        % Richardson 外推公式：消去 h^{2(j-1)} 误差项
        % T_{m,j} = [4^{j-1}·T_{m,j-1} - T_{m-1,j-1}] / (4^{j-1} - 1)
        T(m, j) = (coeff * T(m, j-1) - T(m-1, j-1)) / (coeff - 1);
        % 第 j 列的精度比第 j-1 列高 2 阶
    end

    fprintf('%6d  %-25.15f  %-20.4e\n', m, T(m, m), abs(T(m, m) - I_exact));

    %% -- 2c：收敛判断 --
    if abs(T(m, m) - T(m-1, m-1)) < tol
        fprintf('\n  收敛！在第 %d 层达到精度 %.0e\n', m, tol);
        fprintf('  本次共调用函数 %d 次（新增中点累计）\n', 2^(m-1) + 1);
        converged_level = m;
        break;
    end
    converged_level = m;
end

%% ---- 第三步：输出结果 ----
fprintf('\nRomberg 积分结果：\n');
fprintf('  近似值 = %.15f\n', T(converged_level, converged_level));
fprintf('  精确值 = %.15f\n', I_exact);
fprintf('  误  差 = %.4e\n\n', abs(T(converged_level, converged_level) - I_exact));

%% ---- 打印完整 Romberg 表 ----
fprintf('完整 Romberg 表（T_{m,j}，行=对分层数，列=外推次数）：\n');
fprintf('（对角线 T_{m,m} 收敛最快，每列精度提升2阶）\n');
for row = 1 : converged_level
    fprintf('m=%d: ', row);
    for col = 1 : row
        fprintf('%15.10f  ', T(row, col));
    end
    fprintf('\n');
end
fprintf('\n');

%% ---- 绘图：对角线收敛过程 ----
diag_errs = zeros(1, converged_level);
for m = 1 : converged_level
    diag_errs(m) = abs(T(m, m) - I_exact);
end

figure('Name', 'Romberg 积分收敛过程', 'NumberTitle', 'off');
semilogy(1:converged_level, diag_errs, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8);
% semilogy：x轴普通刻度，y轴对数刻度
xlabel('对分层数 m');
ylabel('|T_{m,m} - 精确值|（对数刻度）');
title('Romberg 对角线收敛过程（精度随层数指数级提升）');
grid on;

fprintf('【关键观察】\n');
fprintf('  Romberg 只需很少层（约8-10层）就能达到双精度极限（1e-15）\n');
fprintf('  同等精度下，Romberg 的函数求值次数远少于固定步长的梯形/Simpson\n');
fprintf('  精度提升极快：每增加一层，误差大约减小 10^3 到 10^4 倍！\n');
