%% =========================================================
%  模板T08：自适应 Simpson 积分法
%
%  【核心思想】
%  不用均匀步长，而是根据函数变化自动决定在哪里细分：
%  - 函数变化剧烈的地方：自动加密节点（小步长）
%  - 函数变化平缓的地方：保持大步长（节省计算）
%  在保证精度的前提下，使函数求值次数最少。
%
%  【误差估计】
%  对子区间 [a,b]，设中点 m=(a+b)/2：
%    整区间 Simpson：    S(a,b)
%    对分后两个Simpson之和：S(a,m) + S(m,b)
%
%  由 Simpson 误差阶数 O(h⁵)，步长减半则误差缩小 2⁵=32 倍：
%    S(a,m)+S(m,b) 的误差 ≈ [S(a,b) 的误差] / 16
%  因此局部误差估计为：
%    |真实误差| ≈ |S(a,m)+S(m,b) - S(a,b)| / 15   （15 = 16-1）
%
%  【接受条件】
%  若 |S(a,m)+S(m,b) - S(a,b)| / 15 ≤ tol，接受结果：
%    I(a,b) ≈ S(a,m)+S(m,b) + [S(a,m)+S(m,b)-S(a,b)]/15
%                              ↑这是 Richardson 校正项，更精确
%  否则：递归对 [a,m] 和 [m,b] 各用 tol/2 的容限继续细化
%  =========================================================

clear; close all; clc;

%% ===== ① 修改这里：定义你的积分 =====
f = @(x) 4 ./ (1 + x.^2);   % 被积函数（逐元素运算！）
a = 0;                        % 积分下限
b = 1;                        % 积分上限
tol = 1e-10;                  % 整体误差容限（允许的最大绝对误差）
I_exact = pi;                 % 精确值（用于验证）
%% =====================================

fprintf('自适应 Simpson 积分法\n');
fprintf('积分区间：[%g, %g]，误差容限：%.0e，精确值：%.15f\n\n', a, b, tol, I_exact);

%% ---- 调用自适应积分 ----
% 先计算初始的整区间 Simpson 值，作为递归的起点
S_init = simp1(f, a, b);
% 然后调用递归函数（定义在文件末尾）
[result, func_count] = adapt_simp(f, a, b, tol, S_init);

fprintf('自适应 Simpson 结果：\n');
fprintf('  近似值   = %.15f\n', result);
fprintf('  精确值   = %.15f\n', I_exact);
fprintf('  误  差   = %.4e\n', abs(result - I_exact));
fprintf('  函数求值 = %d 次（+ 初始3次 = %d 次总计）\n\n', func_count, func_count+3);

%% ---- 与固定步长方法的函数求值次数对比 ----
fprintf('--- 达到同等精度所需的函数求值次数对比 ---\n');

% 梯形公式：找到达到相同精度所需的 n
err_target = abs(result - I_exact);
for n_t = [100, 1000, 10000, 100000]
    h_t = (b-a)/n_t;
    x_t = linspace(a, b, n_t+1);
    fx_t = f(x_t);
    T_t = h_t/2*(fx_t(1)+fx_t(end)+2*sum(fx_t(2:end-1)));
    if abs(T_t - I_exact) <= err_target * 10   % 接近目标精度
        fprintf('  复合梯形达到类似精度需要 n=%d（%d次函数求值）\n', n_t, n_t+1);
        break;
    end
end
fprintf('  自适应Simpson只需 %d 次函数求值\n', func_count+3);
fprintf('  自适应方法更高效！（在函数变化不均匀时优势更明显）\n\n');

%% ---- 测试不同容限 ----
fprintf('--- 不同精度要求下的性能 ---\n');
fprintf('%-12s  %-20s  %-15s  %-12s\n', '误差容限', '计算结果', '实际误差', '函数求值次数');
fprintf('%s\n', repmat('-', 1, 62));
for tol_k = [1e-4, 1e-6, 1e-8, 1e-10, 1e-12]
    S0 = simp1(f, a, b);
    [res_k, cnt_k] = adapt_simp(f, a, b, tol_k, S0);
    fprintf('%-12.0e  %-20.15f  %-15.4e  %-12d\n', ...
        tol_k, res_k, abs(res_k-I_exact), cnt_k+3);
end
fprintf('\n');

%% =========================================================
%  局部函数（必须放在脚本文件末尾）
%  =========================================================

function S = simp1(f, a, b)
    % 单区间 Simpson 公式：S = (b-a)/6 · [f(a) + 4f(m) + f(b)]
    % m = (a+b)/2 是区间中点
    % 这里 (b-a)/6 = h/3，其中 h=(b-a)/2（单区间半长）
    m = (a + b) / 2;                   % 区间中点
    S = (b - a) / 6 * (f(a) + 4*f(m) + f(b));
    % 代数精度：3（能精确积分 1, x, x², x³ 这四类函数）
end

function [I, cnt] = adapt_simp(f, a, b, tol, S_ab)
    % 自适应 Simpson 积分（递归实现）
    %
    % 输入：
    %   f     - 被积函数（匿名函数）
    %   a, b  - 当前子区间的两端点
    %   tol   - 本子区间允许的最大误差（整体误差按区间自动分配）
    %   S_ab  - 已计算好的整区间 [a,b] 的 Simpson 值（避免重复计算）
    %
    % 输出：
    %   I   - 本子区间的积分近似值
    %   cnt - 本次递归调用新增的函数求值次数
    %
    % 算法流程：
    % 1. 将 [a,b] 从中点 m 对分
    % 2. 分别计算 S(a,m) 和 S(m,b)
    % 3. 估计误差：|S(a,m)+S(m,b)-S(a,b)| / 15
    % 4. 若误差 ≤ tol：接受（加上校正项）
    %    否则：递归细化 [a,m] 和 [m,b]

    m = (a + b) / 2;   % 将区间从中点对分

    % 计算对分后两个子区间各自的 Simpson 值
    S_am = simp1(f, a, m);   % [a,m] 上的 Simpson
    S_mb = simp1(f, m, b);   % [m,b] 上的 Simpson

    S2 = S_am + S_mb;   % 对分后两个 Simpson 之和

    % 新增函数求值次数：
    % S_am 用到 f(a), f((a+m)/2), f(m)
    % S_mb 用到 f(m), f((m+b)/2), f(b)
    % 其中 f(a), f(b) 上层已知（来自 S_ab），f(m) 算了一次（两边共用）
    % 新增：f((a+m)/2), f(m), f((m+b)/2) 共 3 个
    cnt = 3;

    % 局部误差估计
    err_est = abs(S2 - S_ab) / 15;
    % 15 = 2⁴ - 1：由 Simpson 误差 O(h⁵)，步长减半使误差缩小 32 倍，
    % 故 S2 的误差 ≈ S_ab 误差 / 16，Richardson 外推估计误差 = (S2-S_ab)/15

    if err_est <= tol
        % 误差满足容限：接受结果，并加上 Richardson 校正项
        % 校正 (S2-S_ab)/15 是高阶误差的修正，使精度进一步提升
        I = S2 + (S2 - S_ab) / 15;
    else
        % 误差不满足容限：递归细化两个子区间
        % 误差容限平均分配给两个子区间（各给 tol/2）
        [I_left,  cnt_left]  = adapt_simp(f, a, m, tol/2, S_am);
        [I_right, cnt_right] = adapt_simp(f, m, b, tol/2, S_mb);
        I   = I_left + I_right;
        cnt = cnt + cnt_left + cnt_right;   % 累计函数求值次数
    end
end
