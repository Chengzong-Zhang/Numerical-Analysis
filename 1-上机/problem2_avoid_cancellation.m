%% =========================================================
%  problem2_avoid_cancellation.m
%  功能：设计算法计算 y = x - sin(x)，使有效位丢失最多 1 位
%  =========================================================
%
%  【数值分析背景：相减相消（Cancellation）问题】
%  当 x 接近 0 时，sin(x) ≈ x，两数几乎相等。
%  直接计算 x - sin(x) 会引发"相减相消"（catastrophic cancellation）：
%  两个几乎相等的数相减，虽然绝对误差很小，但结果的相对误差可能极大，
%  导致大量有效位丢失。
%
%  【精度丢失定理（讲义定理 1.3.1）】
%  设 x, y > 0 是规格化二进制浮点数，且 x > y，若
%    2^(-q) ≤ 1 - y/x ≤ 2^(-p)
%  则做减法 x - y 时，至多丢失 q 个、至少丢失 p 个二进制精度位。
%
%  对 f(x) = x - sin(x)（即大数为 x，小数为 sin(x)）：
%    1 - sin(x)/x ≈ x²/6  （x 很小时的 Taylor 近似）
%  因此直接计算约丢失 −log₂(x²/6) = log₂(6/x²) 个二进制精度位。
%
%  【判断"至多丢失 1 位"的阈值】
%  条件：1 - sin(x)/x ≥ 2^(-1) = 0.5
%    => sin(x)/x ≤ 0.5
%    => |x| ≥ x₀，其中 x₀ 约为 1.895（方程 sin(x)/x = 0.5 的根）
%  当 |x| ≥ 1.895 时，直接公式至多丢失 1 个二进制位（精度可接受）。
%  当 |x| < 1.895 时，直接公式丢失精度太多，需用改进算法。
%
%  【改进算法：利用 Taylor 展开避免相减相消】
%  将 sin(x) 的 Taylor 展开代入：
%    sin(x) = x - x³/3! + x⁵/5! - x⁷/7! + ...
%  因此：
%    x - sin(x) = x³/3! - x⁵/5! + x⁷/7! - ...
%               = Σ_{k=1}^{∞} (-1)^{k+1} · x^{2k+1} / (2k+1)!
%  这个级数对所有 x 均收敛，且各项同号不存在相减，无精度丢失问题。
%  对小 x，级数收敛极快（几项即可达到机器精度）。
%  =========================================================

clear; close all; clc;

fprintf('==============================================\n');
fprintf('  题目2：避免相减相消计算 y = x - sin(x)\n');
fprintf('==============================================\n\n');


%% =========================================================
%  第一部分：演示直接公式的相减相消问题
%  =========================================================
%
%  原理：对于小 x，直接公式 x - sin(x) 的结果约为 x³/6，
%  但 x 和 sin(x) 各自的舍入误差约为 x * eps（机器精度量级），
%  当 x 很小时，这个舍入误差可能远大于真实结果 x³/6，导致结果失真。
%  =========================================================

fprintf('--- 第一部分：直接公式的精度损失演示 ---\n\n');

% 选取一系列从大到小的 x 值（涵盖"无问题"到"严重相消"的范围）
x_demo = [1.0, 0.1, 0.01, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8];

fprintf('%-10s  %-22s  %-22s  %-12s  %-8s\n', ...
    'x', '直接公式 x-sin(x)', 'Taylor级数（精确）', '相对误差', '丢失位数');
fprintf('%s\n', repmat('-', 1, 82));

for i = 1:length(x_demo)
    x = x_demo(i);

    % ---- 直接公式（存在相减相消） ----
    y_direct = x - sin(x);

    % ---- Taylor 级数（精确，作为参考真值） ----
    y_taylor = taylor_series(x);

    % ---- 计算相对误差 ----
    % 以 Taylor 级数结果为"真值"，计算直接公式的相对误差
    if y_taylor ~= 0
        rel_err = abs(y_direct - y_taylor) / abs(y_taylor);
    else
        rel_err = 0;
    end

    % ---- 估算丢失的二进制精度位数 ----
    % 由 1 - sin(x)/x ≈ x²/6，丢失位数 ≈ log₂(6/x²)
    if x > 0
        ratio = 1 - sin(x) / x;  % 实际比值（而非近似）
        if ratio > 0
            bits_lost = max(0, -log2(ratio));  % 丢失的二进制位数
        else
            bits_lost = 0;
        end
    else
        bits_lost = 0;
    end

    fprintf('%-10.1e  %-22.15e  %-22.15e  %-12.2e  %.1f\n', ...
        x, y_direct, y_taylor, rel_err, bits_lost);
end

fprintf('\n注：当 x 越小，丢失位数越多，直接公式的相对误差越大。\n');
fprintf('   双精度共 52 个二进制尾数位，丢失 52 位即完全失去精度。\n\n');


%% =========================================================
%  第二部分：改进算法实现
%  =========================================================
%
%  策略：根据 |x| 与阈值的比较，自动选择公式
%    |x| <  THRESHOLD → 使用 Taylor 级数（无相减相消）
%    |x| >= THRESHOLD → 使用直接公式（至多丢失 1 个二进制位）
%
%  阈值 THRESHOLD 的确定：
%    要求直接公式"至多丢失 1 位"，即 1 - sin(x)/x ≥ 1/2
%    数值求解 sin(x)/x = 0.5 得 x ≈ 1.895（在 (0,π) 内唯一）
%  =========================================================

fprintf('--- 第二部分：改进算法——阈值切换 ---\n\n');

% 精确阈值：数值求解 sin(x)/x = 0.5
% 等价于 sin(x) - 0.5*x = 0，在 (1, 2) 内用二分法求解
THRESHOLD = find_threshold();
fprintf('精确阈值（sin(x)/x = 0.5 的解）: x₀ = %.6f\n', THRESHOLD);
fprintf('对比 √3 ≈ %.6f（Taylor 近似 x²/6 = 1/2 时的阈值）\n\n', sqrt(3));

% 验证改进算法的效果：与 Taylor 级数（参考真值）对比
x_test = [1e-10, 1e-8, 1e-5, 1e-3, 0.1, 0.5, 1.0, 1.5, THRESHOLD, 2.0, pi, 10.0];

fprintf('%-12s  %-8s  %-22s  %-22s\n', 'x', '使用公式', '改进算法', '直接公式');
fprintf('%s\n', repmat('-', 1, 70));

for i = 1:length(x_test)
    x = x_test(i);

    if abs(x) < THRESHOLD
        method = 'Taylor';
        y_improved = taylor_series(x);
    else
        method = '直接';
        y_improved = x - sin(x);
    end

    y_direct = x - sin(x);

    fprintf('%-12.4e  %-8s  %-22.15e  %-22.15e\n', ...
        x, method, y_improved, y_direct);
end

fprintf('\n结论：对于小 x，改进算法（Taylor）与大 x 的直接公式结果一致，\n');
fprintf('      但直接公式在小 x 时精度明显较差（最后几位数字不同）。\n\n');


%% =========================================================
%  第三部分：量化验证——比较两种算法的精度
%  =========================================================
%
%  验证方法：
%  由于 Taylor 级数可以计算到机器精度（各项均为有效加法，无精度损失），
%  我们以 Taylor 级数结果为"参考真值"，计算直接公式的相对误差。
%  =========================================================

fprintf('--- 第三部分：精度量化验证 ---\n\n');

% 取 1000 个对数均匀分布的 x 值
% 范围 1e-12 到 3：覆盖"完全相消"区（左端）到"阈值右侧"区（右端）
% 使得阈值 1.895 恰好落在图的右侧，可在双对数图中清晰看到误差曲线与阈值线交叉
N_pts = 1000;
x_range = logspace(-12, 0.5, N_pts);  % logspace(a,b,n)：在 [10^a, 10^b] 内取 n 个对数均匀点

% 预分配结果数组（避免循环内动态扩容，提升性能）
y_direct_all  = zeros(1, N_pts);  % 直接公式结果
y_improved_all = zeros(1, N_pts); % 改进算法结果
rel_err_direct = zeros(1, N_pts); % 直接公式相对误差

for i = 1:N_pts
    x = x_range(i);
    y_ref = taylor_series(x);          % Taylor 级数作为参考真值
    y_direct_all(i)   = x - sin(x);   % 直接公式
    y_improved_all(i) = improved_formula(x, THRESHOLD);  % 改进算法

    if y_ref ~= 0
        rel_err_direct(i) = abs(y_direct_all(i) - y_ref) / abs(y_ref);
    end
end

% 找最大相对误差发生在哪里
[max_err, idx_max] = max(rel_err_direct);
fprintf('直接公式的最大相对误差：%.2e，发生在 x = %.2e\n', max_err, x_range(idx_max));
fprintf('对应 Taylor 真值：%.6e\n', taylor_series(x_range(idx_max)));
fprintf('对应直接公式值：%.6e\n\n', y_direct_all(idx_max));


%% =========================================================
%  第四部分：图形可视化
%  =========================================================

fprintf('--- 第四部分：图形可视化 ---\n\n');

% figure：创建图形窗口；Position：[左边距, 下边距, 宽, 高]（像素）
figure('Name', '题目2：y=x-sin(x) 算法比较', 'NumberTitle', 'off', 'Position', [100 100 1000 450]);

% ---- 左图：两种算法的计算结果（loglog 双对数坐标） ----
subplot(1, 2, 1);

% 直接公式在完全相消区（x 极小时）给出 0，loglog 中 log(0) 无法绘制，
% 将这些点替换为 NaN 使折线在该区间断开，直观反映"结果为零"的失效现象
y_direct_plot = abs(y_direct_all);
y_direct_plot(y_direct_plot == 0) = NaN;  % NaN 点在 loglog 中自动跳过

% loglog：双对数坐标
loglog(x_range, y_direct_plot, 'r-', 'LineWidth', 1.5, 'DisplayName', '直接公式 x-sin(x)');
hold on;
loglog(x_range, abs(y_improved_all), 'b--', 'LineWidth', 1.5, 'DisplayName', '改进算法（阈值切换）');

% 标记阈值位置
% xline：在指定 x 值处画垂直线，在 loglog 轴上同样有效
xline(THRESHOLD, 'k:', 'LineWidth', 1.8, 'DisplayName', sprintf('阈值 x=%.3f', THRESHOLD));

xlabel('x（对数刻度）');
ylabel('|y| = |x - sin(x)|（对数刻度）');
title('两种算法的计算结果');
legend('Location', 'NorthWest');
grid on;

% ---- 右图：直接公式的相对误差（loglog 双对数坐标） ----
%
%  【为什么必须用 loglog，不能用 semilogx（线性 y 轴）？】
%  相对误差的理论近似为 12·eps/x²，跨越从 1（完全相消）到 eps（~1e-16）
%  共约 16 个数量级。线性 y 轴下，阈值 1.895 处误差 ≈ 3·eps ≈ 7e-16，
%  与 y=0 无法区分，阈值线悬在空白区域，与曲线毫无交集。
%  用 loglog 后，误差曲线呈斜率 -2 的直线，阈值线恰好与曲线相交，
%  直观展示"在 x₀ 处误差降到约 1 个二进制位"。
subplot(1, 2, 2);

% 将恰好为 0 的误差替换为 NaN（loglog 不能绘制 0）
rel_err_plot = rel_err_direct;
rel_err_plot(rel_err_plot == 0) = NaN;

% loglog：双对数坐标——误差曲线 ~ 12·eps/x² 呈斜率 -2 的直线
loglog(x_range, rel_err_plot, 'r-', 'LineWidth', 1.5, 'DisplayName', '直接公式相对误差');
hold on;

% 阈值线：xline 在 loglog 轴上画垂直线，此时它与误差曲线在 y ≈ 3·eps 处相交
xline(THRESHOLD, 'k--', 'LineWidth', 1.8, 'DisplayName', sprintf('阈值 x=%.3f', THRESHOLD));

% 1 个二进制位损失对应的误差上界：约 2·eps（即 1 ULP，Unit in the Last Place）
% 2^(-1) × 2 × eps = eps 量级；直接公式在 x = THRESHOLD 处误差 ≈ 3·eps
% 在图上画参考线，直观标定"至多 1 bit 丢失"的误差水平
yline(2*eps, 'g--', 'LineWidth', 1.2, 'DisplayName', '2·eps（≈1 bit 丢失上界）');
yline(eps,   'b:',  'LineWidth', 1.2, 'DisplayName', 'eps（机器精度）');

% 在曲线上标注阈值点：计算阈值处的实际相对误差值，在图上打标记
% （该点是 "至多 1 bit 丢失" 的临界点，直接公式误差 ≈ 3·eps）
err_at_threshold = 12 * eps / THRESHOLD^2;  % 理论近似值
loglog(THRESHOLD, err_at_threshold, 'ko', 'MarkerSize', 10, ...
    'MarkerFaceColor', 'k', 'HandleVisibility', 'off');
% text：在图中指定坐标处添加文字标注
text(THRESHOLD * 1.1, err_at_threshold * 3, ...
    sprintf('(%.3f, %.1f\\cdoteps)', THRESHOLD, err_at_threshold/eps), ...
    'FontSize', 8, 'Color', 'k');

xlabel('x（对数刻度）');
ylabel('直接公式相对误差（对数刻度）');
title({'直接公式精度损失（理论值 ≈ 12·eps/x²）'; ...
       '——误差曲线斜率 -2，阈值处交于 ~3·eps'});
legend('Location', 'SouthWest');
grid on;

% 整体标题
sgtitle('避免相减相消：改进算法计算 y = x - sin(x)');
% sgtitle：为整个 figure 添加总标题

fprintf('图形已生成，请查看弹出的图形窗口。\n');
fprintf('\n程序运行完毕。\n');


%% =========================================================
%  局部函数（Local Functions）—— 必须置于脚本末尾
%  =========================================================

function y = taylor_series(x)
    % taylor_series: 用 Taylor 级数计算 x - sin(x)，无相减相消
    %
    % 级数公式：
    %   x - sin(x) = x³/3! - x⁵/5! + x⁷/7! - ...
    %              = Σ_{k=1}^{∞} (-1)^{k+1} · x^{2k+1} / (2k+1)!
    %
    % 递推关系（相邻两项的比值）：
    %   term_{k+1} = term_k · (-x²) / ((2k+2)(2k+3))
    %   分析：x^{2k+3}/(2k+3)! ÷ x^{2k+1}/(2k+1)! = x²/((2k+2)(2k+3))
    %         符号交替故乘以 -1
    %
    % 收敛判据：当下一项的绝对值 < 当前和的绝对值 × eps 时停止

    if x == 0
        y = 0;
        return;
    end

    % 第一项：x³/3! = x³/6（k=1 时的项）
    term = x^3 / 6;
    y = term;
    k = 1;  % k：当前项的序号

    % 循环累加级数，直到收敛或超过最大迭代次数
    while abs(term) > abs(y) * eps
        % 由 term_k 递推 term_{k+1}：乘以 -x² / ((2k+2)(2k+3))
        % 当 k=1 时：分母为 4×5=20，即 -x⁵/120 = -x⁵/5!（正确）
        term = term * (-x^2) / ((2*k + 2) * (2*k + 3));
        y = y + term;
        k = k + 1;

        % 安全出口：防止对极大 x 无限循环（实际此函数仅用于小 x）
        if k > 200
            break;
        end
    end
end


function y = improved_formula(x, threshold)
    % improved_formula: 改进算法，根据 |x| 与阈值自动切换公式
    %
    % 输入：
    %   x         - 计算点（标量）
    %   threshold - 切换阈值（|x| < threshold 时用 Taylor 级数）
    % 输出：
    %   y         - x - sin(x) 的精确近似值

    if abs(x) < threshold
        % 小 x：用 Taylor 级数，彻底消除相减相消
        y = taylor_series(x);
    else
        % 大 x：直接公式（sin(x)/x ≤ 0.5，至多丢失 1 个二进制位）
        y = x - sin(x);
    end
end


function x0 = find_threshold()
    % find_threshold: 二分法数值求解 sin(x)/x = 0.5，即 sin(x) - 0.5x = 0
    %
    % 函数 g(x) = sin(x) - 0.5x 的性质：
    %   g(1) = sin(1) - 0.5 ≈ 0.341 > 0
    %   g(2) = sin(2) - 1.0 ≈ -0.091 < 0
    %   由零点定理，根在 (1, 2) 内
    %
    % 二分法：每次将区间折半，取中点判断符号，缩小含根区间

    a = 1.0;   % 区间左端点（g(a) > 0）
    b = 2.0;   % 区间右端点（g(b) < 0）

    % g(x) = sin(x) - 0.5*x
    % 迭代直到区间长度小于机器精度
    while (b - a) > eps * abs(a)
        mid = (a + b) / 2;             % 区间中点
        g_mid = sin(mid) - 0.5 * mid;  % 在中点处计算 g

        if g_mid > 0
            a = mid;  % 根在 [mid, b] 内
        else
            b = mid;  % 根在 [a, mid] 内
        end
    end

    x0 = (a + b) / 2;  % 返回最终近似根
end
