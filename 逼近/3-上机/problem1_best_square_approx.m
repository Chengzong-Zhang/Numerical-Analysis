%% =========================================================
%  problem1_best_square_approx.m
%  功能：求 f(x) = e^x 在区间 [0,1] 上的 n 次最佳平方逼近多项式
%        φ(x) = Σ_{k=0}^{n} a_k φ_k(x)
%
%  题目四种情形：
%    (1) n=5,  φ_k(x) = x^k（幂函数基）
%    (2) n=5,  φ_k(x) = k 次 Legendre 多项式
%    (3) n=10, φ_k(x) = x^k（幂函数基）
%    (4) n=10, φ_k(x) = k 次 Legendre 多项式
%
%  【数值分析背景：最佳平方逼近（教材§3.4）】
%  在带权 W(x)=1 的内积空间 C[0,1] 中，内积定义为：
%    (f, g) := ∫_0^1 f(x)g(x)dx
%  最佳平方逼近 φ*(x) ∈ Φ_{n+1} = span{φ_0,...,φ_n} 满足（教材式3.4.1）：
%    ||f - φ*||_2^2 = min_{φ∈Φ_{n+1}} ∫_0^1 [f(x)-φ(x)]^2 dx
%
%  由极值条件（对各系数 c_j 求偏导令其为0，教材式3.4.4→3.4.5），得法方程组：
%    G c = f_vec
%  其中 Gram 矩阵 G_{ij} = (φ_i, φ_j) = ∫_0^1 φ_i(x)φ_j(x)dx
%  右端向量 f_vec_j = (f, φ_j) = ∫_0^1 e^x φ_j(x)dx
%
%  两种基函数的关键对比（见教材§3.4理论）：
%  (A) 幂函数基 φ_k(x)=x^k：G 为 Hilbert 矩阵（著名病态矩阵），
%      条件数随 n 指数增长，n=10 时条件数约 10^13，数值解极不稳定
%  (B) Legendre 正交基 φ_k(x)=P_k(2x-1)：G 为对角阵，条件数=1，
%      可直接逐项计算系数，无需解方程组
%  =========================================================

clear; close all; clc;

fprintf('==============================================\n');
fprintf('  第三章上机作业1：f(x)=e^x 的最佳平方逼近\n');
fprintf('  φ(x) = Σ_{k=0}^n a_k φ_k(x)\n');
fprintf('==============================================\n\n');

% --- 定义被逼近函数 ---
% @(x) exp(x)：匿名函数（Anonymous Function）语法
% 语法：f = @(参数列表) 表达式
% 特点：支持向量化——当 x 是向量时，exp(x) 对每个分量分别计算
% 数值分析含义：f(x) = e^x 是光滑函数，可用多项式在 L2 意义下任意精度逼近
f = @(x) exp(x);

% --- ||f||_2^2 的解析值（用于误差公式） ---
% ||f||_2^2 = ∫_0^1 e^{2x} dx = [e^{2x}/2]_0^1 = (e^2-1)/2
% 数值分析含义：这是误差公式 ||f-φ*||_2^2 = ||f||_2^2 - Σ c_k^2/(2k+1) 中的第一项
f_L2sq = (exp(2) - 1) / 2;   % 解析值，无数值误差

% --- 绘图用的密集点（不参与计算，只用于可视化） ---
% linspace(a, b, N)：在 [a,b] 内均匀生成 N 个点（含两端点）
% 语法：x = linspace(起点, 终点, 点数)
x_plot = linspace(0, 1, 500);  % 500 个点足够光滑
y_exact = f(x_plot);            % 精确值向量，用于对比绘图

% 存储各情形结果（用于最后统一绘图）
results = struct();
idx_result = 0;

%% =========================================================
%  情形(1)(3)：幂函数基 φ_k(x) = x^k
%  Gram 矩阵 = Hilbert 矩阵（高度病态）
%  =========================================================

fprintf('===== 幂函数基 φ_k(x) = x^k =====\n');
fprintf('  理论背景：\n');
fprintf('  G_{ij} = (x^i, x^j) = ∫_0^1 x^{i+j}dx = 1/(i+j+1)\n');
fprintf('  → G 是 Hilbert 矩阵，条件数 κ(G) 随 n 指数增长\n');
fprintf('  → 高 n 时法方程组数值求解极不稳定（损失有效数字）\n\n');

for n = [5, 10]
    idx_result = idx_result + 1;
    label_str = sprintf('(%d) n=%d, 幂函数基', 2*idx_result-1, n);
    fprintf('--- %s ---\n', label_str);

    % ---- 构造 Gram 矩阵（Hilbert 矩阵） ----
    % hilb(m)：MATLAB 内置函数，生成 m×m Hilbert 矩阵
    % 语法：H = hilb(m)
    % 定义：H(i,j) = 1/(i+j-1)，i,j = 1,...,m
    % 与 Gram 矩阵的关系：G(i+1,j+1) = (x^i, x^j) = 1/(i+j+1) = H(i+1,j+1)
    % 因此 hilb(n+1) 恰好给出 (n+1)×(n+1) 的 Gram 矩阵
    G = hilb(n+1);
    % 注意：hilb(n+1) 的 (i,j) 元素 = 1/(i+j-1)，
    % 对应 G_{i-1,j-1} = ∫_0^1 x^{i-1} x^{j-1} dx = 1/(i+j-1) ✓

    % ---- 构造右端向量 f_vec ----
    % f_vec(j+1) = (f, φ_j) = ∫_0^1 e^x · x^j dx，j=0,1,...,n
    % integral(fun, a, b)：自适应数值积分（Gauss-Kronrod 求积法）
    % 语法：Q = integral(@(x) 被积函数, 下限, 上限)
    % 数值分析含义：计算 f 与各基函数的内积，构成法方程组右端项
    % 注意：x.^j 表示向量 x 的各元素取 j 次方（逐元素运算符 .^）
    f_vec = zeros(n+1, 1);    % zeros(m,n)：生成 m×n 全零矩阵，预分配内存
    for j = 0:n
        f_vec(j+1) = integral(@(x) exp(x) .* x.^j, 0, 1);
        % .* 是逐元素乘法：当 x 为向量时，对应位置相乘
    end

    % ---- 条件数分析 ----
    % cond(A)：矩阵 A 的 2-范数条件数 κ(A) = σ_max(A)/σ_min(A)
    % 语法：k = cond(A)
    % 数值分析含义：求解 Gc=f_vec 时，系数的相对误差 ≤ κ(G) × 机器精度(eps)
    % 即：若 κ(G) ≈ 10^k，则求解过程约损失 k 位有效数字
    kappa_G = cond(G);
    fprintf('  Gram 矩阵（%d×%d Hilbert 矩阵）条件数：\n', n+1, n+1);
    fprintf('  cond(G) = %.4e\n', kappa_G);
    fprintf('  （机器精度 eps≈2.2e-16，损失约 %.0f 位有效数字）\n', log10(kappa_G));
    % log10(x)：以 10 为底的对数，用来估计损失的十进制有效位数

    % ---- 求解法方程组 G·c = f_vec ----
    % A \ b（左除，反斜杠算子）：求解线性方程组 Ax=b
    % 语法：x = A \ b
    % 数值分析含义：等价于 x = inv(A)*b，但比先算 inv(A) 再相乘更数值稳定
    % MATLAB 内部对正定对称矩阵使用 Cholesky 分解
    % 当 n=10 时，G 极度病态，c 的精度严重下降（这正是问题的核心）
    c_mono = G \ f_vec;
    % c_mono(k+1) = a_k，即 φ*(x) 中 x^k 的系数

    % ---- 在绘图点上求值 φ*(x) = Σ_{k=0}^n c_k · x^k ----
    % 使用 Horner 嵌套乘法（秦九韶算法），减少浮点运算次数
    % Horner 变形：φ*(x) = (...((c_n · x + c_{n-1})·x + c_{n-2})·x + ... + c_0)
    phi_mono = eval_monomial(c_mono, x_plot);
    % eval_monomial：本文件末定义的局部函数

    % ---- 计算 L2 误差和最大误差 ----
    % L2 误差：||f-φ*||_2 = sqrt(∫_0^1 [e^x - φ*(x)]^2 dx)
    % 用 integral 数值积分计算
    L2_err = sqrt(integral(@(x) (exp(x) - eval_monomial(c_mono, x)).^2, 0, 1));
    % max(v)：向量 v 的最大元素；abs(v)：向量逐元素取绝对值
    max_err = max(abs(y_exact - phi_mono));

    fprintf('  L2 误差  ||f-φ*||_2 = %.4e\n', L2_err);
    fprintf('  最大误差 ||f-φ*||_∞ = %.4e\n\n', max_err);

    % 存储结果
    results(idx_result).label   = label_str;
    results(idx_result).phi     = phi_mono;
    results(idx_result).L2_err  = L2_err;
    results(idx_result).max_err = max_err;
    results(idx_result).kappa   = kappa_G;
end

%% =========================================================
%  情形(2)(4)：k 次 Legendre 多项式 φ_k(x) = P_k(2x-1)
%  Gram 矩阵为对角阵（正交基），系数可逐项显式计算
%  =========================================================

fprintf('===== Legendre 正交基 φ_k(x) = P_k(2x-1) =====\n');
fprintf('  理论背景（教材§3.3 和§3.4）：\n');
fprintf('  Legendre 多项式 P_k(t) 在 [-1,1] 上关于权函数 W=1 正交：\n');
fprintf('    ∫_{-1}^1 P_m(t)P_k(t)dt = 2/(2k+1)·δ_{mk}\n');
fprintf('  变量代换 t=2x-1（把[0,1]映射到[-1,1]），令 φ_k(x)=P_k(2x-1)：\n');
fprintf('    ∫_0^1 φ_i(x)φ_j(x)dx = δ_{ij}/(2i+1)  ← Gram矩阵是对角阵!\n');
fprintf('  因此系数公式（教材式3.4.7）直接给出：\n');
fprintf('    c_k = (f,φ_k)/(φ_k,φ_k) = (2k+1)·∫_0^1 e^x·P_k(2x-1)dx\n\n');

for n = [5, 10]
    idx_result = idx_result + 1;
    label_str = sprintf('(%d) n=%d, Legendre 基', 2*idx_result-3, n);
    fprintf('--- %s ---\n', label_str);

    % ---- 逐项计算各 Legendre 展开系数 ----
    % 由正交性，Gram 矩阵对角元 (φ_k,φ_k) = 1/(2k+1)（对角阵，条件数=1）
    % 故 c_k = (f,φ_k)/(φ_k,φ_k) = (2k+1)·∫_0^1 e^x·P_k(2x-1)dx
    c_leg = zeros(n+1, 1);     % 预分配系数向量
    for k = 0:n
        % integral：自适应数值积分
        % legendre_on_01(k, x)：本文件末局部函数，计算 P_k(2x-1) 的向量值
        inner_fk = integral(@(x) exp(x) .* legendre_on_01(k, x), 0, 1);
        % 内积 (f, φ_k) = ∫_0^1 e^x · P_k(2x-1) dx
        norm_sq_k = 1 / (2*k + 1);    % (φ_k, φ_k) = 1/(2k+1)
        c_leg(k+1) = inner_fk / norm_sq_k;   % c_k = 内积 / 模长平方
        % 等价写法：c_leg(k+1) = (2*k+1) * inner_fk;
    end

    fprintf('  Gram 矩阵（对角阵）条件数 = 1（正交基的最优性）\n');

    % ---- 在绘图点上求值 φ*(x) = Σ_{k=0}^n c_k·P_k(2x-1) ----
    phi_leg = zeros(size(x_plot));   % 预分配输出，初始化为零
    for k = 0:n
        % legendre_on_01(k, x_plot)：计算 P_k(2x-1) 在密集绘图点上的值
        phi_leg = phi_leg + c_leg(k+1) * legendre_on_01(k, x_plot);
        % + 是逐元素加法（向量加法），* 是标量乘以向量
    end

    % ---- 利用公式(3.4.9)计算 L2 误差 ----
    % ||f-φ*||_2^2 = ||f||_2^2 - Σ_{k=0}^n c_k^2·(φ_k,φ_k)
    %              = (e^2-1)/2 - Σ_{k=0}^n c_k^2/(2k+1)
    % 数值分析含义：利用 Bessel 不等式，L2 误差随 n 增大单调减小
    norms_sq = 1 ./ (2*(0:n)' + 1);    % 列向量 [1, 1/3, 1/5, ..., 1/(2n+1)]
    % (0:n)'：行向量 [0,1,...,n] 的转置，变成列向量
    % ./ 是逐元素除法
    L2_err_sq = f_L2sq - sum(c_leg.^2 .* norms_sq);
    % c_leg.^2：系数向量逐元素平方；.* norms_sq：逐元素乘以各模长平方
    % sum(v)：向量 v 所有元素之和
    % 理论上 L2_err_sq ≥ 0；若因数值误差略为负数，取 max(0,...)
    L2_err_sq = max(0, L2_err_sq);
    L2_err = sqrt(L2_err_sq);
    max_err = max(abs(y_exact - phi_leg));

    fprintf('  L2 误差  ||f-φ*||_2 = %.4e\n', L2_err);
    fprintf('  最大误差 ||f-φ*||_∞ = %.4e\n\n', max_err);

    results(idx_result).label   = label_str;
    results(idx_result).phi     = phi_leg;
    results(idx_result).L2_err  = L2_err;
    results(idx_result).max_err = max_err;
    results(idx_result).kappa   = 1;    % 对角 Gram 矩阵，条件数 = 1
end

%% =========================================================
%  汇总对比表
%  =========================================================

fprintf('==============================================\n');
fprintf('  四种情形误差与条件数汇总\n');
fprintf('==============================================\n');
fprintf('%-28s  %-12s  %-12s  %-12s\n', '情形', 'L2 误差', '最大误差', 'cond(G)');
fprintf('%s\n', repmat('-', 1, 72));
% repmat('-', 1, 72)：将字符 '-' 重复 72 次，生成分隔线字符串
for k = 1:4
    fprintf('%-28s  %-12.4e  %-12.4e  %-12.4e\n', ...
        results(k).label, results(k).L2_err, results(k).max_err, results(k).kappa);
    % %-28s：左对齐字符串，宽度28；%-12.4e：左对齐科学计数法，宽度12，4位小数
end
fprintf('\n');
fprintf('观察：\n');
fprintf('  1. 幂函数基（Hilbert矩阵）与Legendre基逼近效果相近（理论上等价）\n');
fprintf('  2. Legendre基条件数=1，幂函数基条件数极大，数值稳定性天壤之别\n');
fprintf('  3. n越大，幂函数基条件数越大，Legendre基始终保持完美数值稳定\n\n');

%% =========================================================
%  绘图：四种情形的逼近效果
%  =========================================================

% figure：创建新图形窗口
% 'Name'：窗口标题；'NumberTitle','off'：不显示默认"Figure 1"等编号
% 'Position'：[左边距 下边距 宽度 高度]（单位：像素）
figure('Name', 'f(x)=e^x 的最佳平方逼近', 'NumberTitle', 'off', ...
    'Position', [50, 50, 1200, 900]);

% 四个子图的标题字符串
% 子图标题与 results(1..4) 一一对应
subplot_titles = {
    '(1) n=5,  幂函数基 x^k（Hilbert矩阵）', ...   % results(1)
    '(3) n=10, 幂函数基 x^k（极度病态）', ...        % results(2)
    '(2) n=5,  Legendre 正交基', ...                % results(3)
    '(4) n=10, Legendre 正交基'                     % results(4)
};
% 绘图顺序：上行为幂函数基（n=5,n=10），下行为Legendre基（n=5,n=10）
% results(1)=幂n=5, results(2)=幂n=10, results(3)=Legn=5, results(4)=Legn=10
plot_order = [1, 2, 3, 4];   % results 的下标顺序（按存储顺序直接对应子图）

for k = 1:4
    r = plot_order(k);        % 取对应结果
    % subplot(m, n, k)：将图窗划分为 m 行 × n 列，当前绘第 k 个子图
    subplot(2, 2, k);

    % plot(x, y, 线型颜色, 属性名, 值, ...)：绘制二维曲线
    % 'k-'：黑色(k)实线(-); 'LineWidth'：线宽; 'DisplayName'：图例标签
    plot(x_plot, y_exact, 'k-', 'LineWidth', 2, 'DisplayName', 'f(x)=e^x 精确值');
    % hold on：保持当前坐标轴，后续绘图叠加在上面
    hold on;
    % 'b--'：蓝色(b)虚线(--)
    plot(x_plot, results(r).phi, 'b--', 'LineWidth', 1.5, 'DisplayName', '逼近多项式 φ*(x)');

    % xlabel/ylabel：设置坐标轴标签
    xlabel('x', 'FontSize', 10);
    ylabel('y', 'FontSize', 10);
    % title：设置子图标题，cell 数组 {} 中多个字符串各占一行
    % sprintf：格式化字符串（类似C语言的printf，返回字符串而不输出）
    title({subplot_titles{k}, ...
           sprintf('L2误差=%.2e, cond(G)=%.2e', results(r).L2_err, results(r).kappa)}, ...
          'FontSize', 9);
    % legend：添加图例；'Location','best'：自动选择最佳位置
    legend('Location', 'northwest', 'FontSize', 8);
    % grid on：显示网格线
    grid on;
    hold off;   % 释放 hold，下一个 subplot 重新开始
end

% sgtitle：为整个 figure 添加总标题（Subplot Group Title，R2018b+）
sgtitle('f(x) = e^x 在 [0,1] 上的 n 次最佳平方逼近多项式比较', 'FontSize', 13);

fprintf('图形已生成，请查看弹出的图形窗口。\n\n');
fprintf('程序运行完毕。\n');


%% =========================================================
%  局部函数（Local Functions）—— 必须置于脚本末尾
%  =========================================================

function vals = eval_monomial(c, x)
    % eval_monomial：用 Horner 嵌套法（秦九韶算法）计算幂函数基多项式的值
    %
    % 【数值分析：Horner 嵌套乘法（教材第二章）】
    % 多项式 φ*(x) = c_0 + c_1·x + c_2·x^2 + ... + c_n·x^n
    % 直接计算需 n 次乘方（x^k）+ n 次乘法 + n 次加法
    % Horner 变形：φ*(x) = (...((c_n·x + c_{n-1})·x + c_{n-2})·x + ... + c_0)
    % 只需 n 次乘法 + n 次加法，计算量减半，且数值误差更小
    %
    % 输入：
    %   c - 系数向量（长度 n+1），c(k+1) 是 x^k 的系数 a_k
    %   x - 求值点（可以是向量）
    % 输出：
    %   vals - φ*(x) 在各点的值（与 x 形状相同）

    n = length(c) - 1;               % 多项式次数
    vals = c(n+1) * ones(size(x));   % 从最高次系数开始初始化
    % ones(size(x))：生成与 x 形状相同的全 1 数组，用于标量×向量

    % Horner 递推：从 n-1 次开始向下递推到 0 次
    for k = n : -1 : 1
        % n:-1:1 表示从 n 递减到 1 的整数序列（步长 -1）
        vals = vals .* x + c(k);
        % .* 是逐元素乘法（x 为向量时各元素分别乘），确保向量化
        % + c(k) 是标量加到向量的每个元素上（MATLAB 自动广播）
    end
end


function vals = legendre_on_01(k, x)
    % legendre_on_01：计算第 k 次 Legendre 多项式 P_k 在 t=2x-1 处的值
    %
    % 【数值分析：Legendre 多项式（教材§3.3 (二)）】
    % Legendre 多项式 P_k(t) 定义在 [-1,1] 上，满足（教材式3.3.17）：
    %   P_k(t) = 1/(2^k k!) · d^k/dt^k (t^2-1)^k
    %
    % 三项递推关系（教材式3.3.21，避免逐阶微分，数值效率高）：
    %   P_0(t) = 1,  P_1(t) = t
    %   P_{k+1}(t) = (2k+1)/(k+1) · t · P_k(t) - k/(k+1) · P_{k-1}(t)
    %
    % 令 t = 2x-1（将 x∈[0,1] 映射到 t∈[-1,1]），则 φ_k(x) = P_k(2x-1)
    % 正交性验证：∫_0^1 P_i(2x-1)P_j(2x-1)dx = δ_{ij}/(2i+1)（教材导出）
    %
    % 输入：
    %   k - 多项式阶数（非负整数）
    %   x - 求值点（在 [0,1] 内，可以是向量）
    % 输出：
    %   vals - P_k(2x-1) 在各点的值（与 x 形状相同）

    % 变量代换：[0,1] → [-1,1]，t = 2x-1
    t = 2*x - 1;
    % 逐元素运算：当 x 为向量时，t 也是同形状的向量

    if k == 0
        % P_0(t) = 1（常数多项式）
        vals = ones(size(x));
        % ones(size(x))：生成与 x 等大小的全 1 数组

    elseif k == 1
        % P_1(t) = t = 2x-1（一次多项式）
        vals = t;

    else
        % 三项递推（从 P_0, P_1 出发，推到 P_k）
        p_prev = ones(size(x));   % P_0(t) = 1
        p_curr = t;                % P_1(t) = t
        for j = 1 : k-1
            % 第 j+1 步递推，从 P_j 计算 P_{j+1}
            % 教材式(3.3.21)：P_{j+1}(t) = (2j+1)/(j+1)·t·P_j(t) - j/(j+1)·P_{j-1}(t)
            p_next = ((2*j+1) * t .* p_curr - j * p_prev) / (j+1);
            % .* 是逐元素乘法（t 和 p_curr 都是向量）
            % p_next 与 x 形状相同
            p_prev = p_curr;   % 更新"前一项"
            p_curr = p_next;   % 更新"当前项"
        end
        vals = p_curr;   % 循环结束后 p_curr = P_k(t)
    end
end
