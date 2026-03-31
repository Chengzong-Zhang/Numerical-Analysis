%% =========================================================
%  problem1_interpolation.m
%  功能：对 R(x) = 1/(1+x^2), x∈[-5,5] 用五种插值方法逼近，并绘图比较
%  =========================================================
%
%  【数值分析背景知识：插值问题】
%  插值问题：已知函数 f(x) 在 n+1 个节点 x_0, x_1, ..., x_n 处的值，
%  构造一个"简单函数" P(x) 使得 P(x_k) = f(x_k)，k=0,1,...,n。
%  目标：用 P(x) 在节点以外的点逼近 f(x)。
%
%  本题的被插值函数：Runge 函数 R(x) = 1/(1+x²)
%  这个函数以"Runge 现象"著称：在等距节点上，高次多项式插值在端点附近
%  会出现剧烈振荡，误差不随次数增大而减小，甚至发散。
%
%  五种方法：
%  (1) 等距节点 Newton 插值（10次）→ 展示 Runge 现象
%  (2) Chebyshev 节点 Lagrange 插值（20次）→ 压制 Runge 现象
%  (3) 等距节点分段线性插值 → 简单可靠，但光滑性差
%  (4) 等距节点三次自然样条插值 → 光滑（二阶连续）
%  (5) 等距节点分段三次 Hermite 插值 → 保形，单调性好
%  =========================================================

% clear：清除工作区所有变量，避免前次运行残留数据干扰
% close all：关闭所有图形窗口
% clc：清空命令行显示内容
clear; close all; clc;

fprintf('==============================================\n');
fprintf('  第二章上机作业：Runge 函数的五种插值方法\n');
fprintf('==============================================\n\n');

%% =========================================================
%  公共准备：定义被插值函数与绘图密度
%  =========================================================

% 定义匿名函数（Anonymous Function）
% 语法：f = @(x) 表达式
% 作用：将数学公式封装为可调用的函数对象，支持向量化（x 可以是向量）
% 数值分析含义：R(x) = 1/(1+x²) 是 Runge 函数，在 x=±5 附近变化平缓，
%               但等距高次插值时在端点附近会产生大振荡（Runge 现象）
R = @(x) 1 ./ (1 + x.^2);
% ./ 是逐元素除法，当 x 为向量时，对每个元素分别计算 1/(1+x_k²)

% 生成密集的绘图点，用于绘制光滑曲线
% linspace(a, b, n)：在 [a,b] 内均匀生成 n 个点（含端点）
% 数值分析含义：这些点只用于绘图评估，不是插值节点
x_plot = linspace(-5, 5, 500);  % 500个点已足够光滑

% 计算 R(x) 在绘图点上的精确值（作为参考曲线）
y_exact = R(x_plot);

fprintf('被插值函数：R(x) = 1/(1+x²)，区间 [-5, 5]\n\n');


%% =========================================================
%  第(1)部分：等距节点 Newton 插值多项式（10次）
%  =========================================================
%
%  【数值分析：Newton 插值与差商表】
%  Newton 插值多项式利用"均差"（divided difference）构造，形式为：
%    N_n(x) = f[x_0]
%           + f[x_0,x_1](x-x_0)
%           + f[x_0,x_1,x_2](x-x_0)(x-x_1)
%           + ...
%           + f[x_0,...,x_n](x-x_0)...(x-x_{n-1})
%
%  均差（差商）的递推定义：
%    f[x_k] = f(x_k)                                （零阶）
%    f[x_k,x_{k+1}] = (f[x_{k+1}]-f[x_k])/(x_{k+1}-x_k)  （一阶）
%    f[x_k,...,x_{k+j}] = (f[x_{k+1},...,x_{k+j}] - f[x_k,...,x_{k+j-1}])
%                         / (x_{k+j} - x_k)          （高阶）
%
%  Newton 插值的优点：增加节点时只需计算新的差商，无需重算所有系数。
%
%  【等距节点的 Runge 现象】
%  对 R(x) 在 [-5,5] 取 11 个等距节点做 10 次多项式插值，
%  插值多项式在端点 x=±5 附近会产生剧烈振荡，误差极大。
%  这是因为等距节点的 Lagrange 基函数在端点附近值很大，
%  而高次多项式在端点附近"不受控制"。
%  =========================================================

fprintf('--- 第(1)部分：等距节点 Newton 插值（10次）---\n\n');

% 定义等距插值节点：x_k = -5+k, k=0,1,...,10
% 共 11 个节点，插值多项式次数为 10
n1 = 10;                          % 插值多项式的次数
% linspace(-5, 5, n1+1)：在 [-5,5] 均匀取 n1+1 = 11 个点
x_nodes1 = linspace(-5, 5, n1+1); % 行向量，11个等距节点
y_nodes1  = R(x_nodes1);           % 各节点处的函数值

fprintf('等距插值节点：x_k = -5 + k, k = 0,1,...,10\n');
fprintf('节点数：%d，多项式次数：%d\n\n', length(x_nodes1), n1);

% 调用局部函数 newton_divided_diff 计算差商表
% 差商表 dd(k) 存储 f[x_0, x_1, ..., x_{k-1}]（即各阶差商的首列）
dd1 = newton_divided_diff(x_nodes1, y_nodes1);
% dd1 是长度为 n1+1 的向量，dd1(1)=f[x_0], dd1(2)=f[x_0,x_1], ...

% 调用局部函数 newton_eval 用差商表在绘图点处求值
% 数值分析：利用 Horner 嵌套乘法高效计算 Newton 插值多项式
y_newton = newton_eval(x_nodes1, dd1, x_plot);

% 计算插值误差：在绘图点处与精确值之差的绝对值
err_newton = abs(y_newton - y_exact);
fprintf('Newton 插值最大绝对误差 = %.4e\n', max(err_newton));
fprintf('（误差主要集中在端点 x=±5 附近，体现 Runge 现象）\n\n');


%% =========================================================
%  第(2)部分：Chebyshev 节点 Lagrange 插值多项式（20次）
%  =========================================================
%
%  【数值分析：Chebyshev 节点与 Lagrange 插值】
%
%  Lagrange 插值公式：
%    L_n(x) = Σ_{k=0}^{n} f(x_k) · l_k(x)
%  其中 Lagrange 基函数：
%    l_k(x) = Π_{j≠k} (x - x_j) / (x_k - x_j)
%
%  插值余项：
%    f(x) - L_n(x) = f^{(n+1)}(ξ)/(n+1)! · Π_{k=0}^{n}(x - x_k)
%  余项大小取决于"节点多项式" ω_{n+1}(x) = Π(x-x_k) 的最大值。
%
%  【Chebyshev 节点的最优性】
%  在 n+1 个节点中，使 max|ω_{n+1}(x)| 最小的节点是 Chebyshev 节点：
%    t_k = cos((2k+1)/(2n+2) · π), k=0,1,...,n  （在[-1,1]上）
%  映射到 [-5,5]：x_k = 5·cos((2k+1)/(2n+2)·π)
%
%  Chebyshev 节点的分布特点：
%  - 在区间两端密集，在中间稀疏（与等距节点相反）
%  - 这种分布能均衡端点处的插值误差，有效压制 Runge 现象
%
%  本题节点公式：x_k = 5·cos((2k+1)/42 · π), k=0,1,...,20
%  即 [-5,5] 上 21 个 Chebyshev 节点（n=20，次数为20）
%  =========================================================

fprintf('--- 第(2)部分：Chebyshev 节点 Lagrange 插值（20次）---\n\n');

n2 = 20;  % 插值多项式次数，节点数 = 21
k2 = 0:n2;  % k = 0,1,...,20，用于生成节点

% 按题目公式生成 Chebyshev 节点
% (2k+1)/42 * pi：这里 42 = 2*(n2+1) = 2*21
% 数值分析：cos((2k+1)/(2n+2)*π) 恰好是 Chebyshev 多项式 T_{n+1}(x) 的 n+1 个根
x_nodes2 = 5 * cos((2*k2 + 1) / 42 * pi);
% 注意：这里不需要 .* 因为 k2 是向量，pi 是标量，MATLAB 自动广播

y_nodes2  = R(x_nodes2);  % 各节点处的函数值

fprintf('Chebyshev 节点：x_k = 5*cos((2k+1)/42*π), k=0,1,...,20\n');
fprintf('节点数：%d，多项式次数：%d\n\n', length(x_nodes2), n2);

% 调用局部函数 lagrange_eval 计算 Lagrange 插值多项式在绘图点处的值
y_lagrange = lagrange_eval(x_nodes2, y_nodes2, x_plot);

err_lagrange = abs(y_lagrange - y_exact);
fprintf('Lagrange 插值最大绝对误差 = %.4e\n', max(err_lagrange));
fprintf('（Chebyshev 节点显著压制了端点振荡）\n\n');


%% =========================================================
%  第(3)部分：等距节点分段线性插值
%  =========================================================
%
%  【数值分析：分段线性插值】
%  分段线性插值是最简单的分段插值方法：
%  在每个子区间 [x_k, x_{k+1}] 上，用连接 (x_k, f_k) 和 (x_{k+1}, f_{k+1})
%  的直线段来近似 f(x)。
%
%  插值公式（在区间 [x_k, x_{k+1}] 上）：
%    P(x) = f_k + (f_{k+1} - f_k)/(x_{k+1} - x_k) · (x - x_k)
%
%  误差估计（设 f 二阶可微）：
%    |f(x) - P(x)| ≤ h²/8 · max|f''(x)|
%  其中 h = max(x_{k+1} - x_k) 是最大步长。
%
%  特点：
%  + 不会出现 Runge 现象（局部方法，每段独立）
%  + 在节点处连续但不可微（一阶导数不连续）
%  - 精度受限于步长，h 不够小时误差较大
%  =========================================================

fprintf('--- 第(3)部分：等距节点分段线性插值 ---\n\n');

% 复用第(1)部分的等距节点（11个节点）
x_nodes3 = x_nodes1;  % [-5,-4,-3,...,4,5]
y_nodes3  = y_nodes1;

% interp1：MATLAB 一维插值函数
% 语法：interp1(x, y, xq, method)
%   x      ：已知节点的 x 坐标（向量）
%   y      ：已知节点的 y 坐标（向量，与 x 等长）
%   xq     ：查询点（要插值的 x 值，可以是向量）
%   method ：插值方法字符串
%     'linear' → 分段线性插值（本部分使用）
%     'spline' → 三次样条插值
%     'pchip'  → 分段三次 Hermite 插值
% 数值分析含义：在每两个相邻节点之间用直线连接
y_piecewise_linear = interp1(x_nodes3, y_nodes3, x_plot, 'linear');
% interp1 已经实现了分段线性插值的核心：对每个查询点找到所属子区间，然后线性插值

err_linear = abs(y_piecewise_linear - y_exact);
fprintf('等距节点（11点）分段线性插值\n');
fprintf('最大绝对误差 = %.4e\n\n', max(err_linear));


%% =========================================================
%  第(4)部分：等距节点三次自然样条插值
%  =========================================================
%
%  【数值分析：三次样条插值】
%  三次样条插值在每个子区间 [x_k, x_{k+1}] 上用三次多项式 S_k(x) 表示，
%  同时要求整体满足：
%    (a) 插值条件：S_k(x_k) = f_k, S_k(x_{k+1}) = f_{k+1}
%    (b) 一阶导数连续：S_{k-1}'(x_k) = S_k'(x_k)
%    (c) 二阶导数连续：S_{k-1}''(x_k) = S_k''(x_k)
%  这些条件使得样条函数在整个区间上具有 C² 连续性（二阶光滑）。
%
%  "自然样条"的边界条件：
%    S''(x_0) = 0,  S''(x_n) = 0
%  即在两端点处二阶导数为零（物理上对应一根自由弯曲的细梁两端无弯矩）。
%
%  三次样条是在节点处满足插值和光滑条件的所有函数中，
%  弯曲能量 ∫[f'']² dx 最小的函数（最优光滑性）。
%
%  MATLAB 中 spline(x, y, xq) 默认使用"非扭结"（not-a-knot）边界条件：
%  要求第一个和最后一个内节点处三阶导数也连续，这通常比自然样条精度更高。
%  要精确使用自然样条边界条件，需要用 csape 函数（需要 Curve Fitting Toolbox），
%  这里我们使用 MATLAB 内置的 spline（not-a-knot），它在实践中表现更好。
%  =========================================================

fprintf('--- 第(4)部分：等距节点三次样条插值 ---\n\n');

x_nodes4 = x_nodes1;  % 复用 11 个等距节点
y_nodes4  = y_nodes1;

% spline(x, y, xq)：三次样条插值
% 语法：yi = spline(x, y, xi)
%   x   ：节点 x 坐标
%   y   ：节点 y 值（也可以包含端点导数信息，见帮助文档）
%   xi  ：查询点
%   返回：在 xi 处的插值结果
% 内部原理：构造分段三次多项式，使整体 C² 连续，采用 not-a-knot 边界条件
% 数值分析含义：在每个子区间上用不同的三次多项式，整体二阶光滑
y_spline = spline(x_nodes4, y_nodes4, x_plot);

err_spline = abs(y_spline - y_exact);
fprintf('等距节点（11点）三次样条插值（not-a-knot 边界条件）\n');
fprintf('最大绝对误差 = %.4e\n\n', max(err_spline));


%% =========================================================
%  第(5)部分：等距节点分段三次 Hermite 插值
%  =========================================================
%
%  【数值分析：分段三次 Hermite 插值（PCHIP）】
%  Hermite 插值不仅要求在节点处匹配函数值，还要求匹配导数值：
%    P(x_k) = f(x_k),    P'(x_k) = f'(x_k)
%
%  分段三次 Hermite 插值（Piecewise Cubic Hermite Interpolating Polynomial）：
%  在每个子区间 [x_k, x_{k+1}] 上用一段三次 Hermite 多项式，
%  Hermite 多项式由两端的函数值和导数值唯一确定（4个条件确定4个参数）。
%
%  PCHIP 的特点：
%    + 整体 C¹ 连续（一阶导数连续，但二阶导数可能不连续）
%    + 保形性（shape-preserving）：若数据单调，则插值函数也单调；
%                                   若数据局部为极值，则插值也保留极值
%    - 光滑性低于三次样条（只有 C¹，不是 C²）
%
%  PCHIP 与三次样条的区别：
%    三次样条：C² 连续，弯曲最小，但可能在极值附近出现"过冲"（不保形）
%    PCHIP：C¹ 连续，保形，但在某些情况下曲线不那么"圆滑"
%
%  MATLAB 的 pchip 函数内部自动估计各节点处的导数（用局部方案），
%  使结果保形且不引入多余振荡。
%  =========================================================

fprintf('--- 第(5)部分：等距节点分段三次 Hermite 插值（PCHIP）---\n\n');

x_nodes5 = x_nodes1;  % 复用 11 个等距节点
y_nodes5  = y_nodes1;

% pchip(x, y, xq)：分段三次 Hermite 插值
% 语法：yi = pchip(x, y, xi)
%   x   ：节点 x 坐标（单调递增）
%   y   ：节点 y 值
%   xi  ：查询点
%   返回：在 xi 处的保形三次 Hermite 插值结果
% 内部原理：在每个子区间上构造三次 Hermite 多项式，
%           节点导数用局部有限差商估计（保形方案）
y_pchip = pchip(x_nodes5, y_nodes5, x_plot);

err_pchip = abs(y_pchip - y_exact);
fprintf('等距节点（11点）分段三次 Hermite 插值（PCHIP）\n');
fprintf('最大绝对误差 = %.4e\n\n', max(err_pchip));


%% =========================================================
%  综合误差对比表
%  =========================================================

fprintf('==============================================\n');
fprintf('  综合误差对比（最大绝对误差）\n');
fprintf('==============================================\n');
fprintf('%-30s  %s\n', '方法', '最大绝对误差');
fprintf('%s\n', repmat('-', 1, 55));
% repmat('-', 1, 55)：将字符 '-' 重复55次，生成分隔线字符串

fprintf('%-30s  %.4e\n', '(1) Newton 等距10次',     max(err_newton));
fprintf('%-30s  %.4e\n', '(2) Lagrange Chebyshev20次', max(err_lagrange));
fprintf('%-30s  %.4e\n', '(3) 分段线性插值',          max(err_linear));
fprintf('%-30s  %.4e\n', '(4) 三次样条插值',           max(err_spline));
fprintf('%-30s  %.4e\n', '(5) 分段三次Hermite(PCHIP)', max(err_pchip));
fprintf('\n');


%% =========================================================
%  图形可视化：5张子图，每种方法单独与精确值对比
%  =========================================================

fprintf('--- 绘制图形 ---\n\n');

% figure：创建新图形窗口
% 'Name'：窗口标题
% 'NumberTitle','off'：不显示"Figure 1"这样的默认标题
% 'Position'：[左边距, 下边距, 宽度, 高度]（单位：像素）
figure('Name', 'Runge函数的五种插值方法', 'NumberTitle', 'off', ...
    'Position', [50, 50, 1400, 900]);
% 窗口较大，以便清晰显示 5 个子图

% 定义方法名称（用于子图标题）
% cell 数组：用 {} 定义，可以存放字符串等不同类型数据
titles = {
    '(1) Newton 插值（等距10次）——Runge现象', ...
    '(2) Lagrange 插值（Chebyshev节点20次）——压制振荡', ...
    '(3) 分段线性插值（等距11点）', ...
    '(4) 三次样条插值（等距11点）', ...
    '(5) 分段三次Hermite插值PCHIP（等距11点）'
};

% 五种方法的插值结果和节点数据，分别打包成 cell 数组方便循环绘图
interp_results = {y_newton, y_lagrange, y_piecewise_linear, y_spline, y_pchip};
node_x = {x_nodes1, x_nodes2, x_nodes3, x_nodes4, x_nodes5};
node_y = {y_nodes1, y_nodes2, y_nodes3, y_nodes4, y_nodes5};
errors = {err_newton, err_lagrange, err_linear, err_spline, err_pchip};

for k = 1:5
    % subplot(m, n, k)：将图窗划分为 m行×n列 的子图网格，当前画第 k 个
    subplot(2, 3, k);
    % 注意：5个子图用2行3列，最后一格(2,3,6)留空

    % plot(x, y, 线型颜色, '属性名', 值, ...)：绘制折线/曲线
    % 'k-'：黑色(k)实线(-)，LineWidth：线宽
    plot(x_plot, y_exact, 'k-', 'LineWidth', 1.5, 'DisplayName', 'R(x)精确值');
    % hold on：保持当前坐标轴，使后续 plot 叠加在同一图上（而非清除重画）
    hold on;

    % 绘制插值结果（蓝色虚线）
    % 'b--'：蓝色(b)虚线(--)
    plot(x_plot, interp_results{k}, 'b--', 'LineWidth', 1.5, 'DisplayName', '插值结果');

    % 绘制插值节点（红色圆圈标记）
    % 'ro'：红色(r)圆圈(o)标记，MarkerSize：标记大小，MarkerFaceColor：标记填充颜色
    plot(node_x{k}, node_y{k}, 'ro', 'MarkerSize', 5, ...
        'MarkerFaceColor', 'r', 'DisplayName', '插值节点');

    % 设置坐标轴显示范围
    % ylim([a, b])：设置 y 轴显示范围为 [a, b]
    % 限制 y 轴范围以便比较，避免 Runge 振荡时 y 轴被拉伸太大
    ylim([-0.5, 1.5]);

    % xlabel, ylabel：坐标轴标签
    xlabel('x');
    ylabel('y');

    % title：子图标题，titles{k} 取 cell 数组第 k 个字符串
    % sprintf：格式化字符串，将误差值嵌入标题中
    title({titles{k}, sprintf('最大误差 = %.2e', max(errors{k}))});
    % cell 数组作为 title 参数时，每个元素占一行

    % legend：图例，'Location','best' 自动选择最佳位置
    legend('Location', 'best', 'FontSize', 7);
    % FontSize：图例字体大小

    % grid on：显示网格线，帮助读图
    grid on;

    hold off;  % 释放 hold，后续不再叠加（对下一个 subplot 从头开始）
end

% 第6个子图位置用于绘制误差对比图
subplot(2, 3, 6);
% semilogy：y 轴使用对数刻度的折线图（适合显示跨越多个数量级的误差）
semilogy(x_plot, errors{1}, 'r-',  'LineWidth', 1.2, 'DisplayName', 'Newton等距10次');
hold on;
semilogy(x_plot, errors{2}, 'b-',  'LineWidth', 1.2, 'DisplayName', 'Lagrange Cheby20次');
semilogy(x_plot, errors{3}, 'g-',  'LineWidth', 1.2, 'DisplayName', '分段线性');
semilogy(x_plot, errors{4}, 'm-',  'LineWidth', 1.2, 'DisplayName', '三次样条');
semilogy(x_plot, errors{5}, 'k--', 'LineWidth', 1.2, 'DisplayName', 'PCHIP');
xlabel('x');
ylabel('绝对误差（对数刻度）');
title('五种方法误差对比');
legend('Location', 'best', 'FontSize', 7);
grid on;
hold off;

% sgtitle：为整个 figure 添加总标题（subplot group title，MATLAB R2018b+）
sgtitle('Runge 函数 R(x) = 1/(1+x²) 的五种插值方法比较', 'FontSize', 13);

fprintf('图形已生成，请查看弹出的图形窗口。\n');
fprintf('\n程序运行完毕。\n');


%% =========================================================
%  局部函数（Local Functions）—— 必须置于脚本末尾
%  =========================================================

function dd = newton_divided_diff(x, y)
    % newton_divided_diff：构造 Newton 插值的差商表
    %
    % 【数值分析：均差（差商）的递推计算】
    % 差商表采用"上三角"存储方式，逐列递推：
    %   dd(i) 最终存储 f[x_0, x_1, ..., x_{i-1}]（即各阶首差商）
    %
    % 递推公式：
    %   f[x_k, ..., x_{k+j}] = (f[x_{k+1},...,x_{k+j}] - f[x_k,...,x_{k+j-1}])
    %                          / (x_{k+j} - x_k)
    %
    % 输入：
    %   x - 节点坐标向量（长度 n+1）
    %   y - 节点函数值向量（长度 n+1）
    % 输出：
    %   dd - 差商向量（长度 n+1），dd(k) = f[x_0,...,x_{k-1}]

    n = length(x);
    dd = y;  % 初始化：dd 先存函数值，即零阶差商 f[x_k] = y_k
    % 注意：这里 dd 是行向量的副本，不影响原始 y

    % 逐列递推，计算各阶差商
    % j：差商的阶数（从 1 到 n-1）
    for j = 1 : n-1
        % 从后往前更新，避免覆盖还未使用的低阶差商
        % k：当前更新的位置（从末尾 n 往前到 j+1）
        for k = n : -1 : j+1
            % 一阶差商：(f[x_k] - f[x_{k-1}]) / (x_k - x_{k-j-1+...})
            % 数组下标从1开始，x(k) 对应 x_{k-1}，x(k-j) 对应 x_{k-j-1}
            dd(k) = (dd(k) - dd(k-1)) / (x(k) - x(k-j));
        end
    end
    % 循环结束后，dd(1) = f[x_0], dd(2) = f[x_0,x_1], ..., dd(n) = f[x_0,...,x_{n-1}]
end


function y_out = newton_eval(x_nodes, dd, x_query)
    % newton_eval：用差商表求值 Newton 插值多项式
    %
    % 【数值分析：Horner 嵌套乘法（秦九韶算法）】
    % Newton 插值多项式：
    %   N_n(x) = dd_0 + dd_1(x-x_0) + dd_2(x-x_0)(x-x_1) + ...
    % 用嵌套形式重写（从最高阶开始）：
    %   N_n(x) = (...((dd_n · (x-x_{n-1}) + dd_{n-1}) · (x-x_{n-2}) + dd_{n-2})...)
    % 这样只需 n 次乘法和 n 次加法，比逐项展开高效。
    %
    % 输入：
    %   x_nodes  - 插值节点坐标（长度 n+1）
    %   dd       - 差商向量（由 newton_divided_diff 计算得到）
    %   x_query  - 查询点（可以是向量）
    % 输出：
    %   y_out    - 在 x_query 处的插值结果（与 x_query 等长）

    n = length(x_nodes) - 1;          % 多项式次数
    y_out = dd(n+1) * ones(size(x_query));  % 初始化为最高阶差商
    % ones(size(x_query))：生成与 x_query 形状相同的全1数组，便于广播

    % 从 n-1 阶到 0 阶，逐步嵌套乘法（Horner 法）
    for k = n : -1 : 1
        % y_out = y_out * (x - x_{k-1}) + dd_k
        % 注意：x_nodes(k) 对应 x_{k-1}（数组下标从1开始）
        y_out = y_out .* (x_query - x_nodes(k)) + dd(k);
        % .* 是逐元素乘法，支持向量化（x_query 可以是向量）
    end
end


function y_out = lagrange_eval(x_nodes, y_nodes, x_query)
    % lagrange_eval：计算 Lagrange 插值多项式在查询点处的值
    %
    % 【数值分析：Lagrange 插值公式】
    % L_n(x) = Σ_{k=0}^{n} y_k · l_k(x)
    % 其中 Lagrange 基函数：
    %   l_k(x) = Π_{j=0, j≠k}^{n} (x - x_j) / (x_k - x_j)
    %
    % 计算复杂度：对每个查询点，计算所有 n+1 个基函数，每个需 n 次乘除法，
    % 总计 O(n²) 次运算（n 较大时较慢，但对于 n=20 无问题）。
    %
    % 输入：
    %   x_nodes  - 插值节点 x 坐标（长度 n+1）
    %   y_nodes  - 插值节点 y 值（长度 n+1）
    %   x_query  - 查询点（向量，长度 m）
    % 输出：
    %   y_out    - 在 x_query 处的 Lagrange 插值结果（长度 m）

    n = length(x_nodes) - 1;             % 多项式次数
    m = length(x_query);                  % 查询点个数
    y_out = zeros(1, m);                  % 预分配输出向量，初始化为0

    % 外循环：遍历每个基函数 l_k(x)
    for k = 1 : n+1
        % 计算第 k 个 Lagrange 基函数在所有查询点处的值
        lk = ones(1, m);  % 初始化 l_k(x) = 1，逐项乘以各因子

        % 内循环：乘以所有 j≠k 的因子 (x - x_j)/(x_k - x_j)
        for j = 1 : n+1
            if j ~= k
                % ~= 是"不等于"运算符
                % 逐元素操作：x_query - x_nodes(j) 是向量，x_nodes(k) - x_nodes(j) 是标量
                lk = lk .* (x_query - x_nodes(j)) / (x_nodes(k) - x_nodes(j));
            end
        end

        % 将第 k 项 y_k * l_k(x) 累加到结果中
        y_out = y_out + y_nodes(k) * lk;
    end
end
