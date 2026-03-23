%% =========================================================
%  charts_demo.m
%  功能：演示 MATLAB 中最常用的图表类型
%  包含：折线图、散点图、柱状图、直方图、饼图、3D曲面图
%  =========================================================

clear; close all; clc;

%% =========================================================
%  【MATLAB 绘图基础语法速查】
%
%  figure               → 新建一个图形窗口
%  subplot(m, n, k)     → 将图窗划分为 m行×n列，激活第k个子图
%  plot(x, y, '线型')   → 折线图（最基本的绘图函数）
%  xlabel('文字')       → x轴标签
%  ylabel('文字')       → y轴标签
%  title('文字')        → 图标题
%  legend('图例1',...)  → 添加图例
%  grid on/off          → 显示/隐藏网格
%  hold on/off          → 保持/释放当前坐标轴（叠加绘图）
%  xlim([a b])          → 设置x轴范围
%  ylim([a b])          → 设置y轴范围
%
%  常用线型/颜色/标记：
%    颜色：r红 g绿 b蓝 k黑 m紫 c青 y黄
%    线型：-实线  --虚线  :点线  -.点划线
%    标记：o圆  +加  *星  .点  x叉  s方  d菱  ^三角
%    组合：'r--o' 表示红色虚线加圆圈标记
%  =========================================================


%% =========================================================
%  图表1：折线图 (Line Plot)
%  适用场景：展示连续变化的趋势，如函数曲线、时间序列
%  =========================================================

figure('Name', '图1: 折线图', 'NumberTitle', 'off', 'Position', [50 500 700 450]);

% 生成数据
% linspace(a, b, n)：在 [a, b] 区间内均匀生成 n 个点，返回行向量
x = linspace(0, 2*pi, 200);   % 0 到 2π 之间取200个点
% pi：MATLAB 内置常量，= 3.14159...

% MATLAB 支持向量化运算：对向量 x 的每个元素分别计算 sin/cos，结果也是向量
y1 = sin(x);        % 逐元素计算 sin
y2 = cos(x);        % 逐元素计算 cos
y3 = sin(x) .* cos(x);  % .*：逐元素乘法（区别于矩阵乘法 *）

% 绘制三条曲线
plot(x, y1, 'b-',  'LineWidth', 2);   % 蓝色实线，线宽2
hold on;                               % 保持坐标轴，继续叠加绘图
plot(x, y2, 'r--', 'LineWidth', 2);   % 红色虚线
plot(x, y3, 'g:',  'LineWidth', 1.5); % 绿色点线

% 设置坐标轴范围
xlim([0, 2*pi]);
ylim([-1.2, 1.2]);

% 设置 x 轴刻度为 π 的倍数（更美观）
% set(gca, ...)：设置当前坐标轴(gca=get current axes)的属性
set(gca, 'XTick', [0, pi/2, pi, 3*pi/2, 2*pi], ...
         'XTickLabel', {'0', 'π/2', 'π', '3π/2', '2π'}, ...
         'FontSize', 11);
% 'XTick'：刻度位置（数值数组）
% 'XTickLabel'：对应的显示标签（cell数组）

xlabel('x (弧度)');
ylabel('y');
title('折线图：三角函数', 'FontSize', 13);

% legend：添加图例，字符串顺序对应 plot 的顺序
legend('sin(x)', 'cos(x)', 'sin(x)·cos(x)', 'Location', 'best');
% 'Location','best'：自动选择不遮挡曲线的位置

grid on;
hold off;  % 释放坐标轴，恢复正常（新 plot 会覆盖当前图）


%% =========================================================
%  图表2：散点图 (Scatter Plot)
%  适用场景：展示两变量之间的关系/分布，数据点离散
%  =========================================================

figure('Name', '图2: 散点图', 'NumberTitle', 'off', 'Position', [760 500 700 450]);

% 生成随机数据
% randn(m, n)：生成 m行×n列 的标准正态分布随机矩阵（均值0，标准差1）
rng(42);   % rng：设置随机数种子，确保每次运行结果相同（可复现）
n_points = 200;   % 数据点数量

% 生成两组不同分布的点
% 第一组：中心在(1,1)，标准差0.5
x1 = 1 + 0.5 * randn(n_points, 1);   % randn(200,1)：200行1列的列向量
y1 = 1 + 0.5 * randn(n_points, 1);

% 第二组：中心在(3,2)，x方向标准差1，y方向0.3
x2 = 3 + 1.0 * randn(n_points, 1);
y2 = 2 + 0.3 * randn(n_points, 1);

% scatter(x, y, 大小, 颜色)：散点图
% 第三参数：点的大小（面积，单位 points²）
% 第四参数：颜色（可以是颜色字符串，也可以是RGB三元组）
scatter(x1, y1, 30, 'b', 'filled');  % 蓝色实心点，大小30
% 'filled'：填充颜色（否则是空心圆圈）
hold on;
scatter(x2, y2, 30, 'r', 'filled');  % 红色实心点

% 可以用 RGB 三元组指定颜色：[R, G, B]，每个分量取 0~1
% 例：[0.2, 0.8, 0.3] 表示偏绿色

xlabel('X 值');
ylabel('Y 值');
title('散点图：两组正态分布数据', 'FontSize', 13);
legend('第一组 (中心1,1)', '第二组 (中心3,2)', 'Location', 'northwest');
grid on;
hold off;


%% =========================================================
%  图表3：柱状图 (Bar Chart)
%  适用场景：比较离散类别之间的数量大小
%  =========================================================

figure('Name', '图3: 柱状图', 'NumberTitle', 'off', 'Position', [50 30 700 450]);

% 准备数据：某班级各科成绩（行=学生，列=科目）
% 矩阵定义：用 [] 包裹，逗号/空格分隔列，分号分隔行
scores = [85, 92, 78, 95;   % 学生1：语文、数学、英语、物理
          76, 88, 90, 82;   % 学生2
          91, 79, 85, 88;   % 学生3
          68, 95, 72, 91];  % 学生4
% scores 是 4×4 矩阵

% bar(data)：当 data 是矩阵时，每行一组，每列一种颜色（分组柱状图）
bar(scores);

% 设置 x 轴标签为学生名称
set(gca, 'XTickLabel', {'学生1', '学生2', '学生3', '学生4'}, 'FontSize', 11);

xlabel('学生');
ylabel('分数');
title('分组柱状图：各学生各科成绩', 'FontSize', 13);
legend('语文', '数学', '英语', '物理', 'Location', 'best');
ylim([0, 110]);   % y 轴从0到110，留出文字标注空间
grid on;

% 在柱子顶部添加数值标注
% size(scores)：返回矩阵的 [行数, 列数]
[n_students, n_subjects] = size(scores);   % 解包为两个变量
for i = 1:n_students
    for j = 1:n_subjects
        % 计算每根柱子的 x 坐标：需要根据分组位置计算
        % 分组柱状图中，第i组第j根柱的x坐标约为 i + (j - (n_subjects+1)/2)*0.2
        x_pos = i + (j - (n_subjects + 1) / 2) * 0.19;
        % text(x, y, '文字')：在坐标(x,y)处添加文字
        text(x_pos, scores(i,j) + 1.5, num2str(scores(i,j)), ...
            'HorizontalAlignment', 'center', 'FontSize', 8);
        % num2str()：将数字转换为字符串（text函数需要字符串参数）
        % 'HorizontalAlignment','center'：文字水平居中对齐
    end
end


%% =========================================================
%  图表4：直方图 (Histogram)
%  适用场景：展示数据的频率分布/概率密度
%  =========================================================

figure('Name', '图4: 直方图', 'NumberTitle', 'off', 'Position', [760 30 700 450]);

% 生成1000个正态分布随机数
data = 5 + 2 * randn(1000, 1);   % 均值5，标准差2

% histogram(data, n)：绘制直方图，n 为分组数（bin数量）
h = histogram(data, 30);
% h 是返回的图形对象，可以通过它修改属性

% 修改直方图外观
h.FaceColor = [0.3, 0.6, 0.9];   % 填充颜色（浅蓝色），RGB格式
h.EdgeColor = 'white';             % 边框颜色
h.FaceAlpha = 0.8;                 % 透明度：0=全透明，1=不透明

hold on;

% 叠加理论正态分布曲线
x_fit = linspace(min(data)-1, max(data)+1, 300);
% normpdf(x, mu, sigma)：正态分布概率密度函数
% 需要将 pdf 乘以总样本数和 bin 宽度，才能与频数直方图匹配
bin_width = h.BinWidth;   % 获取自动计算的 bin 宽度
y_fit = normpdf(x_fit, 5, 2) * 1000 * bin_width;

plot(x_fit, y_fit, 'r-', 'LineWidth', 2.5);   % 红色曲线

xlabel('数值');
ylabel('频数');
title(sprintf('直方图：正态分布 N(5, 2²)，n=%d', length(data)), 'FontSize', 13);
% sprintf：生成格式化字符串（%d插入变量值），作为标题文字
legend('样本频数', '理论密度曲线', 'Location', 'best');
grid on;
hold off;


%% =========================================================
%  图表5：饼图 (Pie Chart)
%  适用场景：展示各部分占整体的比例
%  =========================================================

figure('Name', '图5: 饼图', 'NumberTitle', 'off', 'Position', [400 200 600 450]);

% 数据：各编程语言使用比例（百分比，总和不必为100，pie会自动归一化）
lang_data   = [35, 25, 18, 12, 10];
lang_labels = {'Python', 'JavaScript', 'Java', 'C++', '其他'};

% explode：控制哪个扇区"爆炸"弹出（突出显示）
% 1=弹出，0=不弹出，与 lang_data 长度相同
explode = [1, 0, 0, 0, 0];   % 只突出显示第1个（Python）

% pie(data, explode, labels)：绘制饼图
pie(lang_data, explode, lang_labels);

title('饼图：编程语言使用比例', 'FontSize', 13);

% 设置配色方案
% colormap：设置整个图形的颜色映射
% parula/jet/hsv/hot/cool 是 MATLAB 内置的颜色映射方案
colormap(gca, parula(length(lang_data)));
% gca：当前坐标轴；parula(n)：生成n种颜色的 parula 配色


%% =========================================================
%  图表6：3D 曲面图 (Surface Plot)
%  适用场景：展示双变量函数 z = f(x, y) 的三维形状
%  =========================================================

figure('Name', '图6: 3D曲面图', 'NumberTitle', 'off', 'Position', [200 100 800 500]);

% 生成网格数据
% meshgrid(x, y)：根据向量 x, y 生成网格矩阵 X, Y
% X 的每一行都是 x，Y 的每一列都是 y
x_range = linspace(-3, 3, 60);   % x 方向60个点
y_range = linspace(-3, 3, 60);   % y 方向60个点
[X, Y] = meshgrid(x_range, y_range);   % 生成 60×60 的网格

% 计算 z 值：二维高斯函数（帽形曲面）
Z = exp(-(X.^2 + Y.^2) / 2);
% .^：逐元素乘方（X.^2 表示矩阵 X 每个元素平方）
% exp()：以 e 为底的指数函数，逐元素计算

% subplot 分两个子图
subplot(1, 2, 1);

% surf(X, Y, Z)：绘制3D曲面（带颜色填充）
surf(X, Y, Z);
% shading interp：平滑颜色过渡（去掉网格线，颜色插值）
shading interp;
colorbar;         % 显示颜色条（色标）
colormap(jet);    % 使用 jet 颜色映射（蓝→绿→红）
xlabel('X');  ylabel('Y');  zlabel('Z');
title('surf：高斯曲面', 'FontSize', 12);

subplot(1, 2, 2);

% contourf(X, Y, Z, n)：等高线填充图（俯视角度）
% n：等高线条数
[C, h_cont] = contourf(X, Y, Z, 15);
% C：等高线数据矩阵；h_cont：图形对象（可用于后续修改）
clabel(C, h_cont, 'FontSize', 8);   % 在等高线上显示数值标签
colorbar;
colormap(jet);
xlabel('X');  ylabel('Y');
title('contourf：等高线图', 'FontSize', 12);

% sgtitle：整个 figure 的总标题
sgtitle('3D曲面图与等高线图：高斯函数', 'FontSize', 14);


%% =========================================================
%  图表7：误差棒图 (Error Bar Plot)
%  适用场景：展示数据的均值和不确定性（标准差/置信区间）
%  在科学实验结果展示中非常常用
%  =========================================================

figure('Name', '图7: 误差棒图', 'NumberTitle', 'off', 'Position', [100 200 600 420]);

% 模拟实验数据：5个测量点，每点多次重复测量
x_exp = 1:5;                              % 实验编号
y_mean = [2.1, 3.5, 2.8, 4.2, 3.9];      % 各点均值（行向量）
y_std  = [0.3, 0.5, 0.2, 0.6, 0.4];      % 各点标准差

% errorbar(x, y, err)：绘制带误差棒的图
% err：误差范围（上下对称），即在 y ± err 处画横线
errorbar(x_exp, y_mean, y_std, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
% 'bo-'：蓝色圆圈实线
% 'MarkerFaceColor','b'：圆圈内填充蓝色

xlabel('实验编号');
ylabel('测量值');
title('误差棒图：实验测量结果（均值 ± 标准差）', 'FontSize', 12);
xlim([0.5, 5.5]);
ylim([0, 6]);
grid on;


fprintf('所有图表已生成，共7个图形窗口。\n');
fprintf('\n【常见图表选择指南】\n');
fprintf('  plot/semilogy  → 折线图：连续趋势、函数曲线\n');
fprintf('  scatter        → 散点图：两变量关系、聚类分布\n');
fprintf('  bar/barh       → 柱状图：离散类别比较（barh为水平柱状）\n');
fprintf('  histogram      → 直方图：频率分布、数据统计\n');
fprintf('  pie            → 饼图：占比关系\n');
fprintf('  surf/mesh      → 3D曲面：双变量函数可视化\n');
fprintf('  contourf       → 等高线：俯视角的曲面\n');
fprintf('  errorbar       → 误差棒：科学实验数据\n');
