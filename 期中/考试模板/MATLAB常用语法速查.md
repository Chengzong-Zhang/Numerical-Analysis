# MATLAB 常用语法与语句速查

这份文档按“语法 + 示例 + 常见注意点”的方式梳理 MATLAB 中最常用的语法、语句和编程结构。适合用于课程作业、数值计算、画图、数据处理和快速复习。

## 目录

- [1. 基本规则](#1-基本规则)
- [2. 变量与数据类型](#2-变量与数据类型)
- [3. 数组、向量与矩阵](#3-数组向量与矩阵)
- [4. 索引与切片](#4-索引与切片)
- [5. 运算符](#5-运算符)
- [6. 常用数学函数](#6-常用数学函数)
- [7. 流程控制语句](#7-流程控制语句)
- [8. 函数、脚本与匿名函数](#8-函数脚本与匿名函数)
- [9. 字符串与文本处理](#9-字符串与文本处理)
- [10. 单元数组、结构体与表格](#10-单元数组结构体与表格)
- [11. 输入输出与文件操作](#11-输入输出与文件操作)
- [12. 绘图语法](#12-绘图语法)
- [13. 数值计算常用语法](#13-数值计算常用语法)
- [14. 符号计算基础](#14-符号计算基础)
- [15. 程序调试与异常处理](#15-程序调试与异常处理)
- [16. 面向对象基础](#16-面向对象基础)
- [17. 性能与代码规范](#17-性能与代码规范)
- [18. 常用命令汇总](#18-常用命令汇总)

## 1. 基本规则

### 1.1 语句结束

语法：

```matlab
表达式
表达式;
```

说明：

- 不加分号：执行后在命令窗口显示结果。
- 加分号：执行后不显示结果，常用于脚本和函数中。

示例：

```matlab
a = 10
b = 20;
c = a + b
```

输出中会显示 `a` 和 `c`，不会显示 `b`。

### 1.2 注释

语法：

```matlab
% 单行注释

%{
多行注释
多行注释
%}
```

示例：

```matlab
% 计算圆面积
r = 3;
area = pi * r^2;

%{
下面这段代码用于测试
暂时不执行
%}
```

### 1.3 续行符

语法：

```matlab
长表达式第一部分 ...
长表达式第二部分
```

示例：

```matlab
y = 1 + 2 + 3 + ...
    4 + 5 + 6;
```

### 1.4 清理环境

常用语法：

```matlab
clear
clc
close all
clearvars
```

示例：

```matlab
clear;      % 清除工作区变量
clc;        % 清空命令窗口
close all;  % 关闭所有图窗
```

常用脚本开头：

```matlab
clear; clc; close all;
```

### 1.5 查看帮助

语法：

```matlab
help 函数名
doc 函数名
lookfor 关键词
```

示例：

```matlab
help plot
doc eig
lookfor interpolation
```

## 2. 变量与数据类型

### 2.1 变量命名

规则：

- 变量名以字母开头。
- 后面可以包含字母、数字、下划线。
- 区分大小写，`A` 和 `a` 是不同变量。
- 不能使用 MATLAB 关键字。

示例：

```matlab
x = 1;
x1 = 2;
my_data = 3;
```

错误示例：

```matlab
% 1x = 5;      % 不能以数字开头
% my-data = 3; % 不能使用减号
% for = 10;    % for 是关键字
```

查看关键字：

```matlab
iskeyword
```

### 2.2 数值类型

常用语法：

```matlab
a = 3.14;        % double，默认数值类型
b = single(3.14);
c = int32(10);
d = uint8(255);
```

查看类型：

```matlab
class(a)
whos
```

示例：

```matlab
x = 10;
y = single(x);
z = int64(x);

class(x)
class(y)
class(z)
```

### 2.3 逻辑类型

语法：

```matlab
tf = true;
tf = false;
```

示例：

```matlab
x = 5;
flag = x > 0;
```

常见逻辑运算结果：

```matlab
3 > 2       % true
3 == 2      % false
3 ~= 2      % true
```

### 2.4 复数

语法：

```matlab
z = a + b*i
z = a + b*1i
z = complex(a, b)
```

示例：

```matlab
z = 3 + 4i;
abs(z)      % 模长
angle(z)    % 辐角
real(z)     % 实部
imag(z)     % 虚部
conj(z)     % 共轭
```

建议使用 `1i` 或 `1j` 表示虚数单位，避免变量 `i`、`j` 被覆盖。

### 2.5 常用特殊值

语法与含义：

```matlab
pi       % 圆周率
inf      % 正无穷
-inf     % 负无穷
NaN      % 非数值
eps      % 浮点相对精度
```

示例：

```matlab
x = 1 / 0;      % inf
y = 0 / 0;      % NaN
isinf(x)
isnan(y)
```

## 3. 数组、向量与矩阵

MATLAB 的核心数据结构是数组。标量、向量、矩阵都可以看作数组。

### 3.1 创建行向量

语法：

```matlab
v = [元素1 元素2 元素3]
v = [元素1, 元素2, 元素3]
```

示例：

```matlab
v1 = [1 2 3 4];
v2 = [1, 2, 3, 4];
```

### 3.2 创建列向量

语法：

```matlab
v = [元素1; 元素2; 元素3]
```

示例：

```matlab
v = [1; 2; 3; 4];
```

### 3.3 冒号生成等差向量

语法：

```matlab
v = 起点:终点
v = 起点:步长:终点
```

示例：

```matlab
a = 1:5;        % [1 2 3 4 5]
b = 0:0.5:2;    % [0 0.5 1 1.5 2]
c = 5:-1:1;     % [5 4 3 2 1]
```

### 3.4 linspace 和 logspace

语法：

```matlab
x = linspace(a, b, n)
x = logspace(a, b, n)
```

示例：

```matlab
x = linspace(0, 1, 5);   % [0 0.25 0.5 0.75 1]
y = logspace(1, 3, 3);   % [10 100 1000]
```

### 3.5 创建矩阵

语法：

```matlab
A = [第一行; 第二行; 第三行]
```

示例：

```matlab
A = [1 2 3;
     4 5 6;
     7 8 9];
```

### 3.6 常用矩阵生成函数

语法：

```matlab
zeros(m, n)
ones(m, n)
eye(n)
rand(m, n)
randn(m, n)
magic(n)
diag(v)
```

示例：

```matlab
Z = zeros(3, 4);   % 3 行 4 列全 0 矩阵
O = ones(2, 3);    % 2 行 3 列全 1 矩阵
I = eye(4);        % 4 阶单位矩阵
R = rand(2, 2);    % 0 到 1 的均匀随机数
N = randn(2, 2);   % 标准正态随机数
M = magic(3);      % 魔方矩阵
D = diag([1 2 3]); % 对角矩阵
```

### 3.7 数组尺寸

语法：

```matlab
size(A)
size(A, dim)
length(A)
numel(A)
ndims(A)
```

示例：

```matlab
A = [1 2 3; 4 5 6];

size(A)       % [2 3]
size(A, 1)    % 行数 2
size(A, 2)    % 列数 3
length(A)     % 最大维度长度 3
numel(A)      % 元素总数 6
```

### 3.8 矩阵转置

语法：

```matlab
A.'
A'
```

区别：

- `A.'`：普通转置。
- `A'`：共轭转置，复数矩阵会取共轭。

示例：

```matlab
A = [1 2; 3 4];
B = A.';
C = A';

z = [1+2i 3+4i];
z1 = z.';   % 只转置
z2 = z';    % 转置并取共轭
```

### 3.9 拼接数组

语法：

```matlab
C = [A B]      % 水平拼接
C = [A; B]     % 垂直拼接
C = cat(dim, A, B)
```

示例：

```matlab
A = [1 2; 3 4];
B = [5 6; 7 8];

C1 = [A B];
C2 = [A; B];
C3 = cat(3, A, B);
```

## 4. 索引与切片

MATLAB 索引从 1 开始。

### 4.1 单个元素索引

语法：

```matlab
A(行, 列)
```

示例：

```matlab
A = [10 20 30;
     40 50 60];

x = A(1, 2);   % 20
y = A(2, 3);   % 60
```

### 4.2 取整行、整列

语法：

```matlab
A(i, :)
A(:, j)
```

示例：

```matlab
A = [1 2 3;
     4 5 6;
     7 8 9];

row2 = A(2, :);   % 第 2 行
col3 = A(:, 3);   % 第 3 列
```

### 4.3 取子矩阵

语法：

```matlab
A(行索引集合, 列索引集合)
```

示例：

```matlab
A = [1 2 3 4;
     5 6 7 8;
     9 10 11 12];

B = A(1:2, 2:4);
C = A([1 3], [2 4]);
```

### 4.4 end 关键字

语法：

```matlab
A(end)
A(end, :)
A(:, end)
```

示例：

```matlab
v = [10 20 30 40];
last = v(end);        % 40
first3 = v(1:end-1);  % [10 20 30]

A = [1 2 3; 4 5 6];
lastCol = A(:, end);
```

### 4.5 线性索引

MATLAB 按列优先顺序存储矩阵。

语法：

```matlab
A(k)
```

示例：

```matlab
A = [1 2 3;
     4 5 6];

A(1)   % 1
A(2)   % 4
A(3)   % 2
```

### 4.6 逻辑索引

语法：

```matlab
A(逻辑条件)
```

示例：

```matlab
A = [1 5 8 2 9];
B = A(A > 5);        % [8 9]
A(A < 3) = 0;        % 小于 3 的元素改为 0
```

二维示例：

```matlab
A = [1 2 3; 4 5 6];
A(A >= 4) = 100;
```

### 4.7 删除元素、行或列

语法：

```matlab
A(index) = []
A(i, :) = []
A(:, j) = []
```

示例：

```matlab
v = [1 2 3 4];
v(2) = [];       % [1 3 4]

A = [1 2 3;
     4 5 6;
     7 8 9];

A(2, :) = [];    % 删除第 2 行
A(:, 1) = [];    % 删除第 1 列
```

## 5. 运算符

### 5.1 算术运算符

| 运算符 | 含义 | 示例 |
| --- | --- | --- |
| `+` | 加 | `a + b` |
| `-` | 减 | `a - b` |
| `*` | 矩阵乘法 | `A * B` |
| `/` | 右除 | `A / B` |
| `\` | 左除 | `A \ b` |
| `^` | 矩阵幂 | `A ^ 2` |
| `.*` | 元素乘法 | `A .* B` |
| `./` | 元素右除 | `A ./ B` |
| `.\` | 元素左除 | `A .\ B` |
| `.^` | 元素乘方 | `A .^ 2` |

### 5.2 矩阵运算与逐元素运算

示例：

```matlab
A = [1 2; 3 4];
B = [5 6; 7 8];

C1 = A * B;    % 矩阵乘法
C2 = A .* B;   % 对应元素相乘

D1 = A ^ 2;    % A * A
D2 = A .^ 2;   % 每个元素平方
```

注意：

- 线性代数意义的乘法用 `*`。
- 对每个元素分别运算时，通常要加点号，例如 `.*`、`./`、`.^`。

### 5.3 左除与解线性方程组

语法：

```matlab
x = A \ b
```

含义：求解线性方程组 `A*x = b`。

示例：

```matlab
A = [2 1; 1 3];
b = [1; 2];
x = A \ b;
```

不要写成：

```matlab
% x = inv(A) * b;
```

一般推荐使用 `A \ b`，数值稳定性和效率更好。

### 5.4 关系运算符

| 运算符 | 含义 |
| --- | --- |
| `<` | 小于 |
| `<=` | 小于等于 |
| `>` | 大于 |
| `>=` | 大于等于 |
| `==` | 等于 |
| `~=` | 不等于 |

示例：

```matlab
x = [1 2 3 4];
mask = x >= 3;       % [false false true true]
y = x(x >= 3);       % [3 4]
```

### 5.5 逻辑运算符

| 运算符 | 含义 | 用法 |
| --- | --- | --- |
| `&` | 逐元素与 | 数组逻辑 |
| `|` | 逐元素或 | 数组逻辑 |
| `~` | 非 | 数组或标量 |
| `&&` | 短路与 | 标量逻辑，常用于 `if` |
| `||` | 短路或 | 标量逻辑，常用于 `if` |

示例：

```matlab
x = [1 2 3 4 5];
y = x(x > 2 & x < 5);   % [3 4]
```

`if` 中常用：

```matlab
a = 3;
b = 5;

if a > 0 && b > 0
    disp('a 和 b 都是正数');
end
```

## 6. 常用数学函数

### 6.1 初等函数

```matlab
abs(x)       % 绝对值或复数模
sqrt(x)      % 平方根
exp(x)       % e^x
log(x)       % 自然对数
log10(x)     % 以 10 为底的对数
sin(x)       % 正弦，弧度制
cos(x)       % 余弦，弧度制
tan(x)       % 正切，弧度制
asin(x)      % 反正弦
acos(x)      % 反余弦
atan(x)      % 反正切
```

示例：

```matlab
x = pi / 6;
y1 = sin(x);
y2 = cos(x);
y3 = exp(1);
```

角度制函数：

```matlab
sind(30)
cosd(60)
tand(45)
```

### 6.2 取整函数

```matlab
round(x)     % 四舍五入
floor(x)     % 向下取整
ceil(x)      % 向上取整
fix(x)       % 向 0 取整
mod(a, b)    % 取模
rem(a, b)    % 取余
```

示例：

```matlab
x = 3.7;
round(x)     % 4
floor(x)     % 3
ceil(x)      % 4
fix(-3.7)    % -3
mod(10, 3)   % 1
```

### 6.3 最大值、最小值和统计

```matlab
max(x)
min(x)
sum(x)
prod(x)
mean(x)
median(x)
std(x)
var(x)
sort(x)
unique(x)
```

示例：

```matlab
x = [3 1 5 2 5];

maxVal = max(x);
minVal = min(x);
avg = mean(x);
s = sort(x);
u = unique(x);
```

矩阵按维度统计：

```matlab
A = [1 2 3; 4 5 6];

sum(A, 1)    % 每列求和
sum(A, 2)    % 每行求和
mean(A, 1)   % 每列均值
mean(A, 2)   % 每行均值
```

### 6.4 查找函数

语法：

```matlab
idx = find(condition)
[row, col] = find(condition)
```

示例：

```matlab
x = [1 5 3 8 2];
idx = find(x > 3);       % [2 4]

A = [1 0 3; 0 5 0];
[r, c] = find(A ~= 0);
```

### 6.5 判断函数

```matlab
isempty(x)
isnan(x)
isinf(x)
isfinite(x)
isreal(x)
isnumeric(x)
islogical(x)
ischar(x)
isstring(x)
iscell(x)
isstruct(x)
istable(x)
```

示例：

```matlab
x = [1 NaN inf];
isnan(x)
isinf(x)
isfinite(x)
```

## 7. 流程控制语句

### 7.1 if 语句

语法：

```matlab
if 条件
    语句
end
```

示例：

```matlab
x = 5;

if x > 0
    disp('正数');
end
```

### 7.2 if-else 语句

语法：

```matlab
if 条件
    语句1
else
    语句2
end
```

示例：

```matlab
x = -2;

if x >= 0
    y = sqrt(x);
else
    y = NaN;
end
```

### 7.3 if-elseif-else 语句

语法：

```matlab
if 条件1
    语句1
elseif 条件2
    语句2
else
    语句3
end
```

示例：

```matlab
score = 86;

if score >= 90
    grade = 'A';
elseif score >= 80
    grade = 'B';
elseif score >= 70
    grade = 'C';
else
    grade = 'D';
end
```

### 7.4 switch 语句

语法：

```matlab
switch 表达式
    case 值1
        语句1
    case 值2
        语句2
    otherwise
        默认语句
end
```

示例：

```matlab
method = 'linear';

switch method
    case 'linear'
        disp('线性方法');
    case 'spline'
        disp('样条方法');
    otherwise
        disp('未知方法');
end
```

一个 `case` 匹配多个值：

```matlab
day = 'Sat';

switch day
    case {'Sat', 'Sun'}
        disp('周末');
    otherwise
        disp('工作日');
end
```

### 7.5 for 循环

语法：

```matlab
for 变量 = 向量
    循环体
end
```

示例：

```matlab
s = 0;
for k = 1:100
    s = s + k;
end
```

遍历数组：

```matlab
x = [3 6 9];
for val = x
    disp(val);
end
```

### 7.6 while 循环

语法：

```matlab
while 条件
    循环体
end
```

示例：

```matlab
n = 1;
s = 0;

while s < 100
    s = s + n;
    n = n + 1;
end
```

### 7.7 break 与 continue

语法：

```matlab
break
continue
```

示例：

```matlab
for k = 1:10
    if k == 5
        break;       % 直接退出循环
    end
    disp(k);
end
```

```matlab
for k = 1:10
    if mod(k, 2) == 0
        continue;    % 跳过偶数
    end
    disp(k);
end
```

### 7.8 return

语法：

```matlab
return
```

示例：

```matlab
function y = safeSqrt(x)
    if x < 0
        y = NaN;
        return;
    end
    y = sqrt(x);
end
```

## 8. 函数、脚本与匿名函数

### 8.1 脚本文件

脚本是 `.m` 文件，里面直接写命令。脚本共享当前工作区变量。

示例：`main.m`

```matlab
clear; clc;

x = 0:0.1:2*pi;
y = sin(x);
plot(x, y);
```

运行：

```matlab
main
```

### 8.2 函数文件

语法：

```matlab
function 输出变量 = 函数名(输入变量)
    函数体
end
```

示例：`squareValue.m`

```matlab
function y = squareValue(x)
    y = x.^2;
end
```

调用：

```matlab
y = squareValue(3);
```

### 8.3 多输入多输出函数

语法：

```matlab
function [out1, out2] = 函数名(in1, in2)
    函数体
end
```

示例：

```matlab
function [s, p] = sumAndProduct(a, b)
    s = a + b;
    p = a * b;
end
```

调用：

```matlab
[s, p] = sumAndProduct(3, 4);
```

只接收部分输出：

```matlab
s = sumAndProduct(3, 4);
[~, p] = sumAndProduct(3, 4);
```

### 8.4 局部函数

一个 `.m` 文件中可以在主函数或脚本后面定义局部函数。

示例：

```matlab
clear; clc;

x = 3;
y = localSquare(x);
disp(y);

function y = localSquare(x)
    y = x.^2;
end
```

### 8.5 匿名函数

语法：

```matlab
f = @(输入参数) 表达式
```

示例：

```matlab
f = @(x) x.^2 + 2*x + 1;
y = f(3);
```

多参数：

```matlab
g = @(x, y) x.^2 + y.^2;
z = g(3, 4);
```

用于绘图：

```matlab
f = @(x) sin(x) ./ x;
fplot(f, [0.1, 10]);
```

### 8.6 函数句柄

语法：

```matlab
h = @函数名
```

示例：

```matlab
h = @sin;
y = h(pi / 2);
```

传递给数值函数：

```matlab
f = @(x) x.^2 - 2;
root = fzero(f, [1, 2]);
```

### 8.7 nargin 与 nargout

用于判断输入、输出参数个数。

示例：

```matlab
function y = powerValue(x, p)
    if nargin < 2
        p = 2;
    end
    y = x.^p;
end
```

调用：

```matlab
powerValue(3)      % 9
powerValue(3, 3)   % 27
```

## 9. 字符串与文本处理

MATLAB 中常见文本类型有字符数组 `char` 和字符串数组 `string`。

### 9.1 字符数组

语法：

```matlab
s = '文本'
```

示例：

```matlab
s = 'hello';
class(s)       % char
```

### 9.2 字符串

语法：

```matlab
s = "文本"
```

示例：

```matlab
s = "hello";
class(s)       % string
```

### 9.3 拼接字符串

示例：

```matlab
name = "MATLAB";
s = "Hello, " + name;
```

字符数组拼接：

```matlab
s1 = 'Hello, ';
s2 = 'MATLAB';
s = [s1 s2];
```

### 9.4 常用文本函数

```matlab
strlength(s)
contains(s, pattern)
startsWith(s, pattern)
endsWith(s, pattern)
replace(s, old, new)
split(s)
join(s)
lower(s)
upper(s)
strtrim(s)
```

示例：

```matlab
s = "  Hello MATLAB  ";
s1 = strtrim(s);
s2 = lower(s1);
tf = contains(s2, "matlab");
```

### 9.5 sprintf 与 fprintf

语法：

```matlab
s = sprintf(格式, 数据)
fprintf(格式, 数据)
```

常用格式：

| 格式 | 含义 |
| --- | --- |
| `%d` | 整数 |
| `%f` | 小数 |
| `%.2f` | 保留 2 位小数 |
| `%e` | 科学计数法 |
| `%s` | 字符串 |
| `\n` | 换行 |

示例：

```matlab
x = pi;
fprintf('x = %.4f\n', x);

s = sprintf('圆周率约为 %.2f', pi);
```

## 10. 单元数组、结构体与表格

### 10.1 单元数组 cell

单元数组可以存放不同类型、不同尺寸的数据。

语法：

```matlab
C = {元素1, 元素2; 元素3, 元素4}
```

示例：

```matlab
C = {1, 'hello'; [1 2 3], magic(3)};
```

访问单元内容与单元本身：

```matlab
C{1, 2}    % 取出内容，结果是 char
C(1, 2)    % 取出单元，结果仍是 cell
```

修改：

```matlab
C{1, 1} = 100;
C{2, 1} = [4 5 6];
```

### 10.2 结构体 struct

结构体使用字段组织数据。

语法：

```matlab
s.字段名 = 值
```

示例：

```matlab
student.name = "Zhang";
student.id = 1001;
student.score = [88 92 95];
student.meanScore = mean(student.score);
```

创建结构体数组：

```matlab
students(1).name = "A";
students(1).score = 90;
students(2).name = "B";
students(2).score = 85;
```

常用函数：

```matlab
fieldnames(student)
isfield(student, 'name')
rmfield(student, 'id')
```

### 10.3 表格 table

表格适合处理带变量名的数据集。

语法：

```matlab
T = table(变量1, 变量2, ...)
```

示例：

```matlab
Name = ["A"; "B"; "C"];
Score = [90; 85; 92];
Age = [20; 21; 19];

T = table(Name, Age, Score);
```

访问表格：

```matlab
T.Score          % 取 Score 这一列，保留原始类型
T(:, "Score")    % 取子表
T(1, :)          % 第 1 行子表
T{:, "Score"}    % 取出原始数组
```

筛选：

```matlab
T2 = T(T.Score >= 90, :);
```

### 10.4 timetable

时间表适合带时间索引的数据。

示例：

```matlab
Time = datetime(2026, 1, 1) + days(0:2)';
Temp = [10; 12; 11];
TT = timetable(Time, Temp);
```

## 11. 输入输出与文件操作

### 11.1 命令行输入

语法：

```matlab
x = input('提示文字')
s = input('提示文字', 's')
```

示例：

```matlab
n = input('请输入 n: ');
name = input('请输入姓名: ', 's');
```

### 11.2 显示输出

语法：

```matlab
disp(x)
fprintf(format, data)
```

示例：

```matlab
x = 3.14159;
disp(x);
fprintf('x = %.2f\n', x);
```

### 11.3 保存和加载 MAT 文件

语法：

```matlab
save 文件名 变量名
load 文件名
```

示例：

```matlab
x = 1:10;
y = x.^2;
save data.mat x y

clear;
load data.mat
```

保存所有变量：

```matlab
save allData.mat
```

### 11.4 读取和写入文本、表格文件

常用语法：

```matlab
T = readtable('file.csv')
writetable(T, 'file.csv')
A = readmatrix('file.csv')
writematrix(A, 'file.csv')
s = readlines('file.txt')
writelines(s, 'file.txt')
```

示例：

```matlab
A = [1 2 3; 4 5 6];
writematrix(A, 'matrix.csv');
B = readmatrix('matrix.csv');
```

表格示例：

```matlab
Name = ["A"; "B"];
Score = [90; 85];
T = table(Name, Score);
writetable(T, 'scores.csv');
T2 = readtable('scores.csv');
```

### 11.5 低级文件读写

语法：

```matlab
fid = fopen(filename, mode)
fprintf(fid, format, data)
line = fgetl(fid)
fclose(fid)
```

常用模式：

| 模式 | 含义 |
| --- | --- |
| `'r'` | 只读 |
| `'w'` | 写入，覆盖原文件 |
| `'a'` | 追加 |
| `'r+'` | 读写 |

写文件示例：

```matlab
fid = fopen('result.txt', 'w');
fprintf(fid, 'x = %.4f\n', pi);
fclose(fid);
```

逐行读取示例：

```matlab
fid = fopen('result.txt', 'r');
line = fgetl(fid);
while ischar(line)
    disp(line);
    line = fgetl(fid);
end
fclose(fid);
```

### 11.6 路径与文件夹

```matlab
pwd                 % 当前文件夹
cd folderName       % 切换文件夹
dir                 % 列出文件
mkdir folderName    % 创建文件夹
exist(name, 'file') % 判断文件是否存在
fullfile(a, b)      % 拼接路径
addpath(pathName)   % 添加路径
rmpath(pathName)    % 移除路径
```

示例：

```matlab
folder = 'data';
if ~exist(folder, 'dir')
    mkdir(folder);
end

file = fullfile(folder, 'result.csv');
```

## 12. 绘图语法

### 12.1 二维折线图 plot

语法：

```matlab
plot(x, y)
plot(x, y, '样式')
```

示例：

```matlab
x = linspace(0, 2*pi, 100);
y = sin(x);

plot(x, y, 'r-', 'LineWidth', 2);
xlabel('x');
ylabel('sin(x)');
title('正弦函数');
grid on;
```

常用线型和标记：

| 写法 | 含义 |
| --- | --- |
| `'r-'` | 红色实线 |
| `'b--'` | 蓝色虚线 |
| `'ko'` | 黑色圆点 |
| `'g:*'` | 绿色点线加星号 |

### 12.2 多条曲线

示例：

```matlab
x = linspace(0, 2*pi, 100);

plot(x, sin(x), 'r-', x, cos(x), 'b--');
legend('sin(x)', 'cos(x)');
grid on;
```

使用 `hold on`：

```matlab
plot(x, sin(x), 'r-');
hold on;
plot(x, cos(x), 'b--');
hold off;
legend('sin', 'cos');
```

### 12.3 图窗与子图

语法：

```matlab
figure
subplot(m, n, k)
tiledlayout(m, n)
nexttile
```

示例：

```matlab
x = linspace(0, 2*pi, 100);

figure;
subplot(2, 1, 1);
plot(x, sin(x));
title('sin');

subplot(2, 1, 2);
plot(x, cos(x));
title('cos');
```

推荐新语法：

```matlab
tiledlayout(2, 1);

nexttile;
plot(x, sin(x));
title('sin');

nexttile;
plot(x, cos(x));
title('cos');
```

### 12.4 散点图、柱状图、直方图

```matlab
scatter(x, y)
bar(x, y)
histogram(data)
stem(x, y)
stairs(x, y)
```

示例：

```matlab
x = randn(100, 1);
y = 2*x + randn(100, 1);

scatter(x, y, 'filled');
grid on;
```

### 12.5 三维绘图

常用语法：

```matlab
plot3(x, y, z)
mesh(X, Y, Z)
surf(X, Y, Z)
contour(X, Y, Z)
contourf(X, Y, Z)
```

示例：

```matlab
[X, Y] = meshgrid(-2:0.1:2, -2:0.1:2);
Z = X.^2 + Y.^2;

surf(X, Y, Z);
xlabel('x');
ylabel('y');
zlabel('z');
title('z = x^2 + y^2');
shading interp;
colorbar;
```

### 12.6 图形属性

常用语法：

```matlab
xlabel('文本')
ylabel('文本')
zlabel('文本')
title('标题')
legend('曲线1', '曲线2')
grid on
axis equal
axis tight
xlim([xmin xmax])
ylim([ymin ymax])
set(gca, 'FontSize', 12)
```

示例：

```matlab
plot(x, sin(x), 'LineWidth', 2);
xlabel('x');
ylabel('y');
title('y = sin(x)');
legend('sin(x)', 'Location', 'best');
grid on;
xlim([0, 2*pi]);
set(gca, 'FontSize', 12);
```

### 12.7 保存图片

语法：

```matlab
saveas(gcf, 'figure.png')
exportgraphics(gcf, 'figure.png', 'Resolution', 300)
```

示例：

```matlab
plot(x, sin(x));
exportgraphics(gcf, 'sin.png', 'Resolution', 300);
```

## 13. 数值计算常用语法

### 13.1 多项式

MATLAB 用系数向量表示多项式，按降幂排列。

例如：

```matlab
% p(x) = x^2 - 3x + 2
p = [1 -3 2];
```

常用函数：

```matlab
polyval(p, x)    % 计算多项式值
roots(p)         % 求根
polyfit(x, y, n) % n 次多项式拟合
conv(p, q)       % 多项式乘法
deconv(p, q)     % 多项式除法
```

示例：

```matlab
x = 0:0.1:5;
y = 2*x.^2 - 3*x + 1 + randn(size(x));

p = polyfit(x, y, 2);
yfit = polyval(p, x);

plot(x, y, 'o', x, yfit, '-');
legend('data', 'fit');
```

### 13.2 解方程

单变量非线性方程：

```matlab
f = @(x) x.^2 - 2;
root = fzero(f, [1, 2]);
```

多变量非线性方程通常使用 Optimization Toolbox 的 `fsolve`：

```matlab
F = @(x) [x(1)^2 + x(2)^2 - 1;
          x(1) - x(2)];
x0 = [0.5; 0.5];
sol = fsolve(F, x0);
```

### 13.3 数值积分

语法：

```matlab
q = integral(f, a, b)
q = integral2(f, xmin, xmax, ymin, ymax)
q = trapz(x, y)
```

示例：

```matlab
f = @(x) sin(x);
q = integral(f, 0, pi);
```

离散数据积分：

```matlab
x = linspace(0, pi, 100);
y = sin(x);
q = trapz(x, y);
```

### 13.4 数值微分

常用语法：

```matlab
diff(y)
gradient(y, h)
```

示例：

```matlab
x = linspace(0, 2*pi, 100);
y = sin(x);
h = x(2) - x(1);

dy1 = diff(y) / h;
dy2 = gradient(y, h);
```

### 13.5 插值

语法：

```matlab
yq = interp1(x, y, xq, method)
vq = interp2(X, Y, V, Xq, Yq, method)
```

常用方法：

| 方法 | 含义 |
| --- | --- |
| `'nearest'` | 最近邻 |
| `'linear'` | 线性插值 |
| `'spline'` | 三次样条 |
| `'pchip'` | 分段三次 Hermite |

示例：

```matlab
x = 0:1:10;
y = sin(x);
xq = 0:0.1:10;

yq = interp1(x, y, xq, 'spline');
plot(x, y, 'o', xq, yq, '-');
```

### 13.6 常微分方程 ODE

语法：

```matlab
[t, y] = ode45(odefun, tspan, y0)
```

示例：求解 `y' = -2y, y(0)=1`

```matlab
f = @(t, y) -2*y;
tspan = [0 5];
y0 = 1;

[t, y] = ode45(f, tspan, y0);
plot(t, y);
```

二阶方程转一阶系统示例：

方程：

```text
y'' + y = 0
```

令：

```text
u1 = y
u2 = y'
```

则：

```text
u1' = u2
u2' = -u1
```

代码：

```matlab
f = @(t, u) [u(2); -u(1)];
tspan = [0, 10];
u0 = [1; 0];

[t, u] = ode45(f, tspan, u0);
plot(t, u(:, 1));
```

### 13.7 线性代数

常用函数：

```matlab
det(A)       % 行列式
rank(A)      % 秩
trace(A)     % 迹
inv(A)       % 逆矩阵
eig(A)       % 特征值
[V, D] = eig(A)  % 特征向量和特征值
norm(A)      % 范数
cond(A)      % 条件数
qr(A)        % QR 分解
lu(A)        % LU 分解
svd(A)       % 奇异值分解
```

示例：

```matlab
A = [2 1; 1 2];
b = [1; 0];

x = A \ b;
lambda = eig(A);
[V, D] = eig(A);
```

### 13.8 随机数

常用语法：

```matlab
rng(seed)
rand(m, n)
randn(m, n)
randi([a b], m, n)
randperm(n)
```

示例：

```matlab
rng(1);
A = rand(3, 3);
B = randn(3, 3);
C = randi([1, 10], 2, 4);
p = randperm(10);
```

## 14. 符号计算基础

需要 Symbolic Math Toolbox。

### 14.1 定义符号变量

语法：

```matlab
syms x y
```

示例：

```matlab
syms x
f = x^2 + 2*x + 1;
```

### 14.2 符号化简、展开、因式分解

```matlab
simplify(expr)
expand(expr)
factor(expr)
collect(expr)
```

示例：

```matlab
syms x
f = (x + 1)^2;
expand(f)
factor(x^2 - 1)
```

### 14.3 求导、积分、极限

语法：

```matlab
diff(f, x)
int(f, x)
int(f, x, a, b)
limit(f, x, a)
```

示例：

```matlab
syms x
f = sin(x)^2;

df = diff(f, x);
F = int(f, x);
L = limit(sin(x)/x, x, 0);
```

### 14.4 解符号方程

语法：

```matlab
sol = solve(equation, variable)
```

示例：

```matlab
syms x
sol = solve(x^2 - 2 == 0, x);
```

方程组：

```matlab
syms x y
[solx, soly] = solve([x + y == 1, x - y == 3], [x, y]);
```

## 15. 程序调试与异常处理

### 15.1 try-catch

语法：

```matlab
try
    可能出错的语句
catch ME
    出错后的处理
end
```

示例：

```matlab
try
    A = readmatrix('not_exist.csv');
catch ME
    fprintf('读取失败: %s\n', ME.message);
end
```

### 15.2 error、warning、assert

语法：

```matlab
error('错误信息')
warning('警告信息')
assert(条件, '错误信息')
```

示例：

```matlab
function y = positiveSqrt(x)
    assert(x >= 0, 'x 必须非负');
    y = sqrt(x);
end
```

### 15.3 调试命令

```matlab
dbstop if error      % 出错时停止
dbclear all          % 清除断点
dbcont               % 继续运行
dbquit               % 退出调试
keyboard             % 在代码中手动进入调试状态
```

常用界面操作：

- 在编辑器左侧点击可设置断点。
- `F10` 单步执行。
- `F11` 进入函数。
- `Shift + F11` 跳出函数。

## 16. 面向对象基础

MATLAB 支持使用 `classdef` 定义类。

### 16.1 类的基本语法

示例：`Point2D.m`

```matlab
classdef Point2D
    properties
        X
        Y
    end

    methods
        function obj = Point2D(x, y)
            obj.X = x;
            obj.Y = y;
        end

        function d = distanceToOrigin(obj)
            d = sqrt(obj.X^2 + obj.Y^2);
        end
    end
end
```

调用：

```matlab
p = Point2D(3, 4);
d = p.distanceToOrigin();
```

### 16.2 属性访问控制

语法：

```matlab
properties
end

properties (Access = private)
end
```

示例：

```matlab
classdef Counter
    properties (Access = private)
        Value = 0
    end

    methods
        function obj = add(obj)
            obj.Value = obj.Value + 1;
        end

        function v = getValue(obj)
            v = obj.Value;
        end
    end
end
```

### 16.3 值类与句柄类

默认类是值类，赋值时会复制对象。继承 `handle` 后是句柄类，赋值时共享对象。

语法：

```matlab
classdef MyClass < handle
end
```

示例：

```matlab
classdef MyHandleClass < handle
    properties
        Value
    end
end
```

## 17. 性能与代码规范

### 17.1 预分配数组

不推荐：

```matlab
for k = 1:10000
    x(k) = k^2;
end
```

推荐：

```matlab
x = zeros(1, 10000);
for k = 1:10000
    x(k) = k^2;
end
```

### 17.2 向量化计算

不推荐：

```matlab
x = 0:0.01:10;
y = zeros(size(x));
for k = 1:length(x)
    y(k) = sin(x(k)) + x(k)^2;
end
```

推荐：

```matlab
x = 0:0.01:10;
y = sin(x) + x.^2;
```

### 17.3 避免覆盖内置函数名

不推荐：

```matlab
sum = 10;
mean = 3;
plot = 5;
```

推荐：

```matlab
totalValue = 10;
avgValue = 3;
plotIndex = 5;
```

如果不小心覆盖了内置函数，可以清除变量：

```matlab
clear sum
```

### 17.4 脚本与函数的选择

建议：

- 临时计算、课程演示：可以用脚本。
- 可复用代码、复杂程序：建议写函数。
- 主程序可以写脚本，核心算法写函数。

### 17.5 常见代码结构

脚本模板：

```matlab
clear; clc; close all;

%% 参数设置
n = 100;

%% 数据生成
x = linspace(0, 1, n);
y = sin(2*pi*x);

%% 计算
result = mean(y);

%% 绘图
plot(x, y);
grid on;
title('Result');
```

函数模板：

```matlab
function output = myFunction(input)
%MYFUNCTION 简要说明函数功能
%   output = MYFUNCTION(input) 详细说明输入输出

    arguments
        input double
    end

    output = input.^2;
end
```

## 18. 常用命令汇总

### 18.1 工作区与命令窗口

```matlab
clear          % 清除变量
clearvars      % 清除工作区变量
clc            % 清空命令窗口
who            % 显示变量名
whos           % 显示变量详细信息
format long    % 长格式显示
format short   % 短格式显示
```

### 18.2 文件与路径

```matlab
pwd
cd
dir
ls
mkdir
rmdir
delete
copyfile
movefile
fullfile
addpath
genpath
savepath
```

### 18.3 数组与矩阵

```matlab
zeros
ones
eye
rand
randn
randi
size
length
numel
reshape
repmat
permute
squeeze
diag
tril
triu
```

示例：

```matlab
A = reshape(1:12, 3, 4);
B = repmat([1 2], 3, 1);
C = squeeze(rand(1, 3, 4));
```

### 18.4 排序与集合

```matlab
sort
sortrows
unique
union
intersect
setdiff
ismember
```

示例：

```matlab
a = [1 2 3 4];
b = [3 4 5 6];

u = union(a, b);
i = intersect(a, b);
d = setdiff(a, b);
tf = ismember(a, b);
```

### 18.5 日期与时间

```matlab
datetime
duration
calendarDuration
dateshift
year
month
day
hours
minutes
seconds
```

示例：

```matlab
t = datetime('now');
t2 = t + days(7);
y = year(t);
m = month(t);
d = day(t);
```

### 18.6 常见一行写法

条件筛选：

```matlab
x = [1 5 2 8 3];
y = x(x > 3);
```

替换异常值：

```matlab
x(isnan(x)) = 0;
```

归一化：

```matlab
xNorm = (x - min(x)) / (max(x) - min(x));
```

均值方差标准化：

```matlab
xStd = (x - mean(x)) / std(x);
```

按行求最大值：

```matlab
rowMax = max(A, [], 2);
```

按列求最大值：

```matlab
colMax = max(A, [], 1);
```

找到最大值位置：

```matlab
[maxVal, idx] = max(x);
```

矩阵中满足条件的元素个数：

```matlab
count = nnz(A > 0);
```

## 附录：常见错误与修正

### 错误 1：矩阵乘法和元素乘法混淆

错误：

```matlab
x = 0:0.1:1;
y = x^2;
```

修正：

```matlab
y = x.^2;
```

### 错误 2：索引从 0 开始

错误：

```matlab
x = [10 20 30];
first = x(0);
```

修正：

```matlab
first = x(1);
```

### 错误 3：if 条件不是标量逻辑

错误：

```matlab
x = [1 2 3];
if x > 0
    disp('positive');
end
```

修正：

```matlab
if all(x > 0)
    disp('all positive');
end
```

或：

```matlab
if any(x > 0)
    disp('at least one positive');
end
```

### 错误 4：字符串比较使用 `==`

对于 `string` 可以使用 `==`：

```matlab
s = "linear";
if s == "linear"
    disp('ok');
end
```

对于字符数组推荐使用 `strcmp`：

```matlab
s = 'linear';
if strcmp(s, 'linear')
    disp('ok');
end
```

### 错误 5：函数文件名与函数名不一致

函数：

```matlab
function y = mySquare(x)
    y = x.^2;
end
```

文件名应为：

```text
mySquare.m
```

### 错误 6：忘记预分配导致程序变慢

错误：

```matlab
for k = 1:100000
    a(k) = k;
end
```

修正：

```matlab
a = zeros(1, 100000);
for k = 1:100000
    a(k) = k;
end
```

### 错误 7：用 `inv(A)*b` 解线性方程

不推荐：

```matlab
x = inv(A) * b;
```

推荐：

```matlab
x = A \ b;
```

## 附录：一个完整小例子

下面的示例综合使用变量、向量化、函数、绘图、文件保存等常用语法。

```matlab
clear; clc; close all;

%% 生成数据
x = linspace(0, 2*pi, 200);
y = sin(x) + 0.1 * randn(size(x));

%% 平滑处理
windowSize = 5;
ySmooth = movmean(y, windowSize);

%% 计算误差
yTrue = sin(x);
rmse = sqrt(mean((ySmooth - yTrue).^2));
fprintf('RMSE = %.4f\n', rmse);

%% 绘图
figure;
plot(x, y, '.', 'DisplayName', 'Noisy data');
hold on;
plot(x, ySmooth, 'r-', 'LineWidth', 2, 'DisplayName', 'Smoothed');
plot(x, yTrue, 'k--', 'LineWidth', 1.5, 'DisplayName', 'True');
hold off;

xlabel('x');
ylabel('y');
title('Signal smoothing example');
legend('Location', 'best');
grid on;

%% 保存结果
T = table(x(:), y(:), ySmooth(:), yTrue(:), ...
    'VariableNames', {'x', 'Noisy', 'Smooth', 'True'});
writetable(T, 'signal_result.csv');
exportgraphics(gcf, 'signal_plot.png', 'Resolution', 300);
```

