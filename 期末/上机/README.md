# 第五章：数值 ODE 上机文件索引

本目录依据 `数值ODE/第五章.pdf` 整理。旧文件仅移动分类，未修改内容；新增 MATLAB 模板均可直接运行。

## 最快使用方式

1. 打开对应 `.m` 文件。
2. 搜索 `【必须替换】`，只修改题目给出的函数、区间、初值、步长或容差。
3. 先运行默认例题确认模板正常，再替换为考试题。
4. 隐式方法若不收敛，减小 `h`、增大 `maxIter`，或改用 Newton 迭代模板。

## 算法与模板对应表

| 教材算法 | 模板 |
|---|---|
| Taylor 级数法 | `02_MATLAB算法模板/01_单步法/T01_Taylor_Series.m` |
| 向前 Euler、向后 Euler、梯形、改进 Euler、中点、二阶 Heun | `02_MATLAB算法模板/01_单步法/T02_Euler_RK2_Family.m` |
| 一般显式 RK、三阶 Heun、三阶 RK、经典 RK4 | `02_MATLAB算法模板/02_Runge_Kutta/T03_Explicit_Runge_Kutta.m` |
| 一级二阶、二级四阶、三级六阶隐式 Gauss-RK | `02_MATLAB算法模板/02_Runge_Kutta/T04_Implicit_Gauss_RK.m` |
| 自适应 Runge-Kutta-Fehlberg 4(5) | `02_MATLAB算法模板/02_Runge_Kutta/T05_Adaptive_RKF45.m` |
| Adams-Bashforth 1 至 6 步 | `02_MATLAB算法模板/03_线性多步法/T06_Adams_Bashforth.m` |
| Adams-Moulton 1 至 6 阶、Adams PECE | `02_MATLAB算法模板/03_线性多步法/T07_Adams_Moulton_PECE.m` |
| Milne、Hamming、Nyström | `02_MATLAB算法模板/03_线性多步法/T08_Milne_Hamming_Nystrom.m` |
| BDF/Gear 1 至 6 阶 | `02_MATLAB算法模板/03_线性多步法/T09_BDF_Gear.m` |
| 方程组、高阶方程化一阶系统、RK4 | `02_MATLAB算法模板/04_方程组与高阶方程/T10_ODE_System_Higher_Order_RK4.m` |
| 不动点迭代、Newton 迭代 | `02_MATLAB算法模板/04_方程组与高阶方程/T11_Implicit_Equation_Iterators.m` |
| 最大误差、收敛阶、Richardson 外推、稳定函数 | `02_MATLAB算法模板/05_误差与稳定性/T12_Error_Order_Stability.m` |
| 一般线性多步法阶数条件、根条件、绝对稳定判别 | `02_MATLAB算法模板/05_误差与稳定性/T13_Linear_Multistep_Order_Root_Check.m` |

运行 `02_MATLAB算法模板/verify_all_templates.m` 可一次验证全部模板的默认例题。

## 文件结构

- `01_原有资料`：原有 Python 模板、讲义、图片和编译中间文件。
- `02_MATLAB算法模板`：本次新增的第五章全部算法模板。
