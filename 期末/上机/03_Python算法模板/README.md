# 第五章数值 ODE Python 算法模板

本目录与 `02_MATLAB算法模板` 一一对应。每个 Python 模板均可独立运行，每个可执行代码行都附有中文语法和数值分析注释。

## 使用方式

1. 打开对应 `.py` 文件。
2. 搜索 `【必须替换】`，修改题目给出的函数、区间、初值、步长或容差。
3. 在 PowerShell 中运行：

```powershell
python ".\期末\上机\03_Python算法模板\01_单步法\T02_Euler_RK2_Family.py"
```

4. 批量验证全部默认示例：

```powershell
python ".\期末\上机\03_Python算法模板\verify_all_templates.py"
```

## 对应关系

| 教材算法 | Python 模板 |
|---|---|
| Taylor 级数法 | `01_单步法/T01_Taylor_Series.py` |
| Euler 与二阶单步法 | `01_单步法/T02_Euler_RK2_Family.py` |
| 显式 RK、Heun3、RK3、RK4 | `02_Runge_Kutta/T03_Explicit_Runge_Kutta.py` |
| 隐式 Gauss-RK | `02_Runge_Kutta/T04_Implicit_Gauss_RK.py` |
| 自适应 RKF45 | `02_Runge_Kutta/T05_Adaptive_RKF45.py` |
| Adams-Bashforth | `03_线性多步法/T06_Adams_Bashforth.py` |
| Adams-Moulton 与 PECE | `03_线性多步法/T07_Adams_Moulton_PECE.py` |
| Milne、Hamming、Nyström | `03_线性多步法/T08_Milne_Hamming_Nystrom.py` |
| BDF/Gear | `03_线性多步法/T09_BDF_Gear.py` |
| 方程组与高阶 ODE | `04_方程组与高阶方程/T10_ODE_System_Higher_Order_RK4.py` |
| 隐式方程迭代器 | `04_方程组与高阶方程/T11_Implicit_Equation_Iterators.py` |
| 误差、收敛阶、Richardson、稳定函数 | `05_误差与稳定性/T12_Error_Order_Stability.py` |
| 一般线性多步法阶数与根条件 | `05_误差与稳定性/T13_Linear_Multistep_Order_Root_Check.py` |
