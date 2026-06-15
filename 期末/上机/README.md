# 第五章：数值 ODE 上机文件索引

本目录依据 `数值ODE/第五章.pdf` 整理。旧文件仅移动分类，未修改内容；新增 MATLAB 与 Python 模板均可直接运行。

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
- `03_Python算法模板`：与 MATLAB 模板一一对应的逐行详注 Python 模板。

## 本地大模型命令行使用指南

本机通过 **Ollama** 运行本地大模型。以下命令均在 PowerShell 中执行。

### 1. 打开命令行

任选一种方式：

1. 按 `Win + X`，选择“终端”或“终端（管理员）”。
2. 在 VS Code 中按 `` Ctrl + ` `` 打开集成终端。
3. 在当前文件夹地址栏输入 `powershell` 后按 Enter。

先检查 Ollama 和本地模型：

```powershell
ollama --version
ollama list
```

截至 2026-06-15，本机已安装并验证可运行的模型为：

```text
qwen35-9b-uncensored:latest
deepseek-r1:7b
```

### 2. 打开 Qwen

进入持续对话模式：

```powershell
ollama run qwen35-9b-uncensored:latest
```

出现 `>>>` 后直接输入问题并按 Enter。结束对话时输入：

```text
/bye
```

只提问一次并在回答后自动返回 PowerShell：

```powershell
ollama run qwen35-9b-uncensored:latest "请用 MATLAB 编写经典四阶 Runge-Kutta 方法"
```

### 3. 打开 DeepSeek

进入持续对话模式：

```powershell
ollama run deepseek-r1:7b
```

只提问一次：

```powershell
ollama run deepseek-r1:7b "解释 Adams-Bashforth 方法的基本思想"
```

不显示模型的思考过程：

```powershell
ollama run deepseek-r1:7b "解释 Adams-Bashforth 方法的基本思想" --hidethinking
```

### 4. 关于 gamma / gemma

当前 `ollama list` 中没有名为 `gamma` 的本地模型，因此不能直接执行 `ollama run gamma`。

如果所指的是 Google 的 **Gemma** 模型，可先安装一个 Gemma 模型，再运行它：

```powershell
ollama pull gemma3:4b
ollama run gemma3:4b
```

下载会占用网络流量和磁盘空间。安装后应以 `ollama list` 显示的完整名称为准。如果 `gamma` 是以后创建的自定义模型，同样先用 `ollama list` 查到它的准确名称，再执行：

```powershell
ollama run <模型完整名称>
```

### 5. 常用管理命令

```powershell
# 查看已经下载的模型
ollama list

# 查看当前载入内存的模型
ollama ps

# 停止指定模型，释放内存/显存
ollama stop qwen35-9b-uncensored:latest
ollama stop deepseek-r1:7b

# 查看模型信息
ollama show qwen35-9b-uncensored:latest
ollama show deepseek-r1:7b
```

模型首次启动通常较慢，因为需要载入内存或显存。若同时运行多个模型导致电脑变慢，先用 `ollama ps` 查看，再用 `ollama stop <模型名>` 停止不用的模型。

### 6. 常见故障

#### 提示 `ollama` 不是命令

关闭并重新打开 PowerShell。若仍失败，直接运行：

```powershell
& "$env:LOCALAPPDATA\Programs\Ollama\ollama.exe" list
```

#### 提示无法连接 Ollama

先从 Windows 开始菜单启动 Ollama；也可以在一个单独的 PowerShell 窗口执行以下命令，并保持该窗口开启：

```powershell
ollama serve
```

然后在另一个 PowerShell 窗口中执行 `ollama list` 或 `ollama run ...`。

#### 提示找不到模型

模型名称必须与 `ollama list` 完全一致，包括冒号后的标签。例如本机 Qwen 的正确名称不是 `qwen`，而是：

```powershell
ollama run qwen35-9b-uncensored:latest
```
