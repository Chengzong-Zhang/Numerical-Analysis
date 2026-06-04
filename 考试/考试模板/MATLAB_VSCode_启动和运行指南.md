# MATLAB 与 VS Code 启动运行指南

本文档适用于当前项目：

```text
C:\coding\matlab
```

核心原则只有一个：先看清楚你现在在哪个窗口输入命令。

## 1. 先分清两个命令窗口

### PowerShell / VS Code 终端

如果你看到类似下面的提示符：

```powershell
(base) PS C:\coding\matlab>
```

说明你在 PowerShell 里。这里可以输入系统命令，例如：

```powershell
matlab -batch "disp('ok')"
```

### MATLAB 命令窗口

如果你看到类似下面的提示符：

```matlab
>>
```

说明你已经在 MATLAB 里面。这里不能输入 `matlab -batch ...`，只能输入 MATLAB 语句，例如：

```matlab
disp('ok')
```

不要在 `>>` 后面输入：

```matlab
matlab -batch "..."
```

这样会报错，因为 `matlab -batch` 是 PowerShell 命令，不是 MATLAB 语句。

## 2. 检查 MATLAB 能不能从命令行启动

打开 VS Code 终端或 PowerShell，确认提示符类似：

```powershell
(base) PS C:\coding\matlab>
```

输入：

```powershell
matlab -batch "disp('ok')"
```

正常结果应该是：

```text
ok
```

如果能输出 `ok`，说明 MATLAB 命令行启动正常。

## 3. 在 PowerShell 里运行某个 MATLAB 文件

假设要运行这个文件：

```text
C:\coding\matlab\数值积分\4-上机\problem1_pi_integration.m
```

在 PowerShell 或 VS Code 终端里输入：

```powershell
matlab -batch "run('C:\coding\matlab\数值积分\4-上机\problem1_pi_integration.m')"
```

注意：

- 这条命令必须在 `(base) PS ...>` 这种 PowerShell 提示符下输入。
- 不要在 MATLAB 的 `>>` 里输入这条命令。
- 路径要写真实文件路径，不能写 `...\你的文件.m` 这种占位符。

## 4. 在 MATLAB 图形界面里运行文件

先单独打开 MATLAB 桌面版。

### 方法一：直接 run 完整路径

在 MATLAB 命令窗口里，也就是 `>>` 后面输入：

```matlab
run('C:\coding\matlab\数值积分\4-上机\problem1_pi_integration.m')
```

### 方法二：先切换目录，再运行脚本名

在 MATLAB 命令窗口里输入：

```matlab
cd('C:\coding\matlab\数值积分\4-上机')
problem1_pi_integration
```

这里第二行不用写 `.m`。

如果要运行另一个文件，例如：

```text
C:\coding\matlab\逼近\3-上机\problem2_least_squares_fit.m
```

就在 MATLAB 里输入：

```matlab
cd('C:\coding\matlab\逼近\3-上机')
problem2_least_squares_fit
```

## 5. 在 VS Code 里运行当前 MATLAB 文件

当前项目已经配置了 VS Code 任务。

步骤：

1. 用 VS Code 打开一个 `.m` 文件。
2. 按 `Ctrl + Shift + B`。
3. 选择 `Run current MATLAB file`。

这个任务会自动执行类似下面的命令：

```powershell
matlab -batch "run('当前打开的完整文件路径')"
```

所以它不依赖 MATLAB 当前目录，也不需要你手动 `cd`。

## 6. 如果你点右上角三角形运行按钮

右上角三角形按钮可能来自不同扩展。

### Code Runner 的三角形按钮

当前项目已经在 `.vscode/settings.json` 里配置了 Code Runner：

```json
"matlab": "matlab -batch \"run('$fullFileName')\""
```

所以如果按钮来自 Code Runner，它应该会用完整路径运行当前文件。

### MathWorks MATLAB 扩展的运行按钮

如果按钮来自 MathWorks MATLAB 扩展，它可能会连接 MATLAB 会话，并在 MATLAB 当前目录中按文件名运行。

如果出现：

```text
File is not found in the current folder or on the MATLAB path
```

说明 MATLAB 当前目录不在脚本所在文件夹。解决方法：

```matlab
cd('脚本所在文件夹')
脚本名
```

例如：

```matlab
cd('C:\coding\matlab\数值积分\4-上机')
problem1_pi_integration
```

或者改用 `Ctrl + Shift + B` 的 `Run current MATLAB file`。

## 7. 常见错误与解决办法

### 错误一：不支持将脚本 matlab 作为函数执行

常见原因：你在 MATLAB 的 `>>` 里输入了：

```matlab
matlab -batch "..."
```

解决办法：去 PowerShell 或 VS Code 终端里输入这条命令。

### 错误二：File is not found in the current folder or on the MATLAB path

常见原因：MATLAB 当前目录不是 `.m` 文件所在目录。

解决办法一，在 MATLAB 里使用完整路径：

```matlab
run('C:\coding\matlab\数值积分\4-上机\problem1_pi_integration.m')
```

解决办法二，先 `cd` 到脚本目录：

```matlab
cd('C:\coding\matlab\数值积分\4-上机')
problem1_pi_integration
```

### 错误三：输入了 `...\你的文件.m` 后找不到文件

`...\你的文件.m` 只是示例占位符，不能原样复制。

要换成真实路径，例如：

```matlab
run('C:\coding\matlab\数值积分\4-上机\problem1_pi_integration.m')
```

### 错误四：路径里有中文会不会出问题

MATLAB R2024a 通常可以处理中文路径。建议始终把路径放在单引号里：

```matlab
run('C:\coding\matlab\数值积分\4-上机\problem1_pi_integration.m')
```

## 8. 推荐使用顺序

平时最推荐：

```text
VS Code 打开 .m 文件 -> Ctrl + Shift + B -> Run current MATLAB file
```

如果 VS Code 按钮不稳定，就用 PowerShell：

```powershell
matlab -batch "run('完整文件路径.m')"
```

如果已经打开 MATLAB 图形界面，就用 MATLAB 命令窗口：

```matlab
run('完整文件路径.m')
```

## 9. 当前项目常用运行示例

运行数值积分第一题：

```powershell
matlab -batch "run('C:\coding\matlab\数值积分\4-上机\problem1_pi_integration.m')"
```

在 MATLAB 图形界面中运行：

```matlab
run('C:\coding\matlab\数值积分\4-上机\problem1_pi_integration.m')
```

运行数值积分总脚本：

```powershell
matlab -batch "run('C:\coding\matlab\数值积分\4-上机\run_and_save.m')"
```

在 MATLAB 图形界面中运行：

```matlab
run('C:\coding\matlab\数值积分\4-上机\run_and_save.m')
```

## 10. 一句话记忆

PowerShell 里输入：

```powershell
matlab -batch "run('完整路径.m')"
```

MATLAB 的 `>>` 里输入：

```matlab
run('完整路径.m')
```

