# Ollama + DeepSeek R1 本地启动指南

本文档用于在 Windows + VS Code 环境下启动和使用本地 Ollama 的 `deepseek-r1:7b` 模型。适合断网时使用，前提是模型已经提前下载到本机。

## 一、确认模型是否已经在本地

打开 PowerShell，运行：

```powershell
ollama list
```

如果能看到类似下面内容，说明模型已经在本地：

```text
NAME              ID              SIZE
deepseek-r1:7b    755ced02ce7b    4.7 GB
```

如果没有这个模型，并且当前有网络，可以下载：

```powershell
ollama pull deepseek-r1:7b
```

断网时不能下载新模型，只能使用已经存在于本机的模型。

## 二、命令行方式启动 DeepSeek

最简单的方式是直接运行：

```powershell
ollama run deepseek-r1:7b
```

成功后会进入对话模式，可以直接输入问题，例如：

```text
你好，只回复 OK
```

退出对话可以按：

```text
Ctrl + D
```

或者直接关闭当前终端。

## 三、启动 Ollama 本地服务

很多 VS Code 插件不是直接调用 `ollama run`，而是连接 Ollama 本地服务：

```text
http://127.0.0.1:11434
```

通常 Ollama 会自动在后台启动服务。你可以用下面命令测试服务是否可用：

```powershell
Invoke-RestMethod http://127.0.0.1:11434/api/tags
```

如果能返回模型列表，说明服务正常。

也可以手动启动服务：

```powershell
ollama serve
```

如果出现：

```text
Error: listen tcp 127.0.0.1:11434: bind: Only one usage of each socket address...
```

这通常不是坏事，而是说明 `11434` 端口已经有 Ollama 服务在运行。此时直接使用即可。

## 四、图形化窗口方式启动 Ollama

如果你安装的是 Windows 版 Ollama，一般可以通过图形化方式启动：

1. 打开开始菜单
2. 搜索 `Ollama`
3. 点击启动
4. 启动后 Ollama 会在后台运行

然后打开 PowerShell 测试：

```powershell
ollama list
ollama run deepseek-r1:7b
```

如果这两个命令正常，说明图形化启动成功。

## 五、在 VS Code 中使用本地 DeepSeek

推荐使用已经安装的 `Ollama Local Chat` 插件。

### 1. 确认插件已安装

```powershell
code --list-extensions | Select-String -Pattern "ollama" -CaseSensitive:$false
```

如果看到：

```text
maurokrekels.ollama-chat-vscode
```

说明插件已经安装。

如果没安装，联网时运行：

```powershell
code --install-extension MauroKrekels.ollama-chat-vscode
```

断网环境不能从 Marketplace 安装插件，需要提前准备 `.vsix` 文件。

### 2. 确认 VS Code 配置

配置文件位置：

```text
C:\Users\你的用户名\AppData\Roaming\Code\User\settings.json
```

其中应该包含：

```json
"ollama.baseUrl": "http://127.0.0.1:11434",
"ollama.model": "deepseek-r1:7b"
```

你的机器上当前配置应为：

```json
"ollama.baseUrl": "http://127.0.0.1:11434",
"ollama.model": "deepseek-r1:7b"
```

### 3. 打开插件聊天窗口

在 VS Code 中按：

```text
Ctrl + Shift + P
```

搜索并运行：

```text
Open Ollama Chat
```

如果没有这个命令，就搜索：

```text
Ollama
```

打开聊天窗口后，输入问题即可使用本地 DeepSeek。

### 4. 用命令 ID 启动 Ollama Chat

当前安装的插件是：

```text
maurokrekels.ollama-chat-vscode
```

它注册的 VS Code 命令 ID 是：

```text
ollama.openChat
```

图形化启动方法：

1. 在 VS Code 中按 `Ctrl + Shift + P`
2. 输入 `Open Ollama Chat`
3. 回车

如果中文界面搜不到英文命令，可以直接搜索：

```text
Ollama
```

也可以给它绑定快捷键。打开：

```text
C:\Users\你的用户名\AppData\Roaming\Code\User\keybindings.json
```

加入：

```json
[
  {
    "key": "ctrl+alt+o",
    "command": "ollama.openChat"
  }
]
```

之后在 VS Code 中按：

```text
Ctrl + Alt + O
```

即可打开 Ollama Chat。

### 5. 当前插件是否有侧边栏

当前安装的 `Ollama Local Chat` 插件没有注册左侧侧边栏入口。它只提供：

- 命令：`Open Ollama Chat`
- 配置项：`ollama.baseUrl`
- 配置项：`ollama.model`

也就是说，它更像是一个通过命令打开的聊天窗口，而不是固定在 VS Code 左侧活动栏里的侧边栏插件。

如果想要真正显示为 VS Code 侧边栏，可以改用其他支持 Ollama 的插件，例如：

```powershell
code --install-extension ollacoder.ollacoder
```

`OllaCoder` 支持通过 Ollama 使用本地模型，并提供侧边栏入口。安装后常见启动方式是：

```text
Ctrl + Shift + L
```

或者在 VS Code 左侧活动栏点击 `OllaCoder` 图标。

常用配置项：

```json
"ollacoder.ollama.baseUrl": "http://localhost:11434",
"ollacoder.models.chat": "deepseek-r1:7b"
```

注意：`deepseek-r1:7b` 更适合聊天和推理；如果要在侧边栏里使用 Agent、自动改代码、补全代码，更推荐额外提前下载代码模型，例如：

```powershell
ollama pull qwen2.5-coder:7b
```

断网使用前，插件和模型都必须提前安装/下载好。

## 六、常见问题处理

### 1. Error: upgrade in progress

这是 Ollama 正在升级或升级状态卡住。

先结束 Ollama 进程：

```powershell
taskkill /F /IM ollama.exe /T
```

再删除可能残留的升级锁：

```powershell
Remove-Item "$env:USERPROFILE\.ollama\update.lock" -ErrorAction SilentlyContinue
Remove-Item "$env:USERPROFILE\.ollama\.ollama-upgrade" -ErrorAction SilentlyContinue
Remove-Item "$env:TEMP\ollama*" -Recurse -ErrorAction SilentlyContinue
```

然后重新测试：

```powershell
ollama list
ollama run deepseek-r1:7b
```

### 2. ollama serve 提示端口被占用

错误示例：

```text
listen tcp 127.0.0.1:11434: bind: Only one usage of each socket address
```

这表示 Ollama 服务已经在运行。直接测试：

```powershell
ollama list
Invoke-RestMethod http://127.0.0.1:11434/api/tags
```

如果能返回模型列表，就不用管这个报错。

### 3. 代理导致 Ollama 报错

如果日志里出现类似：

```text
proxyconnect tcp: dial tcp 127.0.0.1:7890: connectex
```

说明当前终端设置了代理，但代理软件没有运行。

临时清除当前 PowerShell 的代理：

```powershell
Remove-Item Env:HTTP_PROXY -ErrorAction SilentlyContinue
Remove-Item Env:HTTPS_PROXY -ErrorAction SilentlyContinue
Remove-Item Env:http_proxy -ErrorAction SilentlyContinue
Remove-Item Env:https_proxy -ErrorAction SilentlyContinue
```

本地使用已经下载好的模型时，一般不需要联网，也不需要代理。

### 4. VS Code 插件连不上 Ollama

先确认 Ollama 服务可用：

```powershell
Invoke-RestMethod http://127.0.0.1:11434/api/tags
```

如果不通，先启动 Ollama：

```powershell
ollama serve
```

如果提示端口占用，说明服务已启动。

再确认 VS Code 设置：

```json
"ollama.baseUrl": "http://127.0.0.1:11434",
"ollama.model": "deepseek-r1:7b"
```

修改设置后重启 VS Code：

```powershell
code --reuse-window .
```

## 七、断网使用 checklist

断网前确认：

- Ollama 已安装
- `deepseek-r1:7b` 已经通过 `ollama list` 显示在本地
- VS Code 插件已经安装
- VS Code `settings.json` 已配置 `ollama.baseUrl` 和 `ollama.model`

断网后启动顺序：

```powershell
ollama list
ollama run deepseek-r1:7b
```

如果命令行能对话，再打开 VS Code：

```powershell
code --reuse-window .
```

然后运行：

```text
Open Ollama Chat
```

## 八、推荐使用方式

如果只是问问题、解释代码、复习知识：

```powershell
ollama run deepseek-r1:7b
```

如果想在 VS Code 里边看代码边聊天：

```text
Open Ollama Chat
```

如果要让 AI 自动改文件，`deepseek-r1:7b` 不一定稳定。更推荐额外准备代码模型，例如：

```powershell
ollama pull qwen2.5-coder:7b
```

但断网前必须提前下载。
