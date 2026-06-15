# 批量验证全部 Python ODE 模板。  # 每个模板在独立子进程中运行，避免变量和图窗相互影响。
import os  # 标准库 os 用于复制并修改子进程环境变量。
import subprocess  # 标准库 subprocess 用于独立运行每一个模板。
import sys  # 标准库 sys 提供当前 Python 解释器完整路径。
from pathlib import Path  # pathlib 提供面向对象的跨平台路径操作。

root = Path(__file__).resolve().parent  # __file__ 是当前脚本路径；resolve 得到绝对路径。
files = sorted(path for path in root.rglob("*.py") if path.name != Path(__file__).name)  # rglob 递归查找，并排除验证器自身。
environment = os.environ.copy()  # 复制当前环境变量，避免破坏父进程环境。
environment["MPLBACKEND"] = "Agg"  # 使用无窗口绘图后端，使 plt.show() 在批量验证时不弹窗。
failures = []  # 创建空列表保存运行失败的模板。
for path in files:  # 逐个验证模板。
    print(f"RUN {path}")  # 输出当前正在验证的文件路径。
    result = subprocess.run([sys.executable, str(path)], env=environment, capture_output=True, text=True, encoding="utf-8", errors="replace")  # 使用同一 Python 解释器运行；errors=replace 兼容 Windows 子进程混合编码输出。
    if result.returncode != 0:  # 子进程返回码非零表示模板运行失败。
        failures.append((path, result.stderr))  # 保存失败路径与错误信息。
print(f"TOTAL={len(files)}, FAILURES={len(failures)}")  # 输出批量验证汇总。
for path, error in failures:  # 遍历显示全部失败项。
    print(f"FAIL {path}\n{error}")  # 打印失败文件与完整错误信息。
assert not failures, "存在运行失败的 Python 模板。"  # 有失败项时使验证器返回失败状态。
