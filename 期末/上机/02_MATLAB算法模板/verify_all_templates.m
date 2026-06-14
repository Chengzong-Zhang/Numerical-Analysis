function verify_all_templates()                            % 函数文件名与函数名一致，便于在 MATLAB 命令窗直接调用。
root = fileparts(mfilename('fullpath'));                   % mfilename('fullpath') 获取本验证器的完整路径。
items = dir(fullfile(root,'**','*.m'));                    % ** 表示递归查找全部子目录中的 MATLAB 文件。
items = items(~strcmp({items.name},'verify_all_templates.m')); % 排除验证器自身，避免递归运行自己。
failures = strings(0,1);                                   % 创建空 string 列向量，用于收集失败文件。
set(groot,'defaultFigureVisible','off');                   % 验证时隐藏图窗，避免批量弹出窗口。
for i = 1:numel(items)                                     % 逐个验证模板。
    file = fullfile(items(i).folder,items(i).name);         % fullfile 跨平台拼接完整路径。
    fprintf('RUN %s\n',file);                              % 输出当前正在验证的文件。
    try                                                    % try/catch 让单个模板失败后仍继续检查其他模板。
        run_one(file);                                     % 在独立函数工作区运行，模板中的 clear 不会清掉主循环变量。
        close all force;                                   % 关闭该模板产生的图窗。
    catch ME                                               % 捕获 MATLAB 异常对象。
        failures(end+1,1)=string(file)+" -> "+string(ME.message); % 保存文件名和错误消息。
    end                                                    % 异常处理结束。
end                                                        % 批量验证结束。
fprintf('TOTAL=%d, FAILURES=%d\n',numel(items),numel(failures)); % 输出汇总。
disp(failures);                                            % 显示全部失败项；空数组表示全部成功。
assert(isempty(failures),'存在运行失败的模板。');           % 有失败项时令自动验证返回非零状态。
end                                                        % 主验证函数结束。

function run_one(file)                                     % 每个模板在这个独立工作区中运行。
run(file);                                                 % run 执行脚本；脚本中的 clear 只影响当前辅助函数工作区。
end                                                        % 辅助函数结束。
