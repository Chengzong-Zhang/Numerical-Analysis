%% =========================================================
%  problem4_fft_compression.m
%  功能：实现改进的 FFT 算法，并用于图像或音频信号的压缩处理
%
%  【数值分析背景：离散 Fourier 变换与 FFT（教材§3.5.2）】
%  离散 Fourier 变换（DFT）定义（教材 p.27 附近）：
%    X[k] = Σ_{n=0}^{N-1} x[n] · e^{-2πi·kn/N},  k=0,...,N-1
%  逆变换（IDFT）：
%    x[n] = (1/N) Σ_{k=0}^{N-1} X[k] · e^{2πi·kn/N}
%
%  DFT 直接计算：O(N^2) 次复数乘法
%  Cooley-Tukey FFT（N=2^p）：O(N log_2 N) 次复数乘法
%  → 当 N=1024，FFT 比 DFT 快约 N/(log_2 N) = 1024/10 ≈ 100 倍！
%
%  FFT 的核心思想：蝴蝶分解（Butterfly Decomposition）
%  将 N 点 DFT 分解为偶数下标和奇数下标的两个 N/2 点 DFT：
%    X[k]     = E[k] + W_N^k · O[k],   k=0,...,N/2-1
%    X[k+N/2] = E[k] - W_N^k · O[k],   k=0,...,N/2-1
%  其中旋转因子 W_N^k = e^{-2πi·k/N}，E[k] 和 O[k] 分别是偶/奇索引的 N/2 点 DFT
%
%  压缩原理：实际信号的能量集中在少数频率分量
%  → 保留幅度最大的 K 个频率系数，将其余置0
%  → 用 IFFT 重建信号，实现有损压缩
%  =========================================================

clear; close all; clc;

fprintf('==============================================\n');
fprintf('  第三章上机作业4：FFT 算法实现与信号压缩\n');
fprintf('==============================================\n\n');

%% =========================================================
%  第一部分：实现并验证 Cooley-Tukey 基-2 FFT 算法
%  =========================================================

fprintf('===== 第一部分：Cooley-Tukey 基-2 FFT 算法验证 =====\n\n');

% ---- 正确性验证：与 MATLAB 内置 fft 结果对比 ----
% 选择 N=64（2^6），包含复数输入以充分测试
N_test = 64;
% randn(1, N)：生成 1×N 的标准正态随机向量（实部）
% 1i：MATLAB 中的虚数单位 i（避免与变量名 i 冲突用 1i）
x_test = randn(1, N_test) + 1i * randn(1, N_test);  % 复数测试信号

% 调用自实现 FFT（局部函数，定义在本文件末）
X_my  = my_fft(x_test);
% 调用 MATLAB 内置 FFT（使用高度优化的 FFTW 库）
% fft(x)：对向量 x 计算离散 Fourier 变换，返回等长复数向量
% 语法：X = fft(x) 或 X = fft(x, N)（第二参数指定 DFT 长度，不足补零）
X_ref = fft(x_test);

% 验证精度：计算两者之差的最大绝对值
% max(abs(v))：向量 v 各元素取模（复数的模），再取最大值
err_verify = max(abs(X_my - X_ref));
fprintf('FFT 算法验证：\n');
fprintf('  测试信号长度 N = %d（= 2^%d）\n', N_test, log2(N_test));
fprintf('  自实现 FFT 与 MATLAB 内置 fft 的最大误差 = %.4e\n', err_verify);
if err_verify < 1e-10
    fprintf('  验证通过：误差在数值精度范围内（< 1e-10）\n\n');
else
    fprintf('  警告：误差偏大，请检查 FFT 实现\n\n');
end

% ---- 运算速度对比（DFT vs FFT）----
fprintf('运算速度对比（DFT vs FFT）：\n');
N_speed = 256;   % 测试规模

% 直接 DFT 实现（O(N^2)）
x_speed = randn(1, N_speed);

% tic/toc：MATLAB 计时函数
% tic：启动秒表计时器
% toc：返回自上次 tic 以来经过的时间（秒）
tic;
X_dft = zeros(1, N_speed);
for k = 0:N_speed-1
    n_vec = 0:N_speed-1;  % 时域下标 n
    % 直接 DFT 定义：X[k] = Σ_n x[n] · e^{-2πi·kn/N}
    X_dft(k+1) = sum(x_speed .* exp(-2i*pi*k*n_vec/N_speed));
    % sum(v)：向量 v 的所有元素之和
    % exp(-2i*pi*...)：计算复数指数（旋转因子），.* 是逐元素乘法
end
time_dft = toc;  % 记录 DFT 耗时

tic;
X_fft = my_fft(x_speed);
time_fft = toc;  % 记录 FFT 耗时

fprintf('  信号长度 N = %d：\n', N_speed);
fprintf('  直接 DFT 耗时 = %.4f 秒（O(N^2) = O(%d)）\n', time_dft, N_speed^2);
fprintf('  自实现 FFT 耗时 = %.6f 秒（O(N·log₂N) = O(%d)）\n', ...
    time_fft, N_speed * log2(N_speed));
fprintf('  加速比 = %.1f 倍\n\n', time_dft / time_fft);

%% =========================================================
%  第二部分：音频信号压缩示例
%  =========================================================

fprintf('===== 第二部分：音频信号的 FFT 频域压缩 =====\n\n');

% ---- 生成合成音频信号（模拟真实音频的多频率叠加特性）----
Fs = 1000;         % 采样频率（Hz）：每秒采样 Fs 个点
% Fs 需满足 Nyquist 采样定理：Fs > 2 × 最高频率（这里最高200Hz，Fs=1000≥400）

% t：时间轴，从 0 到 1 秒，采样间隔 1/Fs
% 语法：start : step : stop（生成等差数列）
t = 0 : 1/Fs : 1 - 1/Fs;   % 长度为 Fs = 1000 的时间向量（秒）
N_sig = length(t);           % 信号点数

% 合成信号：主要由三个正弦波叠加，加少量噪声
% 数值分析含义：模拟真实信号——能量集中在少数频率分量，FFT 后大多数系数接近0
% sin(2*pi*f*t)：频率为 f Hz 的正弦波（2π·f·t 是角频率×时间 = 相位）
x_signal = 2.0 * sin(2*pi*50*t)   ...  % 主频 50 Hz，幅度 2.0
          + 1.0 * sin(2*pi*120*t)  ...  % 120 Hz 分量，幅度 1.0
          + 0.5 * sin(2*pi*200*t)  ...  % 200 Hz 分量，幅度 0.5
          + 0.2 * randn(size(t));        % 加性高斯白噪声（噪声很小）
% randn(size(t))：生成与 t 等大小的标准正态随机向量
% ... 是 MATLAB 续行符（一条语句分多行书写）

fprintf('合成音频信号：\n');
fprintf('  采样频率 Fs = %d Hz\n', Fs);
fprintf('  信号长度 N = %d 点（时长 %.1f 秒）\n', N_sig, N_sig/Fs);
fprintf('  成分：50 Hz（幅度2.0）+ 120 Hz（幅度1.0）+ 200 Hz（幅度0.5）+ 噪声\n\n');

% ---- 补零到最近的2的幂（提高 FFT 效率）----
% nextpow2(n)：返回满足 2^p ≥ n 的最小整数 p
% 语法：p = nextpow2(n)
N_fft = 2^nextpow2(N_sig);
% 补零（Zero-Padding）：在信号末尾填充零到 N_fft 个点
% 数值分析含义：补零不改变信号信息，但提高频率分辨率，并使长度为2的幂（FFT最优）
x_padded = [x_signal, zeros(1, N_fft - N_sig)];
% [A, B]：水平拼接向量 A 和 B；zeros(1, n)：1×n 全零行向量
fprintf('信号补零：%d → %d 点（= 2^%d，满足基-2 FFT 要求）\n\n', ...
    N_sig, N_fft, nextpow2(N_sig));

% ---- 计算 FFT ----
X_freq = my_fft(x_padded);   % 自实现 FFT（复数数组，长度 N_fft）
% X_freq(k+1) 是频率 k·Fs/N_fft Hz 处的复数 Fourier 系数

% ---- 频率轴（单边，只取 0 ~ Fs/2）----
% 对实数信号，FFT 输出具有共轭对称性：X[N-k] = conj(X[k])
% 只需取前半部分（单边谱）：k=0,...,N_fft/2-1，对应频率 0,...,Fs/2 Hz
% 频率分辨率 = Fs / N_fft（Hz/点）
freq_axis = (0 : N_fft/2 - 1) * Fs / N_fft;
% (0:N_fft/2-1)：整数序列 0,1,...,N_fft/2-1
% * Fs/N_fft：换算为 Hz（每个频率间隔 = 采样率/总点数）

% ---- 幅度谱（单边，归一化）----
X_mag = abs(X_freq(1 : N_fft/2));  % 取前半部分的复数模（幅度）
% abs(z)：复数 z 的模（幅度），即 |z| = sqrt(实部^2 + 虚部^2)
% 工程上有时乘以 2/N_fft 归一化，这里只看相对幅度

fprintf('频谱分析（前5个主要频率分量）：\n');
[sorted_mag, sort_idx] = sort(X_mag, 'descend');
% sort(v, 'descend')：对向量 v 降序排列
% 返回两个值：sorted_mag（排序后的值）和 sort_idx（原始下标）
for k = 1:5
    fprintf('  第%d大：f = %6.2f Hz, 幅度 = %.2f\n', k, ...
        freq_axis(sort_idx(k)), sorted_mag(k));
end
fprintf('\n');

%% =========================================================
%  频域截断压缩：保留幅度最大的 K 个系数
%  =========================================================

fprintf('===== 频域截断压缩实验 =====\n\n');

% 测试三种压缩比（保留系数的百分比）
compress_ratios = [0.20, 0.05, 0.01];  % 保留 20%, 5%, 1% 的系数
n_ratios = length(compress_ratios);

% 预分配存储各压缩比的结果
fit_signals = cell(n_ratios, 1);   % cell 数组：存储重建信号
snr_results = zeros(n_ratios, 1);  % 信噪比
rms_results = zeros(n_ratios, 1);  % 重建误差

for r_idx = 1:n_ratios
    ratio = compress_ratios(r_idx);    % 当前保留比例
    K = max(1, round(ratio * N_fft));  % 保留的系数个数（至少1个）
    % round(x)：四舍五入取整；max(a,b)：取两值的较大值

    % ---- 选取幅度最大的 K 个频率系数 ----
    % 对完整双边谱操作（保持实数信号的共轭对称性）
    X_mag_full = abs(X_freq);   % 完整幅度谱（N_fft 个点，双边）
    [~, idx_sorted] = sort(X_mag_full, 'descend');
    % ~ 是占位符（忽略 sort 的第一个返回值，只取下标 idx_sorted）

    % 创建稀疏频谱：只保留最大 K 个系数，其余清零
    X_compressed = zeros(size(X_freq));   % zeros(size(v))：生成与 v 等大的全零数组
    X_compressed(idx_sorted(1:K)) = X_freq(idx_sorted(1:K));
    % 只保留幅度排名前 K 的频率系数，其余为 0

    % ---- 逆 FFT 重建信号 ----
    % my_ifft：利用正 FFT 实现逆 FFT（详见局部函数）
    x_reconstructed = real(my_ifft(X_compressed));
    % real(z)：取复数的实部（消除数值误差引入的微小虚部）
    x_reconstructed = x_reconstructed(1:N_sig);  % 截取到原信号长度

    % ---- 计算重建质量指标 ----
    % 均方根误差（RMS）：衡量重建信号与原信号的差异
    rms_err = sqrt(mean((x_reconstructed - x_signal).^2));

    % 信噪比（Signal-to-Noise Ratio，SNR），单位 dB：
    % SNR = 20·log10(||信号||_2 / ||误差||_2)
    % 越大表示重建质量越好（SNR > 20dB 通常可接受）
    % norm(v)：向量 v 的 2-范数 = sqrt(Σv_i^2)
    % log10(x)：以 10 为底的对数
    snr_db = 20 * log10(norm(x_signal) / (norm(x_reconstructed - x_signal) + 1e-20));
    % + 1e-20：防止分母为零

    rms_results(r_idx) = rms_err;
    snr_results(r_idx) = snr_db;
    fit_signals{r_idx} = x_reconstructed;

    fprintf('保留比例 %5.1f%%（K = %4d 个系数）：RMS = %.4f, SNR = %.1f dB\n', ...
        ratio*100, K, rms_err, snr_db);
end
fprintf('\n');

%% =========================================================
%  第三部分：图像压缩示例（使用二维 FFT）
%  =========================================================

fprintf('===== 第三部分：图像的 FFT 压缩（示例）=====\n\n');

% 生成简单的测试图像（实际使用中可替换为 imread 读入的真实图像）
% 构造 64×64 的合成图像：多个频率的正弦条纹
Nx = 64;   Ny = 64;   % 图像尺寸（必须是2的幂以便 FFT）
[X_img, Y_img] = meshgrid(0:Nx-1, 0:Ny-1);
% meshgrid(x, y)：生成网格坐标矩阵
% X_img(i,j) = x(j)（列方向）；Y_img(i,j) = y(i)（行方向）

% 合成图像：低频和高频正弦条纹的叠加
img_orig = 128 + 60*sin(2*pi*3*X_img/Nx) + 40*cos(2*pi*5*Y_img/Ny) ...
           + 20*sin(2*pi*7*(X_img+Y_img)/Nx) + 10*randn(Ny,Nx);
% 128：灰度中值（使图像在0-255范围内）
img_orig = min(255, max(0, img_orig));  % 裁剪到合法灰度范围 [0, 255]
% min/max：逐元素裁剪（MATLAB 对矩阵操作时逐元素处理）

fprintf('合成测试图像：%d × %d 像素\n', Nx, Ny);

% 二维 FFT（使用 MATLAB 内置 fft2，因为自实现 FFT 按行/列分解即可）
% fft2(img)：二维离散 Fourier 变换，等效于先对每列 FFT，再对每行 FFT
% 数值分析含义：二维 DFT 可分解为多次一维 DFT（可分离性）
% F = fft2(img) 的元素 F(k+1,l+1) 是空间频率 (k,l) 的复数系数
F_img = fft2(img_orig);

% 二维频谱中心化（将零频移到中间，便于可视化）
% fftshift(F)：将 FFT 输出的零频从左上角移到中心
F_img_shift = fftshift(F_img);
% 注意：压缩时对未移位的 F_img 操作，fftshift 只用于可视化

% ---- 图像频域压缩 ----
img_compress_ratios = [0.10, 0.02];   % 保留 10% 和 2% 的频率系数
fprintf('图像 FFT 压缩（保留不同比例频率系数）：\n');

img_reconstructed = cell(length(img_compress_ratios), 1);
img_snr = zeros(length(img_compress_ratios), 1);

for r_idx = 1:length(img_compress_ratios)
    ratio = img_compress_ratios(r_idx);
    K_img = max(1, round(ratio * Nx * Ny));  % 保留的总系数个数
    % Nx * Ny：图像总像素数（= 总频率系数数）

    % 选取幅度最大的 K_img 个二维频率系数
    F_flat = F_img(:);    % (:) 将矩阵展开为列向量
    [~, idx2d] = sort(abs(F_flat), 'descend');   % 按幅度降序排列

    F_compressed = zeros(size(F_img));    % 压缩后的频域（全零初始化）
    F_compressed_flat = F_compressed(:);
    F_compressed_flat(idx2d(1:K_img)) = F_flat(idx2d(1:K_img));
    F_compressed = reshape(F_compressed_flat, Ny, Nx);
    % reshape(v, m, n)：将向量 v 重排为 m×n 矩阵（按列填充）

    % 二维 IFFT 重建图像
    % ifft2(F)：二维逆离散 Fourier 变换
    img_rec = real(ifft2(F_compressed));
    img_rec = min(255, max(0, img_rec));   % 裁剪灰度范围

    % 图像 SNR
    img_snr(r_idx) = 20 * log10(norm(img_orig(:)) / ...
        (norm(img_rec(:) - img_orig(:)) + 1e-20));
    img_reconstructed{r_idx} = img_rec;

    fprintf('  保留 %5.1f%%（%d个系数）：SNR = %.1f dB\n', ...
        ratio*100, K_img, img_snr(r_idx));
end
fprintf('\n');

%% =========================================================
%  绘图：信号压缩效果
%  =========================================================

figure('Name', 'FFT 音频信号压缩效果', 'NumberTitle', 'off', ...
    'Position', [50, 50, 1400, 900]);

% 第1行：原始信号时域和频域
subplot(n_ratios+1, 2, 1);
% plot(t, x)：绘制时域信号
plot(t, x_signal, 'b-', 'LineWidth', 0.8);
xlabel('时间 (s)');  ylabel('幅值');
title(sprintf('原始信号（时域）：50+120+200 Hz 正弦叠加+噪声'), 'FontSize', 9);
xlim([0, 0.1]);   % 只显示前0.1秒，更清晰
grid on;

subplot(n_ratios+1, 2, 2);
% stem(x, y)：用垂直线画离散点图（适合展示离散频谱）
stem(freq_axis, X_mag, 'b', 'Marker', 'none', 'LineWidth', 0.5);
xlabel('频率 (Hz)');  ylabel('幅度');
title('原始信号幅度谱（单边 0~500 Hz）', 'FontSize', 9);
xlim([0, 300]);   % 只显示 0~300 Hz 范围
grid on;

% 第2~4行：各压缩比的重建信号
for r_idx = 1:n_ratios
    ratio = compress_ratios(r_idx);
    K = max(1, round(ratio * N_fft));

    subplot(n_ratios+1, 2, 2*r_idx+1);
    plot(t, x_signal, 'k-', 'LineWidth', 0.8, 'DisplayName', '原始信号');
    hold on;
    plot(t, fit_signals{r_idx}, 'r-', 'LineWidth', 1.2, 'DisplayName', '重建信号');
    xlabel('时间 (s)');  ylabel('幅值');
    title(sprintf('保留 %.0f%% 系数（K=%d），SNR=%.1f dB', ...
        ratio*100, K, snr_results(r_idx)), 'FontSize', 9);
    legend('Location', 'northeast', 'FontSize', 7);
    xlim([0, 0.1]);
    grid on;  hold off;

    subplot(n_ratios+1, 2, 2*r_idx+2);
    % 计算压缩后的幅度谱（重建后的频谱）
    X_rec_fft = my_fft([fit_signals{r_idx}, zeros(1, N_fft - N_sig)]);
    X_rec_mag  = abs(X_rec_fft(1:N_fft/2));
    stem(freq_axis, X_rec_mag, 'r', 'Marker', 'none', 'LineWidth', 0.5);
    xlabel('频率 (Hz)');  ylabel('幅度');
    title(sprintf('压缩后幅度谱（%.0f%% 系数）', ratio*100), 'FontSize', 9);
    xlim([0, 300]);
    grid on;
end

sgtitle('FFT 频域截断压缩：保留不同比例频率系数的重建效果对比', 'FontSize', 12);

% ---- 图像压缩效果图 ----
figure('Name', 'FFT 图像压缩效果', 'NumberTitle', 'off', ...
    'Position', [100, 100, 1000, 600]);

subplot(2, 3, 1);
% imagesc(img)：将矩阵 img 显示为图像，自动拉伸颜色范围
imagesc(img_orig);
% colormap gray：使用灰度色表
colormap(gray);
% colorbar：显示颜色条（像素值刻度）
colorbar;
title('原始图像', 'FontSize', 10);
axis equal tight;  % equal：等比例显示；tight：坐标轴贴合图像边缘

subplot(2, 3, 2);
% fftshift 后的二维频谱对数幅度（取 log 便于显示大动态范围）
imagesc(log(1 + abs(F_img_shift)));
colormap(gray);  colorbar;
title('二维 FFT 幅度谱（log 刻度）', 'FontSize', 10);
% log(1+x)：避免 log(0)=-Inf，加1保证自变量>0
axis equal tight;

for r_idx = 1:length(img_compress_ratios)
    subplot(2, 3, r_idx + 2);
    imagesc(img_reconstructed{r_idx});
    colormap(gray);  colorbar;
    title(sprintf('保留 %.0f%% 系数\nSNR = %.1f dB', ...
        img_compress_ratios(r_idx)*100, img_snr(r_idx)), 'FontSize', 9);
    axis equal tight;
end

% 最后一个子图：SNR vs 压缩比
subplot(2, 3, 5);
compress_pct = img_compress_ratios * 100;
% plot(x, y, '-o')：带圆圈标记的折线图
plot(compress_pct, img_snr, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, ...
    'MarkerFaceColor', 'b');
xlabel('保留系数比例（%）');
ylabel('图像 SNR（dB）');
title('图像质量 vs 压缩比', 'FontSize', 10);
grid on;

sgtitle('二维 FFT 图像压缩效果', 'FontSize', 12);

fprintf('图形已生成，请查看弹出的图形窗口。\n\n');
fprintf('程序运行完毕。\n');


%% =========================================================
%  局部函数（Local Functions）—— 必须置于脚本末尾
%  =========================================================

function X = my_fft(x)
    % my_fft：Cooley-Tukey 基-2 DIT（时域抽取）FFT 递归实现
    %
    % 【数值分析：FFT 算法（教材§3.5.2）】
    % 离散 Fourier 变换（DFT）的定义：
    %   X[k] = Σ_{n=0}^{N-1} x[n] · W_N^{kn},  W_N = e^{-2πi/N}（旋转因子）
    %
    % Cooley-Tukey 基-2 DIT 分解（N=2^p）：
    % 按时域奇偶下标分组：
    %   E[k] = DFT({x[0], x[2], ..., x[N-2]})  （偶下标，N/2 点 DFT）
    %   O[k] = DFT({x[1], x[3], ..., x[N-1]})  （奇下标，N/2 点 DFT）
    % 蝴蝶运算（Butterfly Operation）：
    %   X[k]       = E[k] + W_N^k · O[k],  k=0,...,N/2-1  （上半部分）
    %   X[k + N/2] = E[k] - W_N^k · O[k],  k=0,...,N/2-1  （下半部分）
    % 递归分解：每层将 N 点 DFT 变为两个 N/2 点 DFT，共 log_2(N) 层
    % 计算量：O(N log_2 N) 次复数乘法（比 DFT 的 O(N^2) 快 N/log_2(N) 倍）
    %
    % 输入：x - 输入信号向量（长度 N 必须是 2 的幂，支持复数输入）
    % 输出：X - DFT 结果向量（与 x 等长的复数向量）

    N = length(x);   % 信号长度

    % 递归终止条件：N=1 时，1点 DFT = 原值本身（恒等变换）
    if N == 1
        X = x;
        return;   % return：立即结束函数，返回结果
    end

    % 检查 N 是否为 2 的幂（基-2 FFT 的必要条件）
    if mod(N, 2) ~= 0
        % mod(a, b)：a 除以 b 的余数；~= 是"不等于"逻辑算符
        error('my_fft：输入长度 N=%d 不是2的幂，请预先补零。', N);
        % error('msg')：抛出错误并终止程序，显示错误信息
    end

    % ---- 时域奇偶分解 ----
    x_even = x(1:2:end);   % 取偶数下标元素：x[0], x[2], ..., x[N-2]
    % 1:2:end 是步长为2的整数序列，从1到最后（MATLAB 下标从1开始，对应数学下标 0,2,...）
    x_odd  = x(2:2:end);   % 取奇数下标元素：x[1], x[3], ..., x[N-1]
    % 2:2:end 步长为2，从2开始（对应数学下标 1,3,...）

    % ---- 递归计算 N/2 点 FFT ----
    E = my_fft(x_even);    % 偶数部分的 N/2 点 DFT：E[k] = DFT({x[0],x[2],...})
    O = my_fft(x_odd);     % 奇数部分的 N/2 点 DFT：O[k] = DFT({x[1],x[3],...})
    % 递归调用：每次规模减半，共 log_2(N) 层递归

    % ---- 计算旋转因子向量 W_N^k，k=0,...,N/2-1 ----
    % W_N^k = e^{-2πi·k/N}（复数指数，即单位圆上的旋转）
    k = 0 : N/2 - 1;         % k 的取值范围，生成行向量 [0, 1, 2, ..., N/2-1]
    W = exp(-2i * pi * k / N);  % 旋转因子向量（行向量，每个元素是复数）
    % exp(z)：复数指数函数（当 z=iθ 时，exp(iθ) = cos(θ)+i·sin(θ)，Euler 公式）
    % -2i：乘以 -2i 确保 DFT 定义中的负号（-2πi·kn/N 中的负号来自于此）
    % pi：MATLAB 内置常量 π ≈ 3.14159...

    % ---- 蝴蝶运算：合并两个 N/2 点 DFT 为一个 N 点 DFT ----
    WO = W .* O;   % W_N^k · O[k]：旋转因子与奇数 DFT 的逐元素乘积
    % .* 是逐元素乘法（W 和 O 都是长度 N/2 的行向量）

    % 上半部分 X[k] = E[k] + W_N^k·O[k]，k=0,...,N/2-1
    % 下半部分 X[k+N/2] = E[k] - W_N^k·O[k]，k=0,...,N/2-1
    X = [E + WO, E - WO];   % 水平拼接上下两半，得到长度 N 的输出
    % [A, B]：水平拼接行向量 A 和 B
end


function x = my_ifft(X)
    % my_ifft：逆 FFT（IFFT），利用正 FFT 实现（共轭对称技巧）
    %
    % 【数值分析：IFFT 的对称性推导】
    % IDFT 定义：x[n] = (1/N) Σ_{k=0}^{N-1} X[k] · e^{+2πi·kn/N}
    % 与 DFT：   X[k] = Σ_{n=0}^{N-1} x[n] · e^{-2πi·kn/N}
    % 比较可见：IDFT 中的指数符号与 DFT 相反（+2πi 而非 -2πi）
    %
    % 关键技巧：IFFT 可由 FFT 实现，无需单独编写逆变换：
    %   IFFT(X) = (1/N) · conj(FFT(conj(X)))
    % 证明：设 Y = conj(X)（逐元素取共轭），则
    %   FFT(Y)[k] = Σ_n conj(X[n]) · e^{-2πi·kn/N}
    %             = conj(Σ_n X[n] · e^{+2πi·kn/N})
    %             = conj(N · IDFT(X)[k])
    %   ∴ IDFT(X)[k] = (1/N) · conj(FFT(conj(X))[k])
    %
    % 输入：X - 频域信号向量（长度 N 的复数向量）
    % 输出：x - IDFT 重建的时域信号（复数向量，实信号取 real(x) 即可）

    N = length(X);
    % conj(v)：向量 v 的逐元素共轭（虚部变号）
    % 步骤：① conj(X) → ② my_fft(...) → ③ conj(...) → ④ /N
    x = conj(my_fft(conj(X))) / N;
    % 除以 N：IDFT 的归一化因子（DFT 和 IDFT 之间的 1/N 约定）
end
