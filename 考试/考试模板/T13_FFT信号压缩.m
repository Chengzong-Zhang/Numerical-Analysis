%% =========================================================
%  模板T13：FFT（快速傅里叶变换）与信号/图像压缩
%
%  【DFT 定义】
%  对长度为 N 的序列 x = [x₀, x₁, ..., x_{N-1}]，
%  离散傅里叶变换（DFT）为：
%    X_k = Σ_{n=0}^{N-1} x_n · e^{-i2πkn/N}，k=0,...,N-1
%
%  【FFT 优势】
%  DFT 直接计算需要 O(N²) 次运算
%  Cooley-Tukey 基-2 FFT 只需 O(N·log₂N) 次运算（N 必须是2的幂次）
%
%  【频域压缩思路】
%  1. 对信号做 FFT（变换到频域）
%  2. 只保留幅值最大的若干个频率系数（其余置零）
%  3. 对截断后的频谱做 IFFT（反变换回时域）
%  4. 保留系数越少 = 压缩比越高，但失真越大
%
%  【信噪比 SNR（Signal-to-Noise Ratio）】
%  SNR = 20·log₁₀(||原始信号||₂ / ||误差||₂)（单位：dB）
%  SNR 越高越好（信号远强于失真噪声）
%  =========================================================

clear; close all; clc;

%% ===================================================
%  第一部分：FFT 算法验证与 DFT vs FFT 速度对比
%  ===================================================
fprintf('===== 第一部分：FFT 验证与速度对比 =====\n\n');

%% ===== ① 修改这里：选择信号长度（必须是2的幂次！）=====
N_test = 1024;   % 信号长度，必须是 2^k（如 256, 512, 1024, 2048）
%% =====================================================

% 生成测试信号（复数随机信号，用于验证）
x_test = rand(1, N_test) + 1i * rand(1, N_test);
% rand(1,N)：N 个 [0,1] 均匀随机数；1i 是虚数单位

% 方法1：MATLAB 内置 fft（高度优化的 FFT）
X_matlab = fft(x_test);

% 方法2：自实现 FFT（Cooley-Tukey 基-2 DIT 递归算法）
X_my = my_fft(x_test);

% 验证两者一致
err_verify = max(abs(X_matlab - X_my));
fprintf('自实现 FFT 与 MATLAB fft 的最大误差：%.4e\n', err_verify);
if err_verify < 1e-10
    fprintf('验证通过！两者结果一致。\n\n');
end

% DFT vs FFT 速度对比
N_sizes = [64, 128, 256, 512, 1024];   % 测试不同长度
fprintf('%-10s  %-15s  %-15s  %-10s\n', 'N', 'DFT 运算量', 'FFT 运算量', '速度比');
fprintf('%s\n', repmat('-', 1, 52));
for N_k = N_sizes
    ops_dft = N_k^2;
    ops_fft = N_k * log2(N_k);
    fprintf('%-10d  %-15.0f  %-15.0f  %-10.1f\n', N_k, ops_dft, ops_fft, ops_dft/ops_fft);
end
fprintf('\n');

%% ===================================================
%  第二部分：一维信号的 FFT 频域压缩
%  ===================================================
fprintf('===== 第二部分：一维信号 FFT 频域压缩 =====\n\n');

%% ===== ② 修改这里：定义你的信号 =====
fs = 1000;                    % 采样频率（Hz）
t_end = 1;                    % 信号时长（秒）
t = 0 : 1/fs : t_end - 1/fs; % 时间轴（避免末端重复，共 fs*t_end 个点）
N = length(t);                % 信号长度（= fs * t_end）

% 合成测试信号：三个频率分量叠加（可改成任意信号）
f1 = 50;  A1 = 1.0;    % 50 Hz，幅值1.0
f2 = 120; A2 = 0.6;    % 120 Hz，幅值0.6
f3 = 200; A3 = 0.3;    % 200 Hz，幅值0.3
noise_level = 0.1;     % 加性白噪声的标准差

% 原始信号 = 三个正弦 + 噪声
x_signal = A1*sin(2*pi*f1*t) + A2*sin(2*pi*f2*t) + A3*sin(2*pi*f3*t);
x_noisy  = x_signal + noise_level * randn(size(t));
% randn(size(t))：与 t 同大小的标准正态随机数（均值0，标准差1）
%% =====================================================

% 补零到最近的 2 的幂（使 FFT 最高效）
N_fft = 2^nextpow2(N);   % nextpow2(N)：满足 2^k >= N 的最小 k
% 例：N=1000 → N_fft=1024（2^10）
x_padded = [x_noisy, zeros(1, N_fft - N)];   % 在末尾补零
fprintf('信号长度 N=%d → 补零到 N_fft=%d（最近的2的幂）\n\n', N, N_fft);

% 对补零后的信号做 FFT
X = fft(x_padded);   % 复数频谱，长度为 N_fft
% fft(x)：计算 x 的 DFT，支持任意长度（非2的幂时用混合基 FFT）

% 计算幅值谱（只看正频率半边，因为实信号的频谱是共轭对称的）
freq_axis = (0 : N_fft/2) * fs / N_fft;   % 频率轴（Hz）
amplitude = abs(X(1 : N_fft/2 + 1)) / N_fft;   % 幅值（归一化）
% abs(X)：复数的模（幅值）；/ N_fft 归一化为真实幅值

fprintf('频谱峰值位置：\n');
[~, peak_idx] = sort(abs(X(1:N_fft/2)), 'descend');
% sort(v, 'descend')：降序排列，返回排序后的值和原始下标
for i = 1:5
    freq_peak = (peak_idx(i)-1) * fs / N_fft;
    fprintf('  峰值%d：频率 %.1f Hz，幅值 %.4f\n', i, freq_peak, abs(X(peak_idx(i)))/N_fft);
end
fprintf('\n');

% 三种压缩比的对比
compress_ratios = [0.20, 0.05, 0.01];   % 保留20%、5%、1%的频率系数
rms_errors = zeros(size(compress_ratios));
snr_values  = zeros(size(compress_ratios));

fprintf('%-15s  %-20s  %-15s  %-10s\n', '保留系数比例', '恢复信号 RMS 误差', '信噪比 SNR(dB)', '被置零系数数');
fprintf('%s\n', repmat('-', 1, 64));

for k = 1 : length(compress_ratios)
    ratio = compress_ratios(k);
    keep = round(ratio * N_fft);   % 保留的系数个数

    % 压缩：只保留幅值最大的 keep 个频率系数，其余置零
    X_compressed = zeros(size(X));   % 初始化为全零
    [~, sort_idx] = sort(abs(X), 'descend');   % 按幅值降序排列
    X_compressed(sort_idx(1:keep)) = X(sort_idx(1:keep));   % 保留前 keep 个

    % 重建：对压缩后的频谱做 IFFT
    x_recovered = real(ifft(X_compressed));
    % ifft：逆 FFT；real：取实部（消去因舍入误差产生的极小虚部）
    x_recovered = x_recovered(1:N);   % 截回原始长度

    % 计算误差
    rms_err = sqrt(mean((x_noisy - x_recovered).^2));
    snr_db  = 20 * log10(norm(x_noisy) / norm(x_noisy - x_recovered));
    % norm(v)：向量的2-范数 = sqrt(Σv_i²)
    % 20*log10(信号功率/误差功率) = SNR（分贝）

    rms_errors(k) = rms_err;
    snr_values(k) = snr_db;
    fprintf('%-15.0f%%  %-20.6f  %-15.2f  %-10d\n', ...
        ratio*100, rms_err, snr_db, N_fft - keep);
end

%% ===================================================
%  第三部分：二维 FFT 图像压缩
%  ===================================================
fprintf('\n===== 第三部分：二维 FFT 图像压缩 =====\n\n');

%% ===== ③ 修改这里：可改用实际图像 =====
img_size = 64;   % 图像尺寸（修改这里；用实际图片时替换下面的 img 生成）

% 生成合成图像（实际考题可能用 imread 读入真实图片）
[X_grid, Y_grid] = meshgrid(1:img_size, 1:img_size);
% meshgrid(x,y)：生成网格坐标矩阵
img = sin(2*pi*X_grid/8) + cos(2*pi*Y_grid/16) + ...
      0.5*sin(2*pi*(X_grid+Y_grid)/12);
% 实际用真实图像时用：img = double(imread('filename.png')); 或类似

% 读入真实图像的替代方法（取消注释即可）:
% img = double(imread('lena.png'));   % 读入 PNG 图像（会是彩色，可能需要转灰度）
% img = rgb2gray(img);               % 彩色转灰度
% img = im2double(img);              % 转为 [0,1] 范围的 double
%% =====================================================

% 二维 FFT
F = fft2(img);   % 二维 DFT，等价于对每行再每列各做一次一维 FFT
F_shift = fftshift(F);
% fftshift：将零频分量移到中心（默认 fft2 把零频放在左上角）

% 压缩：保留二维频谱中幅值最大的一部分系数
keep_ratio = 0.05;   % 保留5%的系数
n_total = img_size^2;
n_keep  = round(keep_ratio * n_total);

F_vec = F(:);   % 将二维频谱拉成列向量
[~, sort_idx_2d] = sort(abs(F_vec), 'descend');
F_compressed = zeros(size(F_vec));
F_compressed(sort_idx_2d(1:n_keep)) = F_vec(sort_idx_2d(1:n_keep));
F_compressed = reshape(F_compressed, size(F));   % 还原为二维矩阵
% reshape(v, m, n)：将向量 v 重新排列为 m×n 矩阵

% 重建图像
img_recovered = real(ifft2(F_compressed));   % 二维 IFFT
% ifft2：二维逆 DFT

% 计算误差
rms_2d = sqrt(mean(mean((img - img_recovered).^2)));
% mean(mean(...))：先对各列求均值再对均值向量求均值（= 所有元素的平均）
snr_2d = 20*log10(norm(img(:)) / norm(img(:) - img_recovered(:)));
fprintf('图像压缩（保留 %.0f%% 系数）：RMS=%.4f，SNR=%.2f dB\n\n', ...
    keep_ratio*100, rms_2d, snr_2d);

%% ---- 绘图 ----
figure('Name', 'FFT 信号与图像压缩', 'NumberTitle', 'off', 'Position', [50,50,1400,800]);

% 子图1：原始信号及其频谱
subplot(2, 4, 1);
plot(t, x_noisy, 'b-', 'LineWidth', 0.8);
xlabel('时间 (s)'); ylabel('幅值');
title('原始含噪信号');
grid on;

subplot(2, 4, 2);
plot(freq_axis, amplitude, 'r-', 'LineWidth', 1.2);
xlabel('频率 (Hz)'); ylabel('幅值');
title('幅值频谱（单边）');
xlim([0, fs/2]);
grid on;

% 子图2-4：三种压缩比的重建结果
for k = 1:length(compress_ratios)
    ratio = compress_ratios(k);
    keep = round(ratio * N_fft);
    X_c = zeros(size(X));
    [~, si] = sort(abs(X), 'descend');
    X_c(si(1:keep)) = X(si(1:keep));
    x_rec = real(ifft(X_c));
    x_rec = x_rec(1:N);

    subplot(2, 4, k+2);
    plot(t, x_noisy, 'b-', 'LineWidth', 0.5, 'DisplayName', '原始');
    hold on;
    plot(t, x_rec, 'r-', 'LineWidth', 1.2, 'DisplayName', '重建');
    hold off;
    xlabel('时间 (s)');
    title(sprintf('保留%.0f%%系数\nSNR=%.1f dB', ratio*100, snr_values(k)));
    legend('Location', 'northeast', 'FontSize', 7);
    grid on;
end

% 子图5-7：图像压缩
subplot(2, 4, 5);
imagesc(img); colormap gray; colorbar;
% imagesc：图像显示（自动缩放颜色范围）；colormap gray：灰度颜色映射
title('原始图像');
axis equal tight;

subplot(2, 4, 6);
F_amp = log(1 + abs(F_shift));   % 对幅值取对数（便于可视化，压缩动态范围）
imagesc(F_amp); colormap jet; colorbar;
% colormap jet：彩色热图（红=大，蓝=小）
title('频谱幅值（对数尺度，中心=零频）');
axis equal tight;

subplot(2, 4, 7);
imagesc(img_recovered); colormap gray; colorbar;
title(sprintf('重建图像（保留%.0f%%）\nSNR=%.1f dB', keep_ratio*100, snr_2d));
axis equal tight;

subplot(2, 4, 8);
imagesc(abs(img - img_recovered)); colormap hot; colorbar;
title(sprintf('重建误差图\nRMS=%.4f', rms_2d));
axis equal tight;

sgtitle('FFT 频域压缩：一维信号 + 二维图像');

%% =========================================================
%  局部函数（必须放在脚本末尾）
%  =========================================================

function X = my_fft(x)
    % Cooley-Tukey 基-2 DIT（时间抽取）FFT 递归实现
    %
    % 输入：x - 长度为 N 的向量（N 必须是2的幂次）
    % 输出：X - x 的 DFT（与 MATLAB 内置 fft 结果相同）
    %
    % 算法思路（蝴蝶分解）：
    % 将 N 点 DFT 分解为两个 N/2 点 DFT：
    %   X_k = DFT{x_偶}[k] + W_N^k · DFT{x_奇}[k]
    %   X_{k+N/2} = DFT{x_偶}[k] - W_N^k · DFT{x_奇}[k]
    % 其中旋转因子 W_N^k = e^{-i2πk/N}
    % 递归到 N=1 时 DFT = 本身，直接返回

    N = length(x);
    if N == 1
        X = x;       % 递归基：1点DFT就是本身
        return;
    end
    if mod(N, 2) ~= 0
        error('my_fft: 输入长度必须是2的幂次！当前 N=%d', N);
    end

    % 分解为偶数下标和奇数下标（0-indexed 的 x₀,x₂,... 和 x₁,x₃,...）
    % MATLAB 下标从1开始，所以：
    %   偶数下标（0,2,4,...）→ x(1), x(3), x(5), ... → x(1:2:end)
    %   奇数下标（1,3,5,...）→ x(2), x(4), x(6), ... → x(2:2:end)
    x_even = x(1:2:end);   % 偶数位元素（DIT分解）
    x_odd  = x(2:2:end);   % 奇数位元素

    % 递归计算两个 N/2 点 FFT
    X_even = my_fft(x_even);   % N/2 点 DFT（递归）
    X_odd  = my_fft(x_odd);

    % 旋转因子：W_N^k = e^{-i2πk/N}，k=0,1,...,N/2-1
    k = 0 : N/2 - 1;   % 行向量 [0,1,...,N/2-1]
    W = exp(-1i * 2*pi*k / N);   % 复指数旋转因子（向量化）
    % -1i 是复数 -i（虚数单位的负值）
    % .* 逐元素乘法

    % 蝴蝶合并（Butterfly Combination）
    X = zeros(1, N);
    X(1 : N/2)   = X_even + W .* X_odd;   % 前半段
    X(N/2+1 : N) = X_even - W .* X_odd;   % 后半段（利用周期性省去重复计算）
end
