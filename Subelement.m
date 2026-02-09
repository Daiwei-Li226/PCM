%% =========================================================================
% Two-Stage Delay-and-Sum (DAS) Beamforming for 16x16 Array
% 目的: 模拟 3x3 虚拟子阵列聚焦，并测试 (10mm, 0mm, 80mm) 处的点目标
% =========================================================================

clear all; close all;

%% 1. 参数设置 (Parameters Setup)

% 阵列参数 (Array Parameters)
N_rows = 16;
N_cols = 16;
N_elements = N_rows * N_cols; % 256
Pitch = 2.0e-3; % 阵元间距 (m)

% 子阵列参数 (Sub-array Parameters)
Sub_size = 3; % 3x3 子阵列
Sub_half = floor(Sub_size / 2);

% 环境和信号参数 (Environment and Signal Parameters)
c = 1540; % 声速 (m/s)
Center_Freq = 5.6e6; % 中心频率 (Hz)
Fs = 4 * Center_Freq; % 采样频率 (Hz), 至少四倍
dt = 1/Fs; % 采样时间间隔

% 阵元中心坐标 (Element Center Coordinates)
x_coords = ((1:N_cols) - (N_cols+1)/2) * Pitch; % x坐标
y_coords = ((1:N_rows) - (N_rows+1)/2) * Pitch; % y坐标
[X_e, Y_e] = meshgrid(x_coords, y_coords);
E_pos = [X_e(:), Y_e(:), zeros(N_elements, 1)]; % 阵元位置 (x, y, z=0)

% 图像空间 (Image Grid)
Z_max = 100e-3; % 最大深度 (m), 确保覆盖 80mm
Nz = 1024; % 深度像素点数
Nx = 128; % 横向像素像素点数
Ny = 64;  % 纵向像素点数

z_grid = linspace(0.5e-3, Z_max, Nz); 
x_grid = linspace(min(x_coords), max(x_coords), Nx); 
y_grid = linspace(min(y_coords), max(y_coords), Ny); 

% 焦点测试位置 (Focus Test Position)
Focus_Test_Point = [10e-3, 0e-3, 80e-3]; % (10mm, 0mm, 80mm)

% 初始化最终图像 (Initialize final image)
B_mode_image = zeros(Nz, Ny, Nx); 

%% 2. 模拟 RF 数据 (Simulated RF Data Generation)

% 如果您有真实的 RF_data，请注释掉这一部分，并加载您的数据。
% 否则，此函数将生成一个在 Focus_Test_Point 处的点目标信号
load('C:\Users\daiwe\Dropbox\PCM_model_based_recon\20251113_PCM_sensitivity_figure\t7.mat');
RF_data = RF;
Time_samples = size(RF_data, 1);

disp(['RF数据已生成/加载，采样点数: ', num2str(Time_samples)]);

%% 3. 两级延迟求和波束形成主循环 (Two-Stage DAS Beamforming)

disp('开始两级 DAS 波束形成...');

for ix = 1:Nx
    x_f = x_grid(ix);
    
    % 使用 parfor (并行计算) 加速 Y 轴循环
    for iy = 1:Ny % 如果不需要并行，将 parfor 改为 for
        y_f = y_grid(iy);
        
        % -----------------------------------------------------------------
        % A. 计算所有阵元到当前焦点 (xf, yf, z_grid) 的距离和延迟
        % -----------------------------------------------------------------
        
        % 阵元到焦点距离矩阵 (Nz x N_elements)
        Dist_mat = zeros(Nz, N_elements);
        for k = 1:N_elements
            R = sqrt((E_pos(k, 1) - x_f)^2 + (E_pos(k, 2) - y_f)^2 + z_grid.^2)';
            Dist_mat(:, k) = R;
        end
        
        % 找到最大距离 (用于计算延迟基准)
        R_max = max(Dist_mat, [], 2); 
        
        % 总延迟矩阵 Tau_total = (R_max - R) / c (Nz x N_elements)
        Tau_total_mat = (R_max - Dist_mat) / c;
        
        % 采样点索引矩阵 Index = floor(Tau / dt) + 1
        Index_mat = floor(Tau_total_mat / dt) + 1;
        Index_mat(Index_mat < 1) = 1;
        Index_mat(Index_mat > Time_samples) = Time_samples; % 避免越界

        % -----------------------------------------------------------------
        % B. 两级延迟求和 (Two-Stage DAS)
        % -----------------------------------------------------------------

        % 初始化子阵列输出 (Nz x N_rows x N_cols)
        % 注意: 我们只在 3x3 块中心 (r,c) 存储有效输出
        Sub_array_output = zeros(Nz, N_rows, N_cols); 

        % 子阵列权重 (Hann 窗加权用于旁瓣抑制)
        W_sub = hann(Sub_size)' .* hann(Sub_size);
        W_sub_vec = W_sub(:); % 9x1 向量
        
        % 遍历子阵列中心点 (r, c) - (从 2 到 15)
        for r = (1 + Sub_half) : (N_rows - Sub_half)
            for c = (1 + Sub_half) : (N_cols - Sub_half)
                
                % 1. 计算 3x3 阵元在 256 阵元中的一维线性索引 (关键修正部分)
                row_idx = (r - Sub_half) : (r + Sub_half);
                col_idx = (c - Sub_half) : (c + Sub_half);
                
                [Row_sub, Col_sub] = meshgrid(row_idx, col_idx);
                % Linear_idx_sub 是 9 个阵元在 256 元素中的一维索引 (9x1)
                Linear_idx_sub = sub2ind([N_rows, N_cols], Row_sub(:), Col_sub(:)); 
                
                % 2. 提取 9 个阵元的总延迟索引和原始 RF 信号
                Idx_sub_mat = Index_mat(:, Linear_idx_sub); % Nz x 9
                RF_sub_mat = RF_data(:, Linear_idx_sub); % Time_samples x 9

                % 3. 子阵列求和 (对每个深度 k_z)
                S_sub = zeros(Nz, 1);
                
                for k_z = 1:Nz
                    % 提取当前深度下的 9 个阵元的信号值
                    Time_indices = Idx_sub_mat(k_z, :); % 1 x 9 (时间采样索引)
                    
                    % 线性索引到 RF_sub_mat: (时间索引, 1:9 列索引)
                    Lin_samples_idx = sub2ind(size(RF_sub_mat), Time_indices', (1:Sub_size^2)');
                    
                    Sample_values_vec = RF_sub_mat(Lin_samples_idx); % 9 x 1 信号值
                    
                    % 应用子阵列加权和求和 (第一级 DAS)
                    S_sub(k_z) = sum(Sample_values_vec .* W_sub_vec);
                end
                
                % 4. 存储子阵列的输出信号
                Sub_array_output(:, r, c) = S_sub;
            end
        end

        % --- 第二级: 全孔径求和 (Full-aperture Summing) ---

        % 在本架构中，由于延迟已在第一级完全施加，第二级仅是简单的求和。
        Beamformed_signal = sum(Sub_array_output, [2, 3]); % 对 R 和 C 维度求和 (Nz x 1)
        
        % -----------------------------------------------------------------
        % C. 结果存储 (Store Results)
        % -----------------------------------------------------------------
        
        B_mode_image(:, iy, ix) = abs(hilbert(Beamformed_signal));
    end
end

disp('波束形成完成。');

%% 4. 图像显示 (Image Display and Verification)

% 找到焦点 X=10mm, Y=0mm 对应的索引
[~, ix_focus] = min(abs(x_grid - Focus_Test_Point(1))); 
[~, iy_focus] = min(abs(y_grid - Focus_Test_Point(2))); 
[~, iz_focus] = min(abs(z_grid - Focus_Test_Point(3))); 

B_mode_db = 20 * log10(abs(B_mode_image) / max(abs(B_mode_image(:))));

figure('Name', '两级DAS波束形成结果');
subplot(2, 2, 1);
% X-Z 截面 (穿过 Y=0mm)
imagesc(x_grid*1000, z_grid*1000, squeeze(B_mode_db(:, iy_focus, :))); 
colormap gray;
colorbar;
axis image;
hold on;
plot(Focus_Test_Point(1)*1000, Focus_Test_Point(3)*1000, 'r+', 'MarkerSize', 15, 'LineWidth', 2); % 标记焦点
hold off;
xlabel('横向位置 X (mm)');
ylabel('深度 Z (mm)');
title('X-Z 截面 (Y=0mm)');
caxis([-60 0]); 

subplot(2, 2, 2);
% Y-Z 截面 (穿过 X=10mm)
imagesc(y_grid*1000, z_grid*1000, squeeze(B_mode_db(:, :, ix_focus))); 
colormap gray;
colorbar;
axis image;
hold on;
plot(Focus_Test_Point(2)*1000, Focus_Test_Point(3)*1000, 'r+', 'MarkerSize', 15, 'LineWidth', 2); % 标记焦点
hold off;
xlabel('纵向位置 Y (mm)');
ylabel('深度 Z (mm)');
title('Y-Z 截面 (X=10mm)');
caxis([-60 0]); 

subplot(2, 2, 3);
% X-Y 截面 (穿过 Z=80mm)
imagesc(x_grid*1000, y_grid*1000, squeeze(B_mode_db(iz_focus, :, :))'); 
colormap gray;
colorbar;
axis image;
hold on;
plot(Focus_Test_Point(1)*1000, Focus_Test_Point(2)*1000, 'r+', 'MarkerSize', 15, 'LineWidth', 2); % 标记焦点
hold off;
xlabel('横向位置 X (mm)');
ylabel('纵向位置 Y (mm)');
title('X-Y 截面 (Z=80mm)');
caxis([-60 0]); 

sgtitle('16x16 阵列两级 DAS 波束形成（点目标在红+处）');