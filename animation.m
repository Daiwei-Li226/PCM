%% 
% clear all
% clc

%%
%% 1. 数据准备 & 物理坐标系定义
% ------------------------------------------------------------
vol_size = [63 63 63];
voxel_size = 0.08; % mm

% 模拟数据
PCM_data = flip(bbb,3); 
SuperRes_points = flip(ppp_1,1); 

% 定义空间参考 (这步仅用于计算，不传给 volshow)
xWorldLimits = [-(vol_size(2)-1)/2, (vol_size(2)-1)/2] * voxel_size;
yWorldLimits = [-(vol_size(1)-1)/2, (vol_size(1)-1)/2] * voxel_size;
zWorldLimits = [-(vol_size(3)-1)/2, (vol_size(3)-1)/2] * voxel_size;
R = imref3d(vol_size, xWorldLimits, yWorldLimits, zWorldLimits);

%% 2. 核心配准：将散点转为密度矩阵
% ------------------------------------------------------------
[subI, subJ, subK] = worldToSubscript(R, SuperRes_points(:,1), SuperRes_points(:,2), SuperRes_points(:,3));

% 过滤越界点
valid_idx = subI >= 1 & subI <= vol_size(1) & ...
            subJ >= 1 & subJ <= vol_size(2) & ...
            subK >= 1 & subK <= vol_size(3);
subI = subI(valid_idx); subJ = subJ(valid_idx); subK = subK(valid_idx);

% 生成密度图
SR_Density = accumarray([subI, subJ, subK], 1, vol_size);

%% 3. 亮度与对比度增强 (解决"太暗"问题)
% ------------------------------------------------------------
% % % % % A. 处理 PCM: 归一化
% % % % PCM_Show = (PCM_data - min(PCM_data(:))) / (max(PCM_data(:)) - min(PCM_data(:)));

V = PCM_data;
low_val = min(V(:));
high_val = prctile(V(:), 99.95); % use 99% as upper limit (can be 99.5 or 99.9)
%  Clipping
V(V > high_val) = high_val; 
% normalize to [0,1]
V = (V - low_val) / (high_val - low_val);
% Gamma calibration (optional，<1 make the darker region brighter)
PCM_Show = V .^ 1.1;
PCM_Show = PCM_Show * 0.95;



% B. 处理 Density: Log 压缩 + Gamma 校正 (关键！)
% 先做 Log 防止高亮太亮、低亮看不见
SR_Show = log(SR_Density + 1); 

% % % % 归一化到 0-1
% % % SR_Show = (SR_Show - min(SR_Show(:))) / (max(SR_Show(:)) - min(SR_Show(:)));
% % % % Gamma 校正 (0.5次方会让中间灰度变亮)
% % % SR_Show = SR_Show .^ 0.5; 


P = SR_Show;
low_val = min(P(:));
high_val = prctile(P(:), 99.9); % use 99% as upper limit (can be 99.5 or 99.9)
%  Clipping
P(P > high_val) = high_val; 
% normalize to [0,1]
P = (P - low_val) / (high_val - low_val);
% Gamma calibration (optional，<1 make the darker region brighter)
SR_Show = P .^ 1.0;
SR_Show = SR_Show .^ 0.5 * 0.95;

% % % % % % %% 4. 添加"绝对白色"边框
% % % % % % % ------------------------------------------------------------
% % % % % % box_mask = false(vol_size);
% % % % % % box_mask([1, end], [1, end], :) = true; 
% % % % % % box_mask([1, end], :, [1, end]) = true;
% % % % % % box_mask(:, [1, end], [1, end]) = true;
% % % % % % box_mask = imdilate(box_mask, strel('cube', 2)); % 边框粗细
% % % % % % 
% % % % % % % 将边框位置设为 1.0 (由于数据最高只有 0.95，1.0 是独一无二的最大值)
% % % % % % PCM_Show(box_mask) = 1.0;
% % % % % % SR_Show(box_mask)  = 1.0;

% 5. 定义自定义色图 (让最大值 = 白色)
% ------------------------------------------------------------
% 为 PCM 创建修改版的 parula
cmap_pcm = hot(256);
cmap_pcm(end, :) = [1 1 1]; % <--- 强制把色图最后一个颜色设为纯白 [R G B]

% 为 Density 创建修改版的 hot (其实 hot 本身最后就是白，但为了保险也加上)
cmap_sr = parula(256);
cmap_sr(end, :) = [1 1 1];

%% 5. 分开展示 (R2022a 修正版)
% ------------------------------------------------------------
% 这里的关键是：
% 1. 去掉了 R 参数
% 2. 加入了 'ScaleFactors'，告诉 MATLAB 像素比例是 0.2mm



% --- 窗口 1: 3D-PCM ---
h1 = volshow(PCM_Show, ...
    'ScaleFactors', [voxel_size voxel_size voxel_size], ... % <--- R2022a 适配
    'Renderer', 'MaximumIntensityProjection', ...
    'Colormap', cmap_pcm, ...
    'BackgroundColor', [0 0 0]);
% 
% % 调整 Alpha 映射提升内部可见度 (防止只有边框亮)
% h1.Alphamap = linspace(0, 1, 256) .^ 0.5; 

% 调整透明度：注意最后一项(边框)要设为不透明
alpha_map = linspace(0, 0.8, 256) .^ 0.5; 
alpha_map(end) = 1.0; % 让边框(值=256的位置)完全不透明
h1.Alphamap = alpha_map';

% --- 窗口 2: Super-Resolution Density ---
h2 = volshow(SR_Show, ...
    'ScaleFactors', [voxel_size voxel_size voxel_size], ... % <--- R2022a 适配
    'Renderer', 'MaximumIntensityProjection', ...
    'Colormap', parula, ...
    'BackgroundColor', [0 0 0]);

% % 设定一个特定的视角
% targetCameraPosition = [30, -10, 30]; % 物理坐标位置
% targetCameraViewAngle = 60;
% 
% % 强行同步两个图的视角
% h1.CameraPosition = targetCameraPosition;
% h1.CameraViewAngle = targetCameraViewAngle;
% 
% h2.CameraPosition = targetCameraPosition;
% h2.CameraViewAngle = targetCameraViewAngle;