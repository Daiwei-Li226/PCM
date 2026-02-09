%% 
% Can 3D rendering 3D VOLUME (double only)
% 2026

%% try plot with box ourside 
V = double(temp_img); 

low_val = min(V(:));
high_val = prctile(V(:), 99.95); % use 99% as upper limit (can be 99.5 or 99.9)

%  Clipping
V(V > high_val) = high_val; 

% normalize to [0,1]
V = (V - low_val) / (high_val - low_val);

% Gamma calibration (optional，<1 make the darker region brighter)
V = V .^ 2; 
% 1. 准备边框 Mask (和之前一样)
vol_size = [67 67 67];
box_mask = false(vol_size);

% 定义 12 条棱
box_mask([1, end], [1, end], :) = true; % Z 轴棱
box_mask([1, end], :, [1, end]) = true; % Y 轴棱
box_mask(:, [1, end], [1, end]) = true; % X 轴棱

% 稍微加粗一点点，保证渲染时能看见
box_mask = imdilate(box_mask, strel('cube', 1));

% 2. 方案：直接将边框"烙印"在数据中

% --- 处理 PCM 数据 ---
PCM_Show = V; % 复制一份用于展示
% 把边框位置的像素设为 PCM 数据的最大值 (通常是 1 或数据的 max)，让它显示为最亮
PCM_Show(box_mask) = max(V(:)); 

% R2022a 展示 PCM (不使用 OverlayData 参数)
h1 = volshow(PCM_Show, ...
    'Renderer', 'VolumeRendering', ...
    'Colormap', hot, ...
    'BackgroundColor', [0 0 0]);

% savefig(gcf, 'Volshow_Volume.fig');  %how to save



%% Volshow (manually adjust view angle)
% V = double(us_tatol_DMAS); 
V = double(bbb); 

low_val = min(V(:));
high_val = prctile(V(:), 99.95); % use 99% as upper limit (can be 99.5 or 99.9)
%  Clipping
V(V > high_val) = high_val; 
% normalize to [0,1]
V = (V - low_val) / (high_val - low_val);
% Gamma calibration (optional，<1 make the darker region brighter)
V = V .^ 1.2;


% 准备边框 Mask 
vol_size = [63 63 63];
box_mask = false(vol_size);
% 定义 12 条棱
box_mask([1, end], [1, end], :) = true; % Z 轴棱
box_mask([1, end], :, [1, end]) = true; % Y 轴棱
box_mask(:, [1, end], [1, end]) = true; % X 轴棱
% 稍微加粗一点点，保证渲染时能看见
box_mask = imdilate(box_mask, strel('cube', 1));

 

PCM_Show = V; % 复制一份用于展示
% 把边框位置的像素设为 PCM 数据的最大值 (通常是 1 或数据的 max)，让它显示为最亮
PCM_Show(box_mask) = max(V(:)); 

alpha = [0 0 0.7 1.0];
color = hot(256);        %gray(256)


% % %% 4. 添加边框 (R2022a "烙印法")
% % % ------------------------------------------------------------
% % box_mask = false(vol_size);
% % box_mask([1, end], [1, end], :) = true; 
% % box_mask([1, end], :, [1, end]) = true;
% % box_mask(:, [1, end], [1, end]) = true;
% % box_mask = imdilate(box_mask, strel('cube', 2)); % 稍微加粗
% % 
% % % 将边框烙印进数据 (设为 1.0，即最亮)
% % PCM_Show(box_mask) = 1.0;
% % SR_Show(box_mask)  = 1.0;

vol = volshow(PCM_Show, 'Colormap', color, 'Renderer', 'MaximumIntensityProjection'); 


%% Previous rendering3D
% Load data
% Load your volume data
data = load('AccMAP3d_3000_Mip.mat');
vol = data.AccMAP3d_3000_Mip;

%%
vol = abs(mean(us_rec,4));
% Normalize volume to [0, 1]
vol = vol - min(vol(:));
vol = vol / max(vol(:));

% Display using volshow with MIP rendering style
% figure;
hVol = volshow(vol, ...
    'RenderingStyle', 'MaximumIntensityProjection', ... % <-- 这是新版的写法
    'Colormap', gray(256), ... % 'BackgroundColor', [0 0 0], ...
    'Alphamap', linspace(0, 1, 256));

viewer = hVol.Parent;

 viewer.BackgroundColor = [0 0 0];
viewer.GradientColor = [0.2 0.2 0.2];


% Optional: adjust view
view(3);
camorbit(30, 15); % change camera angle
title('3D MIP Rendering of Acc3dMAP (hot colormap)', 'Color', 'w');


%% when rending US data
% Get smaller size
vol_2=vol(:,:,30:96);
