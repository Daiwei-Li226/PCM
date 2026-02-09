clear, clc
g = gpuDevice(1);
reset(g);
% 
dirName = 'C:\Users\daiwe\Dropbox\CavitationMapping_2024_Davia_Chaorui\20260204_Random\3_USphantom\';
fileName = 't1__US.mat';
saving_path = [dirName,'DMAS_Results',filesep];
if ~exist(saving_path, 'dir')
    mkdir(saving_path)
end
load([dirName fileName]);

%% Int
% Manually revise some para
Trans.frequency = 4.5; %4.5;   %5.6; 
Trans.elementWidth = 0.5; %0.5;    %1;
Trans.numelements = 256;
Res = 0.12*2;
Res = 0.3;
x_range = [-20,20];
y_range = [-20,20];
z_range = [30,120];
x_img = x_range(1):Res:x_range(2);
y_img = y_range(1):Res:y_range(2);
z_img = z_range(1):Res:z_range(2);

%% reading sensor location; unit: mm
temperature_c= 21;
sos= (1402.4+5.01*temperature_c-0.055*temperature_c^2+0.00022*temperature_c^3);
fs = Receive(1).decimSampleRate*1e6;

% 4.5 MHz random-pattern array
load("C:\Users\daiwe\Dropbox\CavitationMapping_2024_Davia_Chaorui\20260204_Random\run001_eval0000774_centers.mat");
waveLength = (sos/(Trans.frequency*1e6))*1e3; % mm
elemWidthwl = Trans.elementWidth/waveLength;
Trans.ElementPosMm = elem_centers_mm;
Trans.ElementPosMm(:,3)=0;
Trans.ElementPosWL = Trans.ElementPosMm /waveLength;
x_trans = Trans.ElementPosMm(:,1); % mm
y_trans = Trans.ElementPosMm(:,2);
z_trans = Trans.ElementPosMm(:,3);
M = reshape(x_trans, 16, 16);  
N = reshape(y_trans, 16, 16);   
M =rot90(M,-1);  
N =rot90(N,-1);

x_trans = reshape(M, [], 1);
y_trans = reshape(N, [], 1);


% % Re-compute Transducer info when using 5.6Mhz regular rigid array
% Trans2561.name = 'qcr2561';
% Trans2561.units = 'wavelengths';
% Trans2561 =  computeTrans_256_2D_Martix(Trans2561,sos);   
% Trans2562.name = 'qcr2562';
% Trans2562.units = 'wavelengths'; 
% Trans2562 = computeTrans_256_2D_Martix(Trans2562,sos);
% waveLength = (sos/1e3)/Trans.frequency;
% elemWidthwl = Trans.elementWidth/(sos/1e3/Trans.frequency);  
% Trans.ElementPosMm = [Trans2561.ElementPosMm; Trans2562.ElementPosMm];
% Trans.ElementPosWL = Trans.ElementPosMm /waveLength;
% x_trans = Trans.ElementPosMm(:,1);
% y_trans = Trans.ElementPosMm(:,2);
% z_trans = Trans.ElementPosMm(:,3);

figure, scatter3(x_trans,y_trans,z_trans);xlabel('x (mm)');ylabel('y (mm)');zlabel('z (mm)');
figure, hold on
for ielem = 1:256  %size(x_trans,1)
    scatter(x_trans(ielem),y_trans(ielem));
    text(x_trans(ielem),y_trans(ielem),num2str(ielem));
end

%% get RF data
clear RFData
RFData = RcvData{1};
RFData = RFData(1:Receive(end).endSample,:,:);

%% Prepare 3D Recon
savedata = true; 
mat = [x_trans(:), y_trans(:), z_trans(:)];
ptCloud = pointCloud(mat);
normals = pcnormals(ptCloud,3);
nx = normals(:,1);
ny = normals(:,2);
nz = normals(:,3);

%% Reconstruction Loop
clear transSenAll rRCV
% init
if ~exist('transSenAll','var')
    transSenAll = calc_angular_sens(x_trans,y_trans,z_trans,x_img,...
        y_img,z_img,elemWidthwl,Trans,nx,ny,nz);
    transSenAll = gpuArray(single(transSenAll));
end
if ~exist('rRCV','var')
    rRCV = calc_deley_mat(x_trans,y_trans,z_trans,x_img,...
        y_img,z_img,Trans);
end

TXNum = size(TX,2);
RFshift = TW.peak/(Trans.frequency)*Receive(1).decimSampleRate;

% Initialize DAS & DMAS matrix
framenum = size(RFData,3);
us_rec_DAS  = single(zeros(length(x_img),length(y_img),length(z_img),framenum));
us_rec_DMAS = single(zeros(length(x_img),length(y_img),length(z_img),framenum));
tic

% Parameters for sensitivity simulation
num_coefficients = 256; 
mean_val = 1.0;         
std_dev = 0.05;         
seed = 42;              
weight = generateSensitivityCoefficients(num_coefficients, mean_val, std_dev, seed);

for itx = 1:TXNum
    itx
    alineStart = Receive(itx).startSample;
    alineEnd = Receive(itx).endSample;
    
    % Apply weight and shift
    iRF = weight(itx) * single(RFData(alineStart:alineEnd,:,:));
    iRF = imtranslate(iRF,[0 -RFshift],'FillValues',0);
    iRF = gpuArray(single(hilbert(iRF))); % Complex analytical signal
    
    txelem = find(TX(itx).Apod>0);
    
    % Calculate Delay Indices
    idxAll = gpuArray(single(round((repmat(rRCV(:,:,:,txelem),[1 1 1 Trans.numelements])...
        + rRCV - 2*P.startDepth*waveLength)/1e3/sos*fs)));
    
    for iframe = 1:framenum
        % DMAS Initialization  DMAS_term = us_dmas_term = us_sum_gpu.^2 - us_sq_sum_gpu;
        % two accumulaters needed ：one for Sum(us_sq_sum_gpu)，one for Sum(us_sum_gpu^2)      
        us_sum_gpu    = gpuArray(single(zeros(size(us_rec_DAS(:,:,:,1))))); 
        us_sq_sum_gpu = gpuArray(single(zeros(size(us_rec_DAS(:,:,:,1)))));

        for ielem = 1:Trans.numelements
            [us_sum_gpu, us_sq_sum_gpu] = GPUReconLoop_DMAS(iRF, idxAll, iframe, ielem, us_sum_gpu, us_sq_sum_gpu, transSenAll);
        end
        
        % DAS
        us_rec_DAS(:,:,:,iframe) = us_rec_DAS(:,:,:,iframe) + single(gather(us_sum_gpu));
        
        % DMAS 
        % Algebraic Identity: DMAS = (Sum)^2 - Sum(Square)
        us_dmas_term = us_sum_gpu.^2 - us_sq_sum_gpu;
        
%         is_coherent = abs(us_sum_gpu).^2 > abs(us_sq_sum_gpu);
%         us_bmode_dmas_mag = sqrt(max(0, abs(us_sum_gpu).^2 -us_sq_sum_gpu));

        % Dimensionality restoration (Signed Square Root)
        % y = sign(x) * sqrt(abs(x))，sign(x) = x ./ abs(x)
        us_bmode_dmas = (us_dmas_term ./ (abs(us_dmas_term) + 1e-12)) .* sqrt(abs(us_dmas_term));
        
        % Accumulate frames (if simple summing is desired)
        us_rec_DMAS(:,:,:,iframe) = us_rec_DMAS(:,:,:,iframe) + single(gather(us_bmode_dmas));
%         us_rec_DMAS(:,:,:,iframe) = gather(single(us_bmode_dmas_mag));
%         us_rec_DAS(:,:,:,iframe) = gather(single(abs(us_sum_gpu)));

    end
end
toc

%% Display MAP
CR = 15; a = 0;
vol_to_plot = us_rec_DMAS;
xymap = squeeze(max(abs(mean(vol_to_plot(:,:,:,:),4)),[],3))';
xzmap = squeeze(max(abs(mean(vol_to_plot(:,:,:,:),4)),[],1))';
yzmap = squeeze(max(abs(mean(vol_to_plot(:,:,:,:),4)),[],2))';

fig1 = figure();
subplot(1,3,1),imagesc(x_img,y_img,20*log10(xymap/max(xymap(:))))
colormap gray,caxis([-CR a]) 
pbaspect([length(x_img)/length(y_img) 1 1]),xlabel('x [mm]'),ylabel('y [mm]'), title('DMAS XY')
subplot(1,3,2),imagesc(x_img,z_img,20*log10(xzmap/max(xzmap(:))))
colormap gray,caxis([-CR a]) 
pbaspect([length(y_img)/length(z_img) 1 1]),xlabel('y [mm]'),ylabel('z [mm]'), title('DMAS YZ')
subplot(1,3,3),imagesc(y_img,z_img,20*log10(yzmap/max(yzmap(:))))
colormap gray,caxis([-CR a]) 
pbaspect([length(x_img)/length(z_img) 1 1]),xlabel('x [mm]'),ylabel('z [mm]'), title('DMAS XZ')
drawnow

vol_to_plot2 = us_rec_DAS;
xymap2 = squeeze(max(abs(mean(vol_to_plot2(:,:,0:250,:),4)),[],3))';
xzmap2 = squeeze(max(abs(mean(vol_to_plot2(:,:,200:250,:),4)),[],1))';
yzmap2 = squeeze(max(abs(mean(vol_to_plot2(:,:,200:250,:),4)),[],2))';

fig2 = figure();
subplot(1,3,1),imagesc(x_img,y_img,20*log10(xymap2/max(xymap2(:))))
colormap gray,caxis([-CR a]) 
pbaspect([length(x_img)/length(y_img) 1 1]),xlabel('x [mm]'),ylabel('y [mm]'), title('DAS XY')
subplot(1,3,2),imagesc(x_img,z_img,20*log10(xzmap2/max(xzmap2(:))))
colormap gray,caxis([-CR a]) 
pbaspect([length(y_img)/length(z_img) 1 1]),xlabel('y [mm]'),ylabel('z [mm]'), title('DAS YZ')
subplot(1,3,3),imagesc(y_img,z_img,20*log10(yzmap2/max(yzmap2(:))))
colormap gray,caxis([-CR a]) 
pbaspect([length(x_img)/length(z_img) 1 1]),xlabel('x [mm]'),ylabel('z [mm]'), title('DAS XZ')
drawnow

if savedata == true
    fileNameBase = erase(fileName,'.mat');
    exportgraphics(fig1,[saving_path,fileNameBase,'_DMAS.tif']);   
    savefig(fig1, fullfile(saving_path, [fileNameBase,'_DMAS_MAP.fig']));
    exportgraphics(fig2,[saving_path,fileNameBase,'_DAS.tif']);   
    savefig(fig2, fullfile(saving_path, [fileNameBase,'_DAS_MAP.fig']));
    % Save both DAS & DMAS
    save([saving_path,fileNameBase,'_DMAS_recon.mat'],'us_rec_DAS','us_rec_DMAS','x_img','y_img','z_img','Res','fs','temperature_c');
else
    exportgraphics(fig1,[saving_path,fileName,'_DMAS.tif']); 
end

%% Compare plot：DAS vs DMAS (Lateral Profile)

frame_idx = 5; 

vol_DAS_env  = double(abs(us_rec_DAS(:,:,:,frame_idx)));
vol_DMAS_env = double(abs(us_rec_DMAS(:,:,:,frame_idx)));

% Find the maximum point
[~, max_idx] = max(vol_DMAS_env(:));
[idx_x, idx_y, idx_z] = ind2sub(size(vol_DMAS_env), max_idx);

idx_y = 63;
idx_z = 35;

% Get lateral screenshot
profile_DAS  = squeeze(vol_DAS_env(:, idx_y, idx_z));
profile_DMAS = squeeze(vol_DMAS_env(:, idx_y, idx_z));

% Normalize or not ?
profile_DAS  = a;          %profile_DAS;
profile_DMAS = b;            %  profile_DMAS;
% profile_DAS  = 20*log10(profile_DAS / max(profile_DAS(:)));
% profile_DMAS = 20*log10(profile_DMAS / max(profile_DMAS(:)));

% fwhm_das  = calculate_fwhm(x_img, profile_DAS);
% fwhm_dmas = calculate_fwhm(x_img, profile_DMAS);
[fwhm_das, x_L_das, x_R_das, y_half_das]   = calculate_fwhm_stats(x_img, profile_DAS);
[fwhm_dmas, x_L_dmas, x_R_dmas, y_half_dmas] = calculate_fwhm_stats(x_img, profile_DMAS);
fprintf('DAS FWHM: %.3f mm | DMAS FWHM: %.3f mm\n', fwhm_das, fwhm_dmas);

fig2 = figure('Name', 'DAS vs DMAS Profile with FWHM', 'Color', 'w');

yyaxis left
plot(x_img, profile_DAS, 'b--', 'LineWidth', 1.5);
hold on;
plot([x_L_das, x_R_das], [y_half_das, y_half_das], 'b-', 'LineWidth', 2); 
text(x_R_das + 0.5, y_half_das, sprintf('%.2f mm', fwhm_das), 'Color', 'b', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('DAS Amplitude (Linear)', 'Color', 'b');
ylim([0, max(profile_DAS)*1.2]); 
ax = gca; ax.YColor = 'b';

yyaxis right
plot(x_img, profile_DMAS, 'r-', 'LineWidth', 2);
hold on;
plot([x_L_dmas, x_R_dmas], [y_half_dmas, y_half_dmas], 'r-', 'LineWidth', 2);
text(x_R_dmas + 0.5, y_half_dmas, sprintf('%.2f mm', fwhm_dmas), 'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('DMAS Amplitude (Linear)', 'Color', 'r');
ylim([0, max(profile_DAS)*1.2]); 
ax = gca; ax.YColor = 'r';

grid on;
xlabel('Lateral Position x [mm]');
title({['Lateral Profile at z=',num2str(z_img(idx_z)),' mm'], ...
       ['Resolution Improvement: ' sprintf('%.1f%%', (fwhm_das - fwhm_dmas)/fwhm_das * 100)]});

legend({'DAS Signal', 'DAS FWHM', 'DMAS Signal', 'DMAS FWHM'}, 'Location', 'best');

if savedata == true
     exportgraphics(fig2,[saving_path,fileNameBase,'_FWHM_Drawn.tif']);
end

%% Functions

function rRCV = calc_deley_mat(x_trans,y_trans,z_trans,x_img,y_img,z_img,Trans)
    x_img0 = repmat(x_img',[1 length(y_img) length(z_img)]);
    y_img0 = repmat(y_img,[length(x_img) 1 length(z_img)]);
    z_img0 = repmat(z_img',[1 length(x_img) length(y_img)]);
    z_img0 = permute(z_img0,[2 3 1]);

    rRCV = zeros(length(x_img),length(y_img),length(z_img),Trans.numelements);
    wb = waitbar(0,'Calculating RCV dist. map ...');
    for ielem = 1:Trans.numelements
        tdstart = tic;
        rRCV(:,:,:,ielem) = sqrt((x_img0 - x_trans(ielem)).^2 + ...
            (y_img0 - y_trans(ielem)).^2 + ...
            (z_img0 - z_trans(ielem)).^2);
        tdend = toc(tdstart); % Fixed timing logic slightly for simplicity
        if mod(ielem, 10) == 0
             waitbar(ielem/Trans.numelements,wb);
        end
    end
    close(wb)
end

function transSenAll = calc_angular_sens(x_trans,y_trans,z_trans,x_img,y_img,z_img,elemWidthwl,Trans,nx,ny,nz)
    x_img0 = repmat(x_img',[1 length(y_img) length(z_img)]);
    y_img0 = repmat(y_img,[length(x_img) 1 length(z_img)]);
    z_img0 = repmat(z_img',[1 length(x_img) length(y_img)]);
    z_img0 = permute(z_img0,[2 3 1]);

    angleAll = zeros(length(x_img),length(y_img),length(z_img),Trans.numelements);
    transSenAll = angleAll;
    [xsize,ysize,zsize] = size(x_img0);
    x_img1 = reshape(x_img0,[xsize*ysize*zsize,1]);
    y_img1 = reshape(y_img0,[xsize*ysize*zsize,1]);
    z_img1 = reshape(z_img0,[xsize*ysize*zsize,1]);
    
    wb = waitbar(0,'Calculating angular sensitivity map ...');
    for ielem = 1:Trans.numelements
        dotUV = nx(ielem)*(x_img1 - x_trans(ielem)) + ...
            ny(ielem)*(y_img1 - y_trans(ielem)) + ...
            nz(ielem)*(z_img1 - z_trans(ielem));
        normU = norm([nx(ielem) ny(ielem) nz(ielem)]);
        normV = sqrt((x_img1 - x_trans(ielem)).^2 + ...
            (y_img1 - y_trans(ielem)).^2 + ...
            (z_img1 - z_trans(ielem)).^2);
        temp = abs(dotUV./(normU*normV));
        temp(temp>1) = 1;
        theta0 = acos(temp) + 0.0001;
        angleAll(:,:,:,ielem) = reshape(theta0,[xsize,ysize,zsize]);
        transSenAll(:,:,:,ielem) = abs(cos(angleAll(:,:,:,ielem)).*...
            (sin(elemWidthwl*pi*sin(angleAll(:,:,:,ielem))))./...
            (elemWidthwl*pi*sin(angleAll(:,:,:,ielem))));
        if mod(ielem, 10) == 0
            waitbar(ielem/Trans.numelements,wb);
        end
    end
    close(wb)
end

% --- Modified for DMAS ---
function [us_sum, us_sq_sum] = GPUReconLoop_DMAS(rfdataAll, idxAll, iframe, iline, us_sum, us_sq_sum, transSenAll)
    rfdata = rfdataAll(:,iline,iframe); % This is complex (hilbert)
    
    % Get the signal for this element (Delay)
    sig_val = interp1(rfdata, idxAll(:,:,:,iline), 'linear', 0);
    
    % Apply Apodization/Sensitivity
    sig_val = sig_val .* transSenAll(:,:,:,iline);
    
    % Accumulate Sum: Sum(s_i)
    us_sum = us_sum + sig_val;
    
    % Accumulate Square Sum: Sum(s_i^2)
    % Note: For complex numbers z, z.^2 is (x+iy)^2 = x^2-y^2 + i2xy. 
    % This is correct for the algebraic expansion of (Sum(z))^2.
    us_sq_sum = us_sq_sum + sig_val.^2; 
end

function coefficients = generateSensitivityCoefficients(num_coefficients, mean_val, std_dev, seed)
    if nargin < 4
        rng('shuffle'); 
    else
        rng(seed); 
    end
    coefficients = mean_val + std_dev .* randn(1, num_coefficients);
end

function [width, x_left, x_right, half_val] = calculate_fwhm_stats(x_axis, y_signal)

    [max_val, max_idx] = max(y_signal);
    half_val = max_val / 2;
    
    idx_left = find(y_signal(1:max_idx) < half_val, 1, 'last');
    if isempty(idx_left)
        x_left = x_axis(1); 
    else
        % interpolation
        x1 = x_axis(idx_left);      y1 = y_signal(idx_left);
        x2 = x_axis(idx_left+1);    y2 = y_signal(idx_left+1);
        x_left = x1 + (half_val - y1) * (x2 - x1) / (y2 - y1);
    end

    idx_right_temp = find(y_signal(max_idx:end) < half_val, 1, 'first');
    if isempty(idx_right_temp)
        x_right = x_axis(end);
    else
        idx_right = max_idx + idx_right_temp - 1;
        x1 = x_axis(idx_right-1);   y1 = y_signal(idx_right-1);
        x2 = x_axis(idx_right);     y2 = y_signal(idx_right);
        x_right = x1 + (half_val - y1) * (x2 - x1) / (y2 - y1);
    end
    
    width = x_right - x_left;
end