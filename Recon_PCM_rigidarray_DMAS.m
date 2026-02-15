clear, clc
g = gpuDevice(1);
reset(g);

%% load data
dirName = 'C:\Users\daiwe\Dropbox\CavitationMapping_2024_Davia_Chaorui\20250501_5MHz\PCM\t1_glasssurface\';
% saving_path = [dirName,'20260128_DMAS',filesep];
saving_path = 'C:\Users\daiwe\OneDrive\桌面\20260203_PCM_DMAS_test\DMAS\t1_glasssurface';
if ~exist(saving_path, 'dir')
    mkdir(saving_path)
end

dat_file = dir([dirName,'*.dat']);  % PCM raw data
para_file = dir([dirName,'*.mat']);  % PCM raw data
load([dirName,para_file.name]);
 
for i = 11
    fileName = dat_file(i).name;
    f = fopen([dirName,'\',fileName],'r');
    iRFData = fread(f,'int16');
    iRFData = reshape(iRFData,[RFDataSize(1) RFDataSize(2)]);
    fclose(f);
    RFData = iRFData(Receive(1).startSample:Receive(1).endSample,:);
end

video = 0;
figure;plot(RFData(:,120));
%% Initial
Trans.frequency = 5.6; 
Trans.elementWidth = 1;
Trans.numelements = 256;
Res = 0.12*1; 
rough_z = 28;
x_range = [-8,8]; 
y_range = [-8,8]; 
z_range = [20 36];  
x_img = x_range(1):Res:x_range(2);
y_img = y_range(1):Res:y_range(2);
z_img = z_range(1):Res:z_range(2);

%% sensor location & sos
temperature_c= 22;
sos= (1402.4+5.01*temperature_c-0.055*temperature_c^2+0.00022*temperature_c^3);
fs = Receive(1).decimSampleRate*1e6;

Trans2561.name = 'qcr2561'; Trans2561.units = 'wavelengths';
Trans2561 =  computeTrans_256_2D_Martix(Trans2561,sos);   
Trans2562.name = 'qcr2562'; Trans2562.units = 'wavelengths'; 
Trans2562 = computeTrans_256_2D_Martix(Trans2562,sos);
waveLength = (sos/1e3)/Trans.frequency;
Trans.ElementPosMm = [Trans2561.ElementPosMm; Trans2562.ElementPosMm];
x_trans = Trans.ElementPosMm(:,1);
y_trans = Trans.ElementPosMm(:,2);
z_trans = Trans.ElementPosMm(:,3);

%% matrix
fprintf('Calculating One-way Distance Map...\n');
rRCV = calc_one_way_dist(x_trans,y_trans,z_trans,x_img,y_img,z_img,Trans);

%% Sliding Window 
travelingTime = rough_z/sos*fs./1e3;
reconStep = 400;  %500
step_samples = 2; 
start_sample = 10770 -250 -round(travelingTime);
reconTime = 0:step_samples:reconStep;

end_sample   = start_sample + reconStep;
scan_vectors = start_sample : step_samples : end_sample;

fprintf('Sliding Window: %.1f  to %.1f , Total Steps: %d\n', ...
    start_sample, end_sample, length(scan_vectors));

%% 6. PCM Reconstruction Loop (No Angular Sensitivity)

iRF = gpuArray(single(hilbert(RFData))); 

PCM_Map_DAS  = gpuArray(single(zeros(length(x_img),length(y_img),length(z_img))));
PCM_Map_DMAS = gpuArray(single(zeros(length(x_img),length(y_img),length(z_img))));

dist_samples = gpuArray(single(rRCV / 1000 / sos * fs));

tic
wb = waitbar(0, 'PCM Sliding Window Reconstruction (Isotropic)...');

pa_img_DMAS = zeros(length(x_img),length(y_img),...
    length(z_img),length(reconTime));
pa_img_DAS = zeros(length(x_img),length(y_img),...
    length(z_img),length(reconTime));
for i = 1:length(scan_vectors)
    t0_sample = scan_vectors(i); 
    
    img_sum    = gpuArray(double(zeros(size(PCM_Map_DMAS)))); 
    img_sq_sum = gpuArray(double(zeros(size(PCM_Map_DMAS))));
    
    idxAll = t0_sample + dist_samples;
    
    for ielem = 1:Trans.numelements
        [img_sum, img_sq_sum] = GPUReconLoop_PCM_NoAngle(iRF, idxAll, ielem, img_sum, img_sq_sum);
    end
    
    % --- DMAS Calculation ---
    val_dmas_term = img_sum.^2 - img_sq_sum;
    val_dmas = (val_dmas_term ./ (abs(val_dmas_term) + 1e-12)) .* sqrt(abs(val_dmas_term));
    val_das  = img_sum;
    
    pa_img_DMAS(:,:,:,i) = gather(single(abs(val_dmas)));
    pa_img_DAS(:,:,:,i)  = gather(single(abs(img_sum)));
  
    i
    if mod(i, 20) == 0
        waitbar(i/length(scan_vectors), wb, sprintf('PCM Sliding: %d/%d', i, length(scan_vectors)));
    end
   
end
close(wb);
toc

%% Save & video
for a = 1:size(scan_vectors,2)
    temp_amp_DMAS(a) = max(pa_img_DMAS(:,:,:,a),[],'all');
    temp_amp_DAS(a) = max(pa_img_DAS(:,:,:,a),[],'all');    
end
[~,max_frame_DMAS] = max(temp_amp_DMAS);
temp_img_DMAS = abs(pa_img_DMAS(:,:,:,max_frame_DMAS));
[~,max_frame_DAS] = max(temp_amp_DAS);
temp_img_DAS = abs(pa_img_DAS(:,:,:,max_frame_DAS));

 % Store recon data
    paraName = ['T1','_pN11','_parallel','_para_DMAS.mat'];
    save([saving_path,'\',paraName],'temp_amp_DMAS','start_sample','reconTime','fs','step_samples','Res','sos');
    paraName = ['T1','_pN11','_parallel','_para_DAS.mat'];
    save([saving_path,'\',paraName],'temp_amp_DAS','start_sample','reconTime','fs','step_samples','Res','sos');

    saveName = ['T1','_pN11','_parallel_DMAS'];
        f1 = fopen([saving_path,'\',saveName,'.dat'],'w');
%         pa_img_all = gather(single(pa_img));
        fwrite(f1,temp_img_DMAS,'single');
        fclose(f1);
    saveName = ['T1','_pN11','_parallel_DAS'];
        f1 = fopen([saving_path,'\',saveName,'.dat'],'w');
%         pa_img_all = gather(single(pa_img));
        fwrite(f1,temp_img_DAS,'single');
        fclose(f1);

if video
v = VideoWriter('test','MPEG-4');
v.FrameRate = 150;
open(v);
for t = 1:size(scan_vectors,2)
    fig = figure(10);
    temp_img = pa_img_DMAS(:,:,:,t)./550000; %max(pa_img_DMAS(:)); 
%     temp_img = 10*log10(abs(pa_img(:,:,:,t))./max(pa_img(:)));
    subplot(1,3,1);
    XYMAP = squeeze(max(temp_img,[],3));
    imagesc(y_img,x_img,XYMAP);
    colormap("hot"); pbaspect([1 1 1]); 
    subplot(1,3,2);
    XZMAP = squeeze(max(temp_img,[],2));
    imagesc(z_img,x_img,XZMAP); colormap("hot"); pbaspect([1 1 1]);
    title(['DMAS - St Idx = ',num2str(scan_vectors(t)+start_sample)]);
%     caxis([0 1]);
    subplot(1,3,3);
    YZMAP = squeeze(max(temp_img,[],1));
    imagesc(z_img,y_img,YZMAP);
    colormap("hot"); pbaspect([1 1 1]); 
    pause(0.05);
    frame = getframe(fig);
    writeVideo(v,frame);
end
close(v);

k = VideoWriter('test','MPEG-4');
k.FrameRate = 150;
open(k);
for t = 1:size(scan_vectors,2)
    fig = figure(10);
    temp_img = pa_img_DAS(:,:,:,t)./550000; %max(pa_img_DAS(:)); 
%     temp_img = 10*log10(abs(pa_img(:,:,:,t))./max(pa_img(:)));
    subplot(1,3,1);
    XYMAP = squeeze(max(temp_img,[],3));
    imagesc(y_img,x_img,XYMAP);
    colormap("hot"); pbaspect([1 1 1]); 
    subplot(1,3,2);
    XZMAP = squeeze(max(temp_img,[],2));
    imagesc(z_img,x_img,XZMAP); colormap("hot"); pbaspect([1 1 1]);
    title(['DAS - St Idx = ',num2str(scan_vectors(t)+start_sample)]);
%     caxis([0 1]);
    subplot(1,3,3);
    YZMAP = squeeze(max(temp_img,[],1));
    imagesc(z_img,y_img,YZMAP);
    colormap("hot"); pbaspect([1 1 1]); 
    pause(0.05);
    frame = getframe(fig);
    writeVideo(k,frame);
end
close(k);
end

%% MAP
    f1 = figure(5); plot(temp_amp_DMAS)
    f2 = figure(6);  
    temp_img_DMAS = abs(pa_img_DMAS(:,:,:,max_frame_DMAS));

    subplot(1,3,1); 
    XYMAP = squeeze(max(temp_img_DMAS,[],3));
    imagesc(y_img,x_img,XYMAP);
    xlabel('y'); ylabel('x');
    colormap("hot");
    pbaspect([1 1 1]);
%   caxis([0 1]);
    subplot(1,3,2);
    XZMAP = squeeze(max(temp_img_DMAS,[],2));
    imagesc(z_img,x_img,XZMAP);
    xlabel('z'); ylabel('x');    
    colormap("hot");
    pbaspect([1 1 1]);
   title('DMAS');
%   caxis([0 1]);
    subplot(1,3,3);
    YZMAP = squeeze(max(temp_img_DMAS,[],1));
    imagesc(z_img,y_img,YZMAP);
    xlabel('z'); ylabel('y');
    colormap("hot");
    pbaspect([1 1 1]);

savefig(f1, fullfile(saving_path, ['t1','_parallel','_','pN11', '_window' ...
    '_DMAS.fig']));
savefig(f2, fullfile(saving_path, ['t1','_parallel','_','pN11','_MAP_DMAS.fig']));


    f3 = figure(7); plot(temp_amp_DAS)
    f4 = figure(8);  
    temp_img_DAS = abs(pa_img_DAS(:,:,:,max_frame_DAS));

    subplot(1,3,1); 
    XYMAP = squeeze(max(temp_img_DAS,[],3));
    imagesc(y_img,x_img,XYMAP);
    xlabel('y'); ylabel('x');
    colormap("hot");
    pbaspect([1 1 1]);
%   caxis([0 1]);
    subplot(1,3,2);
    XZMAP = squeeze(max(temp_img_DAS,[],2));
    imagesc(z_img,x_img,XZMAP);
    xlabel('z'); ylabel('x');    
    colormap("hot");
    pbaspect([1 1 1]);
   title('DAS');
%   caxis([0 1]);
    subplot(1,3,3);
    YZMAP = squeeze(max(temp_img_DAS,[],1));
    imagesc(z_img,y_img,YZMAP);
    xlabel('z'); ylabel('y');
    colormap("hot");
    pbaspect([1 1 1]);

savefig(f3, fullfile(saving_path, ['t1','_parallel','_','pN11', '_window' ...
    '_DAS.fig']));
savefig(f4, fullfile(saving_path, ['t1','_parallel','_','pN11','_MAP_DAS.fig']));



%% ================= FUNCTIONS =================

function rRCV = calc_one_way_dist(x_trans,y_trans,z_trans,x_img,y_img,z_img,Trans)
    x_img0 = repmat(x_img',[1 length(y_img) length(z_img)]);
    y_img0 = repmat(y_img,[length(x_img) 1 length(z_img)]);
    z_img0 = repmat(z_img',[1 length(x_img) length(y_img)]);
    z_img0 = permute(z_img0,[2 3 1]);

    rRCV = zeros(length(x_img),length(y_img),length(z_img),Trans.numelements);
    for ielem = 1:Trans.numelements
        rRCV(:,:,:,ielem) = sqrt((x_img0 - x_trans(ielem)).^2 + ...
            (y_img0 - y_trans(ielem)).^2 + ...
            (z_img0 - z_trans(ielem)).^2);
    end
end

% [修改点] 移除了 transSenAll 输入，移除了内部的乘法
function [img_sum, img_sq_sum] = GPUReconLoop_PCM_NoAngle(rfdata_gpu, idxAll_gpu, ielem, img_sum, img_sq_sum)
    rf_channel = rfdata_gpu(:, ielem); 
    idx_channel = idxAll_gpu(:,:,:,ielem);
    
    % 线性插值
    sig_val = interp1(rf_channel, idx_channel, 'linear', 0);
    
    
    % DMAS Accumulation
    img_sum    = img_sum + sig_val;
    img_sq_sum = img_sq_sum + sig_val.^2;
end