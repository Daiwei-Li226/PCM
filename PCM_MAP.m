clear all;close all;clc
% 2025/08/26
% use for plot PCM MAP (when only max 3D recon data was saved among sliding-window)

%% Plot PCM sensitivity figure
% load 3D-PCM MAP data
filePath ='C:\Users\daiwe\OneDrive - Duke University\02_PI LAB\02_3D Cavitation_2D array\0614-2024_Bego_Recon_figures\3_Recon\s5_400delay_1.0\2nd\';
name = 's5_400delay_1.0_2nd';
dat_file = dir([filePath,'*.dat']);  % PCM recon data
mat_name = dir([filePath,'*.mat']);  % Recon para

%% Show MAP 
load([filePath,mat_name.name]);
ALL_XYMAP = ones(length(imageArea.x_arr),length(imageArea.y_arr),length(dat_file));

for i = 1:length(dat_file)
    fileName = dat_file(i).name;
    parafileName = mat_name.name;
    load([filePath,parafileName]);
    f = fopen([filePath,'\',fileName],'r');
    temp_img = fread(f,'single');
    temp_img = reshape(temp_img,length(imageArea.x_arr),...
            length(imageArea.y_arr),length(imageArea.z_arr),201);
    fclose(f);

    for a = 1:size(reconTime,2)
    temp_amp(a) = max(temp_img(:,:,:,a),[],'all');
    end
    [~,max_frame] = max(temp_amp);

    temp_img = temp_img(:,:,:,max_frame);

    figure(1);
    subplot(1,3,1);     
    XYMAP = squeeze(max(temp_img,[],1)); 
    imagesc(imageArea.y_arr(:),imageArea.x_arr(:),XYMAP); xlabel('y'); ylabel('x'); colormap("hot"); pbaspect([1 1 1]);
    subplot(1,3,2);     
    XZMAP = squeeze(max(temp_img,[],2));
    imagesc(imageArea.y_arr(:),imageArea.z_arr(:),XZMAP); xlabel('y'); ylabel('z');   colormap("hot"); pbaspect([1 1 1]);
    subplot(1,3,3);     
    YZMAP = squeeze(max(temp_img,[],3)); 
    imagesc(imageArea.x_arr(:),imageArea.z_arr(:),YZMAP);xlabel('x'); ylabel('z'); colormap("hot");  pbaspect([1 1 1]);    

    % max_temp 
    figure(2);plot(temp_amp);

    ALL_XYMAP(:,:,i) = XYMAP;
    
    % whole3dvolume
    ALL_3d_volume(:,:,:,i) = temp_img;
    
    % max amplitude at every pN
    maxamp(i) = max(max(XYMAP));
%     pause(0.5)
end

%% Show MAP for model-based data
% load 3D-PCM MAP data
filePath ='C:\Users\daiwe\Dropbox\PCM_model_based_recon\ModelBased\Davia_modelbased\3_ReconResults\t1\TV\MB_TV\';
% load image area
load('C:\Users\daiwe\Dropbox\PCM_model_based_recon\ModelBased\Davia_modelbased\3_ReconResults\t1\imageRegion.mat');
filePattern = fullfile(filePath, 'MB_Shift_*.mat');
files = dir(filePattern);
filenames = {files.name};
fileNumbers = cellfun(@(x) sscanf(x, 'MB_Shift_%d.mat'), filenames);
[~, sortIdx] = sort(fileNumbers);
sortedFiles = files(sortIdx);

for i = 77 %1:length(sortedFiles)
    baseFileName = sortedFiles(i).name;
    fullFileName = fullfile(filePath, baseFileName);   
    fprintf('Loading: %s\n', baseFileName);
    load(fullFileName);

    temp_img = ModelBaseRecon;

    figure(1);
    subplot(1,3,1);     
    YZMAP = squeeze(max(temp_img,[],3)); 
    imagesc(imageRegion.yArr*1000,imageRegion.xArr*1000,flip(YZMAP));xlabel('y'); ylabel('z'); colormap("hot"); pbaspect([1 1 1]);
    subplot(1,3,2);     
    XZMAP = squeeze(max(temp_img,[],2));
    imagesc(imageRegion.zArr*1000,imageRegion.xArr*1000,flip(XZMAP));xlabel('x'); ylabel('z'); colormap("hot"); pbaspect([1 1 1]);
    subplot(1,3,3);     
    YXMAP = squeeze(max(temp_img,[],1)); 
    imagesc(imageRegion.zArr*1000,imageRegion.yArr*1000,YXMAP);xlabel('x'); ylabel('y'); colormap("hot");  pbaspect([1 1 1]);    
    
    % whole3dvolume
    ALL_3d_volume(:,:,:,i) = temp_img;
    
    % max amplitude at every pN
    maxamp(i) = max(max(XYMAP));
%     pause(0.5)
end
%%
    for a = 1:length(sortedFiles)
    temp_amp(a) = max(ALL_3d_volume(:,:,:,a),[],'all');
    end
    [~,max_frame] = max(temp_amp);
    % max_temp 
    figure(2);plot(temp_amp);


%% Show All pulses XY MAP
accALLXYMAP = squeeze(max(sum(ALL_XYMAP(),3),[],3));
ALLXYMAP = squeeze(max(ALL_XYMAP,[],3));
fig1 = figure;imagesc(imageArea.y_arr(:),imageArea.x_arr(:),accALLXYMAP); xlabel('y'); ylabel('x'); colormap("hot");
title(name,'Accumulated All pulses XY MAP');xlabel('pN');%clim([5e6 17e6])
fig2 = figure;imagesc(imageArea.y_arr(:),imageArea.x_arr(:),ALLXYMAP); xlabel('y'); ylabel('x'); colormap("hot");
title(name,'All pulses XY MAP');xlabel('pN');%clim([1e5 5.5e5])
fig3 = figure;stem(maxamp,'filled');
title(name,'Max amp at individual pN');xlabel('pN');%ylim([0 600000])
% pathfig ='C:\Users\daiwe\OneDrive\桌面\20251219_RT_feedback_RECON\3_recon';
% savefig(fig1,[pathfig,'\',name,'_Accumulated_AllpNXYMAP.fig']);
% savefig(fig2,[pathfig,'\',name,'_AllpNXYMAP.fig']);
% savefig(fig3,[pathfig,'\',name,'_Max_amp.fig']);
% save([pathfig,name,'_results.mat'],'maxamp','ALL_XYMAP');

%% profile
xx = improfile;

figure;
plot(imageArea.x_arr(1:length(xx)),xx);  %imageArea.x_arr(34:length(xx)+33)
hold on;plot(imageArea.y_arr(1:length(yy)),yy)

%% resave the 3Dvolume in case
testName = 'PCM56mhz';
paraName = [testName,'_para.mat'];
% Res = 0.24;
Path = 'C:\Users\daiwe\Dropbox\CavitationMapping_2024_Davia_Chaorui\20250826_PCM_freewater_diffArray';

save([Path,'\',paraName],'imageArea','T','startIdx','reconTime','fs','step','Res','T');

        saveName = [testName];
        f1 = fopen([Path,'\',saveName,'_recon','.dat'],'w');
        fwrite(f1,pa_img_all,'single');
        fclose(f1);




%% Show MAP if ONLY saved all 4D recon data
for i = 1  %:length(dat_file)
    fileName = dat_file(i).name;
    parafileName = mat_name(i).name;
    load([filePath,parafileName]);
    f = fopen([filePath,'\',fileName],'r');

    pa_img_all = fread(f,'single');
    pa_img_all = reshape(pa_img_all,length(imageArea.x_arr),...
            length(imageArea.y_arr),length(imageArea.z_arr),length(reconTime));
    fclose(f);

    for i = 1:size(reconTime,2)
    temp_amp(i) = max(pa_img_all(:,:,:,i),[],'all');
    end
    [~,max_frame] = max(temp_amp);

    temp_img = pa_img_all(:,:,:,max_frame);

     figure;
    subplot(1,3,1);     XYMAP = squeeze(max(temp_img,[],3));
    imagesc(imageArea.y_arr(:),imageArea.x_arr(:),XYMAP); xlabel('y'); ylabel('x');colormap("hot"); pbaspect([1 1 1]);
    subplot(1,3,2);     XZMAP = squeeze(max(temp_img,[],2));
    imagesc(imageArea.z_arr(:),imageArea.x_arr(:),XZMAP); xlabel('z'); ylabel('x');colormap("hot"); pbaspect([1 1 1]);
    subplot(1,3,3);     YZMAP = squeeze(max(temp_img,[],1)); 
    imagesc(imageArea.z_arr(:),imageArea.y_arr(:),YZMAP);xlabel('z'); ylabel('y');colormap("hot");  pbaspect([1 1 1]);    
    
    % max_temp 
    figure(2);plot(temp_amp);

    ALL_XYMAP(:,:,i) = XYMAP;

    pause(0.5)
end


