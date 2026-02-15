clear all;close all;clc
% 2025/08/26
% use for plot PCM MAP (when only max 3D recon data was saved among sliding-window)

%% Plot PCM sensitivity figure
% load 3D-PCM MAP data
filePath ='C:\Users\daiwe\Dropbox\PCM_model_based_recon\ModelBased\Davia_modelbased\3_ReconResults\camera0501\t1\TV\Reg200\MB_Shift\';
name = 'TV_Reg200';
filePattern = fullfile(filePath, 'MB_Shift_*.mat');
files = dir(filePattern);

%% Show MAP if ONLY saved max 3D recon data
% ALL_XYMAP = ones(length(imageArea.x_arr),length(imageArea.y_arr),length(dat_file));

filenames = {files.name};
fileNumbers = cellfun(@(x) sscanf(x, 'MB_Shift_%d.mat'), filenames);
[~, sortIdx] = sort(fileNumbers);
sortedFiles = files(sortIdx);

for i = 1:length(sortedFiles)
    baseFileName = sortedFiles(i).name;
    fullFileName = fullfile(filePath, baseFileName);
    
    fprintf('Loading: %s\n', baseFileName);
    % Load
    data = load(fullFileName);
    temp_img = data.ModelBaseRecon;

    fig1 = figure(1);
    subplot(1,3,1); 
    YZMAP = squeeze(max(temp_img,[],1)); 
    imagesc(YZMAP);xlabel('x'); ylabel('y'); colormap("hot");  pbaspect([1 1 1]);   
    subplot(1,3,2);     
    XZMAP = squeeze(max(temp_img,[],2));    
    imagesc(XZMAP); xlabel('x'); ylabel('x'); colormap("hot"); pbaspect([1 1 1]);title('MB Recon - TV(Reg200)',num2str(i))
    subplot(1,3,3);     
    XYMAP = squeeze(max(temp_img,[],3));
    imagesc(XYMAP); xlabel('y'); ylabel('z'); colormap("hot"); pbaspect([1 1 1]);

    ALL_XYMAP(:,:,i) = XYMAP;
    % max amplitude at every ishift
    maxamp(i) = max(max(XYMAP));
     pause(0.5)

savefig(fig1,[filePath,num2str(i),'_XYMAP.fig']);
end

%% 

fig2 = figure; %stem(maxamp,'filled');
plot(maxamp);
title(name,'Max amp at individual ishift');xlabel('Recon Instance #'); %ylim([0 600000])
savefig(fig2,[filePath,'\',name,'_Max_amp.fig']);

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
    fileName = mat_name(i).name;
    load([filePath,fileName]);
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
    imagesc(imageArea.y_arr(:),imageArea.x_arr(:),XYMAP); xlabel('y'); ylabel('x'); colormap("hot"); pbaspect([1 1 1]);
    subplot(1,3,2);     XZMAP = squeeze(max(temp_img,[],2));
    imagesc(imageArea.z_arr(:),imageArea.x_arr(:),XZMAP); xlabel('z'); ylabel('x');   colormap("hot"); pbaspect([1 1 1]);
    subplot(1,3,3);     YZMAP = squeeze(max(temp_img,[],1)); 
    imagesc(imageArea.z_arr(:),imageArea.y_arr(:),YZMAP);xlabel('z'); ylabel('y'); colormap("hot");  pbaspect([1 1 1]);    
    
    % max_temp 
    figure(2);plot(temp_amp);

    ALL_XYMAP(:,:,i) = XYMAP;

    pause(0.5)
end


