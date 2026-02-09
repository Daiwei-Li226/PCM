%%
clear all; close all; clc; 
g = gpuDevice(1); reset(g);
%% reading RF data
filePath = 'C:\Users\daiwe\Dropbox\CavitationMapping_2024_Davia_Chaorui\20260204_Random\2_PCM\';
testName = 't1';
reconPath = 'C:\Users\daiwe\Dropbox\CavitationMapping_2024_Davia_Chaorui\20260204_Random\PCM_recon';
mkdir(reconPath);

dat_file = dir([filePath,testName,'\','*.dat']);
mat_name = dir([filePath,testName,'\','*.mat']);
load([filePath,testName,'\',mat_name(1).name]);
% load('C:\Users\daiwe\Dropbox\CavitationMapping_2024_Davia_Chaorui\20250807_5MHzV2_PCMsensitivtity\1_PCM\t6_z100_x10\t6_z100_x10_setupParams.mat')
% Manually revise some para (if they were not saved correctly in setupscript)
Trans.frequency = 4.5;  %5.6; 
Trans.elementWidth = 0.5;  %1;

%% Recon 
IFstack = 0;  % means fine recon if 1

for pulse_idx = 20 %length(dat_file)
if IFstack
   load ('C:\Users\daiwe\OneDrive\桌面\20251219_RT_feedback_RECON\2_seg\T10\T10_seg.mat')
end
% idatfilesname = dat_file(size(dat_file,1)).name;
idatfilesname = dat_file(pulse_idx).name;
f = fopen([filePath,testName,'\',idatfilesname],'r');
iRFData = fread(f,'int16');
iRFData = reshape(iRFData,[RFDataSize(1) RFDataSize(2)]);
fclose(f);
RFData = iRFData(Receive(1).startSample:Receive(1).endSample,:);

%   RFData =RFData(6317:7216,:); % minus traveling time: 38.4 mm (496.1779 pts)
% % % % % % % % % % RFData =RFData(1:13000,:);
RFData =RFData(:,:);

% checking RF data
fs = Receive(1).decimSampleRate;
T = 22.0;
sos = (1402.4+5.01*T-0.055*T^2+0.00022*T^3)*1e-3;
timeFs = [1:size(RFData,1)]*(1/fs)*1.5;

if IFstack
else
figure; plot(RFData(:,120)); 
end
% % % sonogram = abs(RFData(:,1:256));
% % % figure
% % % imagesc(1:128,timeFs,sonogram)
% % % colormap("gray")
% % % xlabel('Channel');
% % % ylabel('Traveling dis(mm)')
% % % ylim([0 100])
% % % colorbar

%% init
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

% % 5.6Mhz regular rigid array
% Trans2561.name = 'qcr2561';
% Trans2561.units = 'wavelengths';
% Trans2561 =  computeTrans_256_2D_Martix(Trans2561,sos);    % L7-4 transducer is 'known' transducer so we can use computeTrans.
% Trans2562.name = 'qcr2562';
% Trans2562.units = 'wavelengths'; 
% Trans2562 = computeTrans_256_2D_Martix(Trans2562,sos);
% 
% waveLength = (sos)/Trans.frequency;
% Trans.ElementPosMm = [Trans2561.ElementPosMm; Trans2562.ElementPosMm];
% Trans.ElementPosWL = Trans.ElementPosMm /waveLength;
% elemWidthwl = Trans.elementWidth/(sos/Trans.frequency);  %Used for transSenAll = calc_angular_sens(..) function
% 
% x_trans = Trans.ElementPosMm(:,1);
% y_trans = Trans.ElementPosMm(:,2);
% z_trans = Trans.ElementPosMm(:,3);

figure, hold on
for ielem = 1:256  %size(x_trans,2)
    scatter(x_trans(ielem),y_trans(ielem));
    text(x_trans(ielem),y_trans(ielem),num2str(ielem));
end

%%
sensArr = SensorArray2D(x_trans',y_trans',z_trans');
if IFstack
    factor = 4;
else
    factor = 2;
end
Res = 0.06*factor;    % 60 um- out of memory 
rough_z = 40;
x_range = [-5,5]; 
y_range = [-5,5]; 
z_range = [34,49];   % 4 mm for z axis for fine recon
imageArea = ImageArea2D_matrix(max(x_range), min(x_range), Res, ...
    max(y_range), min(y_range), Res,...
    max(z_range), min(z_range), Res);

% recon
travelingTime = rough_z/sos*fs;
reconStep = 400;  %500

if IFstack
    aastartIdx = startIdx(pulse_idx) - 250-round(travelingTime) ;
else
    aastartIdx = 7037 -250 -round(travelingTime);
end
endIdx = aastartIdx + reconStep;

if IFstack
   step = 3; 
else
   step = 3;  % 4/6 for quick recon
end

reconTime = 0:step:reconStep;

% create sparse-matrix first
rf_temp = RFData(:,:,1);
rf_class = RFclass2D(rf_temp,0,fs,sos);
% mult_mat = IndexMatrix2D_OneWay_AS(transSenAll,sensArr,imageArea,rf_class,"matrix");
if(~exist('mult_mat','var'))
    mult_mat = IndexMatrix2D(sensArr,imageArea,rf_class,"matrix");
    gmul_mat = mult_mat.Send2GPU;
end

%%
RFData = hilbert(RFData);
%%
pa_img = zeros(length(imageArea.x_arr),length(imageArea.y_arr),...
    length(imageArea.z_arr),length(reconTime));
% for iTX = 1:size(RFData,4)
% for iTX = 5:5
%     for iframe = 1:size(RFData,3)
    for iframe = 1:1
        for t = 1:length(reconTime)
            % t
            ishift = step*(t-1)+aastartIdx;
%             rfdata_recon = squeeze(RFData(:,129:256,iframe,iTX));
            rfdata_recon = squeeze(RFData(:,:,iframe));
            rfdata_recon(:,1:128) = 0;                          % set channel 1-128 to zero
            rfdata_recon = circshift(rfdata_recon,[-ishift,0]);
            rf_class = RFclass2D(rfdata_recon,0,fs,sos);

            temp_rec = DAS_mmult(rf_class,gmul_mat);
            temp_rec = reshape(temp_rec,length(imageArea.x_arr),...
                length(imageArea.y_arr),length(imageArea.z_arr));
            pa_img(:,:,:,t) = temp_rec;
            t
        end
    end
% endz

pa_img = abs(pa_img);

pa_img_all = gather(single(pa_img));

for a = 1:size(reconTime,2)
    temp_amp(a) = max(pa_img_all(:,:,:,a),[],'all');
end

[~,max_frame] = max(temp_amp);

    temp_img = abs(pa_img_all(:,:,:,max_frame));


 % Store recon data
    paraName = [testName,'_',num2str(pulse_idx),'_2ndBB','_para.mat'];
    save([reconPath,'\',paraName],'temp_amp','imageArea','T','aastartIdx','reconTime','fs','step','Res','T','sos','fxngen_delay');

        saveName = [testName,'_',num2str(pulse_idx),'_1stBB'];
        f1 = fopen([reconPath,'\',saveName,'.dat'],'w');
%         pa_img_all = gather(single(pa_img));
        fwrite(f1,temp_img,'single');
        fclose(f1);

%%
if IFstack
   video = 0;
else
   video = 1;
end
if video
v = VideoWriter('test','MPEG-4');
v.FrameRate = 150;
open(v);
for t = 1:size(reconTime,2)
    fig = figure(10);
    temp_img = pa_img(:,:,:,t)./max(pa_img(:));  %220000; %
%     temp_img = 10*log10(abs(pa_img(:,:,:,t))./max(pa_img(:)));
    subplot(1,3,1);
    XYMAP = squeeze(max(temp_img,[],3));
    imagesc(imageArea.y_arr,imageArea.x_arr,XYMAP);
    colormap("hot"); pbaspect([1 1 1]); 
    subplot(1,3,2);
    XZMAP = squeeze(max(temp_img,[],2));
    imagesc(imageArea.z_arr,imageArea.x_arr,XZMAP); colormap("hot"); pbaspect([1 1 1]);
    title(['St Idx = ',num2str(reconTime(t)+aastartIdx)]);
%     caxis([0 1]);
    subplot(1,3,3);
    YZMAP = squeeze(max(temp_img,[],1));
    imagesc(imageArea.z_arr,imageArea.y_arr,YZMAP);
    colormap("hot"); pbaspect([1 1 1]); 
    pause(0.05);
    frame = getframe(fig);
    writeVideo(v,frame);
end
close(v);
end
%% MAP
    f1 = figure(5); plot(temp_amp)
    f2=figure(6);
%     temp_img = abs(pa_img(:,:,:,max_frame))./max(pa_img(:));
%   temp_img = 10*log10(abs(pa_img(:,:,:,t))./max(pa_img(:)));   
    temp_img = abs(pa_img(:,:,:,max_frame));

    subplot(1,3,1); 
    XYMAP = squeeze(max(temp_img,[],3));
    imagesc(imageArea.y_arr,imageArea.x_arr,XYMAP);
    xlabel('y'); ylabel('x');
    colormap("hot");
    pbaspect([1 1 1]);
%   caxis([0 1]);
    subplot(1,3,2);
    XZMAP = squeeze(max(temp_img,[],2));
    imagesc(imageArea.z_arr,imageArea.x_arr,XZMAP);
    xlabel('z'); ylabel('x');    
    colormap("hot");
    pbaspect([1 1 1]);
%   title(['St Idx = ',num2str(reconTime(t)+startIdx)]);
%   caxis([0 1]);
    subplot(1,3,3);
    YZMAP = squeeze(max(temp_img,[],1));
    imagesc(imageArea.z_arr,imageArea.y_arr,YZMAP);
    xlabel('z'); ylabel('y');
    colormap("hot");
    pbaspect([1 1 1]);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % length = Res:Res:Res*size(pa_img_all,1);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % x_profile = XYMAP(:,round(end/2));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % y_profile = XYMAP(round(end/2),:);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % figure, hold on
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % plot(length,x_profile); plot(length,y_profile);

% ALL_XYMAP(pulse_idx,:,:) =  XYMAP;
pulse_idx

savefig(f1, fullfile(reconPath, [testName,'_1stBB','_',num2str(pulse_idx), '_window.fig']));
savefig(f2, fullfile(reconPath, [testName,'_1stBB','_',num2str(pulse_idx), '_MAP.fig']));
end



function transSenAll = calc_angular_sens(x_trans,y_trans,z_trans,x_range,y_range,z_range,Res,elemWidthwl,Trans)
x_img = x_range(1):Res:x_range(2);
y_img = y_range(1):Res:y_range(2);
z_img = z_range(1):Res:z_range(2);
x_img0 = repmat(x_img',[1 length(y_img) length(z_img)]);
y_img0 = repmat(y_img,[length(x_img) 1 length(z_img)]);
z_img0 = repmat(z_img',[1 length(x_img) length(y_img)]);
z_img0 = permute(z_img0,[2 3 1]);

% Trans.numelements = 128;    % only half enable in this case

% if ~exist('transSenAll','var')
    angleAll = zeros(length(x_img),length(y_img),length(z_img),Trans.numelements);
    transSenAll = angleAll;
    [xsize,ysize,zsize] = size(x_img0);
    x_img1 = reshape(x_img0,[xsize*ysize*zsize,1]);
    y_img1 = reshape(y_img0,[xsize*ysize*zsize,1]);
    z_img1 = reshape(z_img0,[xsize*ysize*zsize,1]);
    
    taend = 0;
    wb = waitbar(0,'Calculating angular sensitivity map ...');
    for ielem = 1:Trans.numelements
        tastart = tic;
        dotUV = x_trans(ielem)*(x_img1 - x_trans(ielem)) + ...
            y_trans(ielem)*(y_img1 - y_trans(ielem)) + ...
            z_trans(ielem)*(z_img1 - z_trans(ielem));
        normU = norm([x_trans(ielem) y_trans(ielem) z_trans(ielem)]);
        normV = sqrt((x_img1 - x_trans(ielem)).^2 + ...
            (y_img1 - y_trans(ielem)).^2 + ...
            (z_img1 - z_trans(ielem)).^2);
        temp = abs(dotUV./(normU*normV));
        temp(temp>1) = 1;
        %     temp(temp<-1) = -1;
        theta0 = acos(temp);
        angleAll(:,:,:,ielem) = reshape(theta0,[xsize,ysize,zsize]);
        transSenAll(:,:,:,ielem) = abs(cos(angleAll(:,:,:,ielem)).*...
            (sin(elemWidthwl*pi*sin(angleAll(:,:,:,ielem))))./...
            (elemWidthwl*pi*sin(angleAll(:,:,:,ielem))));
        taend = taend + toc(tastart);
        taavg = taend/ielem;
        tarem = (Trans.numelements - ielem)*taavg;
        waitbar(ielem/Trans.numelements,wb,...
            [sprintf('%12.1f',tarem) ' sec remaining calculating angular sensitivity']);
    end
    close(wb)
%     clear x_img1 y_img1 z_img1 temp
% end
end



