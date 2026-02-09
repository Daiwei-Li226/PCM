%%
clear all; close all; clc;
%% Check RF Data
fdir = 'C:\Users\daiwe\Dropbox\00_SHIYAN\20251219_RT_feedback\T15\';
idatfilesname = '_0042.dat';
mat_name = dir([fdir,'*.mat']);
load('C:\Users\daiwe\Dropbox\00_SHIYAN\20251219_RT_feedback\T4\T4_setupParams.mat');
% load([fdir,mat_name.name]);
f = fopen([fdir,idatfilesname],'r');
iRFData = fread(f,'int16');
iRFData = reshape(iRFData,[RFDataSize(1) RFDataSize(2)]);  fclose(f);
RFData = iRFData(Receive(1).startSample:Receive(1).endSample,:);
figure,plot(RFData(:,120));xlabel('points')
% figure,plot([1:length(iRFData(:,1))]/15.625,iRFData(:,1));xlabel('time(us)')
%%% save as .mat file if needed
% % savepath = 'E:\0731-2024_test';
% %      save([savepath,'\','2Darray_t6_0010','.mat'],'iRFData')
 %%
global folderSelection
folderSelection = 'T1';
 
 transParaPath = 'C:\Users\daiwe\Dropbox\CavitationMapping_2024_Davia_Chaorui\20250423_5MHz_PCM&US';
 USParaPath = transParaPath;
 rawDataPath = 'C:\Users\daiwe\Dropbox\00_SHIYAN\20251219_RT_feedback';
 segPath = 'C:\Users\daiwe\OneDrive\桌面\20251219_RT_feedback_RECON\2_seg';
 mkdir(segPath);
 reconPath = 'C:\Users\daiwe\OneDrive\桌面\20251219_RT_feedback_RECON\3_recon';
 mkdir(reconPath);
 bubbleInfoPath = 'C:\Users\daiwe\OneDrive\桌面\20251219_RT_feedback_RECON\4_bubble';
 %mkdir(bubbleInfoPath)
 % L1 set up
 L1FlagMAT = false;
 L1SegPic = true;
 L1SD1 = true;
 resScale = 3;
 fiberPos.zdev = 40;   rough_z = fiberPos.zdev;
 fiberPos.xdev = 0;  %
 fiberPos.ydev = 0; % -2-means[-3.5 1.5]
 zmax = 10+fiberPos.zdev; zmin = -10+fiberPos.zdev; zint = 0.08*resScale;
 xmax = 10+fiberPos.xdev; xmin = -10+fiberPos.xdev; xint = zint;
 ymax = 10+fiberPos.ydev; ymin = -10+fiberPos.ydev; yint = zint;
 imageArea = ImageArea2D(zmax, zmin, zint, xmax, xmin, xint, ymax, ymin, yint);  % Z, X, Y
 step = 5;   % recon step  
 dur = 26.4;   % us, the time duration of recon   500 pts = 22us  1067 pts = 47
 pixelSize = zint;
 %% L1 - segmentation
 if(L1FlagMAT)
     PCM_Signal_Segmentation_MATV2(rawDataPath,segPah,L1SegPic,L1SD1);
 else
      PCM_Signal_Segmentation_DATV3_matrix(rawDataPath,segPath,L1SegPic,L1SD1,fxngen_delay,dur);
 %      PCM_Signal_Segmentation_DATV1(rawDataPath,segPath,L1SegPic,L1SD1);
 end
  %%
 L1Int = false;          % 4x interpolation for better image quality
 L1Overlap = false;       % C1+C2 overlap
 L1MAP = false;           % Only store MAP image
 L1Filter = false;       % 10s->100 pulses
 reconT = tic; 

 if(L1FlagMAT)
     PCM_Recon_MAT(rawDataPath,segPath,reconPath,transParaPath,imageArea,step,dur);
 else
%      PCM_Recon_DATV3(rawDataPath,segPath,reconPath,transParaPath,imageArea,step,dur,L1Int,L1Overlap,L1MAP,L1Filter);
     PCM_Recon_DATV3_matrix(rawDataPath,segPath,reconPath,transParaPath,imageArea,step,dur,L1Int,rough_z);
 end
 toc(reconT)
%%
% L2 set up
L2filterThres = 0.8;
L2TrackingFig = true;
L2AniTrack = false;
trackingT = tic;
% PCM_Bubble_InfoV3(bubbleInfoPath,reconPath,rawDataPath,zint,fiberPos,L2filterThres,L2TrackingFig,L2AniTrack,L1Overlap,L1MAP);
PCM_Bubble_Info_V2_matrix(bubbleInfoPath,reconPath,rawDataPath,zint,L2filterThres,L2TrackingFig,L2AniTrack);
toc(trackingT);

%% MAP
run('PALA_SetUpPaths.m')
% loading data 
clc; close all;
profile on

% re-load L1 set up
  reconPath = 'C:\Users\daiwe\OneDrive\桌面\20250511_Recon\3_Recon';
  bubbleInfoPath = 'C:\Users\daiwe\OneDrive\桌面\20250511_Recon\4_Bubble';
 fiberPos.zdev = 0; 
 fiberPos.xdev = 0;  
 fiberPos.ydev = 0; % -2-means[-3.5 1.5]

MAPSavePath= 'C:\Users\daiwe\OneDrive\桌面\20250511_Recon\5_MAP'; 
mkdir(MAPSavePath)
MAPshowflag= 1; 
% Fiber position correction from B-mode
xx = fiberPos.xdev;  %
yy = fiberPos.ydev;  %
zz = fiberPos.zdev;
dataList = dir(reconPath);
dataLen = length(dataList);
dataListBB = dir(bubbleInfoPath);
dataLenBB = length(dataListBB);
  load('C:\Users\daiwe\OneDrive\桌面\20250511_Recon\3_Recon\S8_5\S8_5_para.mat')

 tLen  = length(reconTime);
%
for dataN = 3                                     % Choose correct .dat reconfile dataset folder 
    testName = dataList(dataN).name;
    datList = dir([reconPath,'\',testName,'\*.dat']);
    datListBB = dir([bubbleInfoPath,'\',testName,'\*.Mat']);
b =1;
    for datN = 2 %1:100 %1:100                                    % pN = 33 
         pulseIdx = datList(datN*2-1).name(1:4);
      pN = str2num(pulseIdx);
%         pN = datN;

         fileNameBB = datListBB(datN*2-1).name;
         load([bubbleInfoPath,'\',testName,'\',fileNameBB]);
          track_tot1 = track_tot_1st;
         % load 2nd bubble 
         fileNameBB2 = datListBB(datN*2).name;
         load([bubbleInfoPath,'\',testName,'\',fileNameBB2]);
         track_tot2 = track_tot_2nd;

        % 1st collapse
        fileName = datList(datN*2-1).name;
        f = fopen([reconPath,'\',testName,'\',fileName],'r');
        pa_img_all_type1 = fread(f,'single');
        pa_img_all_type1 = reshape(pa_img_all_type1,length(imageArea.x_arr),...
            length(imageArea.y_arr),length(imageArea.z_arr),tLen);
        fclose(f);

%         % 2nd collapse
%         fileName2 = datList(datN*2).name;
%         f = fopen([reconPath,'\',testName,'\',fileName2],'r');
%         pa_img_all_type2 = fread(f,'single');
%         pa_img_all_type2 = reshape(pa_img_all_type2,length(imageArea.x_arr),...
%             length(imageArea.y_arr),length(imageArea.z_arr),tLen);
%         fclose(f);       

        if MAPshowflag == 1
           % Plot MAP of each collapse point
           % 1st collapse
            maxInt = max(pa_img_all_type1,[],'all');
            mapXY = zeros(length(imageArea.y_arr),length(imageArea.z_arr));
            y= imageArea.y_arr-yy;
            x= -imageArea.z_arr-xx;

            for i_track = 1:size(track_tot1,1)
                temp = track_tot1{i_track};
                [~,maxTime] = max(temp(:,5));
                t = temp(maxTime,4);
                if t == 0
                   t =1
                end
                for i_track = 1:60
                    t = i_track;
                currentFrame = squeeze (max(pa_img_all_type1(:,:,:,t), [],1));
                Allconverge1(1:84,1:84,i_track) = currentFrame;  %84
% % % % % % % % % % % % % % % % % % % %                   figure; imagesc(currentFrame,[]);
% % % % % % % % % % % % % % % % % % % %                       xlabel('y(mm)');ylabel('x(mm)');
% % % % % % % % % % % % % % % % % % % %                   colormap('hot');
%                  hold on;
%                  plot(temp(maxTime,3),temp(maxTime,2),'blue*');
%                  pos = [temp(maxTime,3),temp(maxTime,2)];
%                  text('Position',pos,'String',num2str(i_track),'FontSize',18,'color','green');
            end

            AllconvergeMAP1 = squeeze (max(Allconverge1,[],3));
% % % % % % % % % % % % % % % % %              figure; imagesc(y,x,AllconvergeMAP1);
% % % % % % % % % % % % % % % % %              title (num2str(pN))
% % % % % % % % % % % % % % % % %              xlabel('y(mm)');ylabel('x(mm)');colormap('hot');
% % % % % % % % % % % % % % % % %              ax = gca; % current axes
% % % % % % % % % % % % % % % % %              ax.LineWidth = 1.5;
% % % % % % % % % % % % % % % % %              ax.FontSize = 12;
% % % % % % % % % % % % % % % % %              ax.FontWeight = 'bold';

%             %% 2nd collapse
%             TF = isempty(track_tot2);
%             if TF == 1
%                AllconvergeMAP2 = zeros(63,63);
%                Allconverge2 = zeros(63,63,1);
%             else
%             if track_tot2{1} == 0
%                AllconvergeMAP2 = zeros(63,63);
%             else
%             maxInt = max(pa_img_all_type2,[],'all');
% 
%             for i_track = 1:size(track_tot2,1)
%                 temp = track_tot2{i_track};
%                 [~,maxTime] = max(temp(:,5));
%                 t = temp(maxTime,4);
%                 currentFrame = squeeze (max(pa_img_all_type2(:,:,:,t), [],1));
%                 Allconverge2(1:63,1:63,i_track) = currentFrame;
% % % % % % % % % % % % % % % % % % % % % % % %                 figure; imagesc(currentFrame,[800000 2000000]);
% % % % % % % % % % % % % % % % % % % % % % % %                   %     xlabel('y(mm)');ylabel('x(mm)');
% % % % % % % % % % % % % % % % % % % % % % % %                   colormap('hot');
% % % % % % % % % % % % % % % % % % % % % % % %                   hold on;
% % % % % % % % % % % % % % % % % % % % % % % %                   plot(temp(maxTime,3),temp(maxTime, 2),'blue*');
% % % % % % % % % % % % % % % % % % % % % % % %                   pos = [temp(maxTime,3),temp(maxTime,2)];
% % % % % % % % % % % % % % % % % % % % % % % %                  text('Position',pos,'String',num2str(i_track),'FontSize',18,'color','green');
%             end
% 
%             AllconvergeMAP2 = squeeze (max(Allconverge2,[],3));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %              figure; imagesc(y,x,AllconvergeMAP2);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %              title (num2str(pN),'- 2nd collapse')
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %             xlabel('y(mm)');ylabel('x(mm)');
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %             colormap('hot');
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %             ax = gca; % current axes
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %              ax.LineWidth = 1.5;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %             ax.FontSize = 12;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %             ax.FontWeight = 'bold';
%             end
%             end
            %% Store the two MAP into one MAP
            a = size(Allconverge1,3);
            Allconverge = zeros(84,84,a);    % x,y,
            Allconverge(:,:,1:size(Allconverge1,3)) = Allconverge1;
%             Allconverge(:,:,size(Allconverge1,3)+1:end) = Allconverge2;

            AllconvergeMAP = squeeze (max(Allconverge,[],3));
               
            fig = figure;
            imagesc(y,x,AllconvergeMAP);  %,[800000 2000000]
            title (num2str(pN),' - 1st&2nd collapse')
             xlabel('y(mm)');ylabel('x(mm)');
             colormap('hot');
             ax = gca; % current axes
             ax.LineWidth = 1.5;
             ax.FontSize = 12;
             ax.FontWeight = 'bold';

            WholePulseMAP(:,:,pN) = AllconvergeMAP;
            a = [];
% % % % % % % % % % % % % % % % % %             a(:,:,:,1:180) =  pa_img_all_type1;
                        a(:,:,:,1:tLen) =  pa_img_all_type1;
%             a(:,:,:,201:400) =  pa_img_all_type2;
            aa = squeeze (max(a,[],4));
            Whole3D100(:,:,:,pN)=aa;
            Whole3Dvolume = squeeze (max(Whole3D100,[],4));
            % Store MAP
 %%%%%%%%%%%           savefig(fig,[MAPSavePath,'\',testName,'_Pulse',num2str(pN),'_MAP.fig']);
            % Store 3D Volume
            save([MAPSavePath,'\',testName,'_Pulse',num2str(pN),'.mat'],'Whole3Dvolume');

    end
end
    end
end
%
Whole3D50 = squeeze (max(Whole3D100(:,:,:,20),[],4));

save([MAPSavePath,'\',testName,'_Pulse50','.mat'],'Whole3D50');

% Put all MAPs together
            WholePulseMAP1 = squeeze (max(Whole3Dvolume,[],1));
            figure; imagesc(y,x,WholePulseMAP1); %,[800000 2000000] y,x,
            title ('100 Pulse MAPs')
            xlabel('y(mm)');ylabel('x(mm)');
            colormap('hot');
            ax = gca; % current axes
            ax.LineWidth = 1.5;
            ax.FontSize = 12;
            ax.FontWeight = 'bold';





%%
% startIdx = startIdx(20,1);
v = VideoWriter('2D_T6_modifyFOV_1','MPEG-4');
v.FrameRate = 20;
open(v);
 for t = 60  % 1:length(tLen)
    fig = figure;
    temp_img = pa_img_all_type1(:,:,:,t)./max(pa_img_all_type1(:));
%     temp_img = 10*log10(abs(pa_img(:,:,:,t))./max(pa_img(:)));
    subplot(1,3,1);
    XYMAP = squeeze(max(temp_img,[],3));
    imagesc(imageArea.y_arr,imageArea.x_arr,XYMAP);
%     imagesc(XYMAP);
    colormap("hot");
%     colorbar;
    pbaspect([1 1 1]);
    caxis([0 0.1]);
%     caxis([-10 0]);
    subplot(1,3,2);
    XZMAP = squeeze(max(temp_img,[],2));
    imagesc(imageArea.z_arr,imageArea.x_arr,XZMAP);
%     imagesc(XZMAP);
    colormap("hot");
    pbaspect([1 1 1]);
%     titlRe(['St Idx = ',num2str(reconTime(t)+startIdx)]);
    caxis([0 0.1]);
%     caxis([-10 0]);
    subplot(1,3,3);
    YZMAP = squeeze(max(temp_img,[],1));
    imagesc(imageArea.z_arr,imageArea.y_arr,YZMAP);
%     imagesc(YZMAP);
    colormap("hot");
    pbaspect([1 1 1]);
%     colorbar;
%     caxis([-10 0]);
    caxis([0 0.1]);
    pause(0.2);
    frame = getframe(fig);
    writeVideo(v,frame);
end
close(v);

%%
for t = 1:201
    t
    temp_amp(t) = max(pa_img_all_type1(:,:,:,t),[],'all');
end
figure, plot(temp_amp)
[~,max_frame] = max(temp_amp);


    fig = figure(8);
    temp_img = abs(pa_img_all_type1(:,:,:,max_frame))./max(pa_img_all_type1(:));
%     temp_img = 10*log10(abs(pa_img(:,:,:,t))./max(pa_img(:)));
    subplot(1,3,1);
    XYMAP = squeeze(max(temp_img,[],3));
    imagesc(imageArea.y_arr,imageArea.x_arr,XYMAP);
%     imagesc(XYMAP);
    xlabel('y'); ylabel('x');
    colormap("hot");
%     colorbar;
    pbaspect([1 1 1]);
    caxis([0 1]);
%     caxis([-10 0]);
    subplot(1,3,2);
    XZMAP = squeeze(max(temp_img,[],2));
    imagesc(imageArea.z_arr,imageArea.x_arr,XZMAP);
%     imagesc(XZMAP);
    xlabel('z'); ylabel('x');    
    colormap("hot");
    pbaspect([1 1 1]);
    title(['St Idx = ',num2str(reconTime(t)+startIdx)]);
    caxis([0 1]);
%     caxis([-10 0]);
    subplot(1,3,3);
    YZMAP = squeeze(max(temp_img,[],1));
    imagesc(imageArea.z_arr,imageArea.y_arr,YZMAP);
%     imagesc(YZMAP);
    xlabel('z'); ylabel('y');
    colormap("hot");
    pbaspect([1 1 1]);
%     colorbar;
%     caxis([-10 0]);
    caxis([0 1]);

length = 0.24:0.24:20.16;
x_profile = XYMAP(:,round(end/2));
y_profile = XYMAP(round(end/2),:);
figure, hold on
plot(length,x_profile); plot(length,y_profile);
