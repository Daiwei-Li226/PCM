%% Check RF data from PCM
fdir = 'C:\Davia\20260204_Random\2_PCM\t3_70mm_y8mm\';
idatfilesname = '_0020.dat';
f = fopen([fdir,idatfilesname],'r');

iRFData = fread(f,'int16');
iRFData = reshape(iRFData,[],256);
fclose(f);
%  figure,plot([1:length(iRFData(:,120))]/15.625,iRFData(:,120));xlabel('time(us)')
   figure,plot(iRFData(:,120));xlabel('points')

%% Check RF data from V8_alternate_MultiFrame
fdir = 'E:\Davia\0712-2024_PlasticPlate\t15_500delay_water\';
idatfilesname = 'FIBER_0006.dat';
f = fopen([fdir,idatfilesname],'r');
iRFData = fread(f,'int16');
    parafile = dir(fullfile(fdir,'*.mat'));
    parafile = parafile.name;
    load([fdir,parafile]);
iRFData0 = reshape(iRFData,[length(iRFData)/256 RFSizeInfo(2) 1]);
fclose(f);

 TX1 = iRFData0(1:11776,:);
 TX2 = iRFData0(11777:12800,:);
 TX3 = iRFData0(12801:13824,:);
  TX31 = iRFData0(41473:42496,:);
   TX32 = iRFData0(42497:43520,:);
    TX33 = iRFData0(43521:44544,:);
      TX31_10 = iRFData0(327169:328192,:);
 figure,plot(TX1(:,1));xlabel('points');title('TX1')
 figure,plot(TX2(:,1));xlabel('points');title('TX2')
 figure,plot([1:length(TX1(:,1))]/15.625,TX1(:,1));xlabel('time(us)')
% % %  figure,plot(TX3(:,1));xlabel('points');title('TX3')
% % %   figure,plot(TX31(:,1));xlabel('points');title('TX31')
% % %    figure,plot(TX32(:,1));xlabel('points');title('TX32')
% % %       figure,plot(TX33(:,1));xlabel('points');title('TX33')
% % %          figure,plot(TX31_10(:,1));xlabel('points');title('31_10')
%% Check RF data from V5_alternate_code
fdir = 'E:\Davia\0531-2024_Bego\3_US_PCM_Alternate\s3_100delay_0.7\';
idatfilesname = 'FIBER_0090.dat';
f = fopen([fdir,idatfilesname],'r');
iRFData = fread(f,'int16');
    parafile = dir(fullfile(fdir,'*.mat'));
    parafile = parafile.name;
    load([fdir,parafile]);
iRFData0 = reshape(iRFData,[length(iRFData)/256 RFSizeInfo(2) 1]);
fclose(f);

 TX1 = iRFData0(1:11776,:);
 TX2 = iRFData0(43521:44544,:);
 figure,plot(TX1(:,1));xlabel('points');title('TX1')
 figure,plot([1:length(TX1)]/15.625,TX1(:,1));xlabel('time(us)');title('TX1')
 figure,plot(TX2(:,1));xlabel('points');title('TX2')

% % RFData = zeros(Receive(TXNum).endSample,256,P.numFrames);
% %             for iframe = 1:P.numFrames
% %                 RFData(:,:,iframe) = ...
% %                     iRFData0(Receive((iframe-1)*TXNum+1).startSample:...
% %                     Receive(iframe*TXNum).endSample,:);
% %             end
% %             
% % leng = size(RFData,1)./31;
% % TX1=RFData(1:leng,1,1);  % First frame
% % TX2=RFData(leng+1:18944,1,1);
% % TX3=RFData(leng*2+1:leng*3,1,1);
% % % figure,plot(TX2);xlabel('points')
% % figure,plot([1:leng]/15.625,TX1);xlabel('time(us)')
% % figure,plot([1:leng]/15.625,TX2);xlabel('time(us)')
% % figure,plot([1:leng]/15.625,TX3);xlabel('time(us)')

