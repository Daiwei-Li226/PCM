%%
% % % % % % % % % % % % % % % clear all; clc; %close all; 
%% reaing RF data
% % % % % % % % % % % % % filePath = 'D:\Davia\20250423\20250423_US\';
% % % % % % % % % % % % % fileName = 'Pulseecho-900us_US.mat';
% matName = dir([filePath '*.mat']);
% load([filePath matName.name]);
% % % % % % load([filePath ...
% % % % % %     fileName]);

RFData1 = RcvData{1};
% for iTX = 1:128
for iTX = 1:256
    RFData(:,:,iTX,:) = double(RFData1(Receive(iTX).startSample:Receive(iTX).endSample,:,:));
end

% checking RF data
fs = Receive(1).decimSampleRate;
T = 21.6;
sos = (1402.4+5.01*T-0.055*T^2+0.00022*T^3)*1e-3;
timeFs = [1:size(RFData,1)]*(1/fs)*1.5/2;
%%
figure
hold on
for iTX = 1:256
%     plot(timeFs,RFData(:,iTX,iTX,1))
    plot(RFData(:,iTX,iTX,1))
end
xlim([0 200]);
ylim([-3000 3000]);
%%
figure
plot(RFData(:,2,2,1))
ylabel('Amplitude')
xlabel('Sample')
%%
figure
hold on
% for iTX = 1:128
for iTX = 1:256
    stIdx = 1300; edIdx =1500;
%     stIdx = 450; edIdx = 550;
%     stIdx = 900; edIdx = 1200;
    rf_tmp = RFData(stIdx:edIdx,iTX,iTX,1);
    pkp_val(iTX) = max(rf_tmp(:)) - min(rf_tmp(:));
%     plot(timeFs,RFData(:,iTX,iTX,1))
end
% pkp_val = pkp_val / max(pkp_val) * 400;
% pkp_val(121:128) = [];
% pkp_val(51:70) = [];

% pkp_val(101:128) = [];

% pkp_val(121:128) = [];
% pkp_val(91:95) = [];
% pkp_val(60:65) = [];
% pkp_val(29:34) = [];
% pkp_val(:,1:3) = [];

plot(pkp_val,'.');
hold on
% yline(mean(pkp_val))
yline(20,'--')
box on;
ylabel('pk-to-pk (unit)');
xlabel('Channel');
xlim([0 128]);
%%
elem_sens = pkp_val.^0.5;
figure;
plot(elem_sens);
hold on;
box on
yline(4.5,'--')
yline(mean(elem_sens),'r-')
ylabel('Sens (unit)');
xlabel('Channel');
% xlim([0 128]);
xlim([0 256]);
%%
numel(find(elem_sens<20))
%%
% figure
% hold on
% % for iTX = 1:100
%     plot(RFData(:,41,41,1))
%     plot(RFData(:,43,43,1))
%     plot(RFData(:,57,57,1))
%     plot(RFData(:,26,26,1))
% % end
% xlim([50 150]);
% ylim([-200 200]);
%%
