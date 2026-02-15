function Integral = saveRcvDataRT_feedback_2D_PCM(ReceiveData)
global Energy_Results
global Peak_Results
global TOI_Used

if isempty(Energy_Results)
   Energy_Results = ones(10,1);
   Peak_Results = ones(10,1);
   TOI_Used = ones(10,1); % [start_sample, end_sample]
end

tic
global RCVdataindex;
global fdir;
global dataname;

% global fs;
% global Trans; %.ElementPosMm
global x_trans;
global y_trans;
global z_trans;
if isempty(x_trans)
load(['C:\Davia\20251125_RT_feedback\','x_trans.mat']);
load(['C:\Davia\20251125_RT_feedback\','z_trans.mat']);
load(['C:\Davia\20251125_RT_feedback\','y_trans.mat']);
end

diginum = 4;
if isempty(RCVdataindex)
    RCVdataindex = 1; 
end

%% Total Integral PCM feedback
x = 0; 
y = 0; 
z = 40; 
K_elements = 100; % find 64 ele cloest to the target (8x8)

distances = sqrt((x_trans - x).^2 + (y_trans - y).^2 + (z_trans - z).^2);
% ROI: find shortest distance elements
[~, sorted_indices] = sort(distances);
roi_indices = sorted_indices(1:K_elements); %Final ROI
if RCVdataindex < 3
fprintf('ROI filtered. Chose the cloest %d elements to (%.1f, %.1f, %.1f) mm \n',  K_elements, x, y, z);
end

W_samples = 500; % TOI window, pts, around 26 us

RF_matrix = ReceiveData(:,:);
            % Find RF peak
            roi_data_full = RF_matrix(:, roi_indices(1)); 
            [max_val_cols, max_idx_cols] = max(abs(roi_data_full));
            [max_peak, i_channel_roi] = max(max_val_cols); % V_Peak
            Peak_Results(RCVdataindex) = max_peak;
            % Find where is RF peak
            N_peak = max_idx_cols(i_channel_roi); 

            % Define TOI
            N_buffer_pre = 100; % Get 100 pts buffer to capture peak
            start_sample = max(1, N_peak - N_buffer_pre); 
            end_sample = min(length(roi_data_full), start_sample + W_samples - 1);
            TOI_Used(RCVdataindex) = start_sample; % Save TOI
            
            % Energy calculation (spatial + temporal integration)
            % Extract data within ROI & TOI 
            roi_data_subset = RF_matrix(start_sample:end_sample, roi_indices); 
            V_squared = roi_data_subset.^2; % Calculate V^2
            energy_per_channel = sum(V_squared, 1) *1;% sampling_period_s; % V^2*s Energy calculation (time-domain integration)
            Total_Acoustic_Energy = sum(energy_per_channel);% Total vacuum radiation energy (space integration)
            Energy_Results(RCVdataindex) = Total_Acoustic_Energy;

            fprintf('Processed pN = %d. E_Total: %.3e V^2*s, V_Peak: %.3f V, TOI: %d', ...
                    RCVdataindex, Total_Acoustic_Energy, max_peak, start_sample);
            if  Peak_Results(RCVdataindex)<50  % If process failed, give NaN
            warning('Error processing pN = %i.', RCVdataindex);
%             Energy_Results(i) = NaN;
%             Peak_Results(i) = NaN;
            end
clear roi_data_full

figure(10);subplot(2,1,1);bar(Energy_Results)
subplot(2,1,2);stem(Peak_Results)

Integral(:,1) = Energy_Results; Integral(:,2) = Peak_Results; Integral(:,3) = TOI_Used;
%'C:\Users\PI-Lab\Desktop\PyDobot\';
Savename = [fdir,'\','E_total.mat'];
save(Savename,'Integral')

%%

PARcvData = ReceiveData(:);

f0 = sprintf(['%0',num2str(diginum),'d'],RCVdataindex);
f1 = fopen([fdir,'\',dataname,'_',num2str(f0),'.dat'],'w');
fwrite(f1,PARcvData,'int16');
fclose(f1);
PARcvData = [];

toc

RCVdataindex = RCVdataindex + 1;