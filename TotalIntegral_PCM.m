% PCM E_total and E_peak based on Spatio-Temporal Integration
% 11/23/2025
% Davia 
clear all; close all;

data_dir = 'C:\Users\daiwe\Dropbox\CavitationMapping_2024_Davia_Chaorui\20250612_Pig\stonw2';
MATfile = false;
fs = 19.23e6; 
% % % % % % % % % % % % % % T = 22.0;
% % % % % % % % % % % % % % sos = (1402.4+5.01*T-0.055*T^2+0.00022*T^3); % m/s
% Target position (mm)
x0 = 15.72; 
y0 =12.38; 
z0 = 38.10; 
K_elements = 100; % find 64 ele cloest to the target (8x8)

%% Element position
num_rows = 16;num_cols = 16;
pitch_mm = 2.0; % mm

load('C:\Users\daiwe\Dropbox\CavitationMapping_2024_Davia_Chaorui\20250807_5MHzV2_PCMsensitivtity\4_PCM_RF_modelbase\x_trans.mat')
load('C:\Users\daiwe\Dropbox\CavitationMapping_2024_Davia_Chaorui\20250807_5MHzV2_PCMsensitivtity\4_PCM_RF_modelbase\y_trans.mat')
load('C:\Users\daiwe\Dropbox\CavitationMapping_2024_Davia_Chaorui\20250807_5MHzV2_PCMsensitivtity\4_PCM_RF_modelbase\z_trans.mat')
%figure, hold on
for ielem = 1:256  %size(x_trans,2)
    scatter(x_trans(ielem),y_trans(ielem));
    text(x_trans(ielem),y_trans(ielem),num2str(ielem));
end
% Get ROI
% distance 
distances = sqrt((x_trans - x0).^2 + (y_trans - y0).^2 + (z_trans - z0).^2);
% ROI: find shortest distance elements
[~, sorted_indices] = sort(distances);
roi_indices_spatial = sorted_indices(1:K_elements); %Final ROI
fprintf('ROI filtered. Chose the cloest %d elements to (%.1f, %.1f, %.1f) mm \n',  K_elements, x0, y0, z0);
%% Batch processing function
%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % function [Energy_Results, Peak_Results, TOI_Used] = process_cavitation_data_focused(data_folder, roi_indices, fs)
    data_folder = data_dir;
    roi_indices = roi_indices_spatial;
%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % .mat or .dat
    if MATfile
    file_list = dir(fullfile(data_dir, '*.mat'));
    else
    file_list = dir(fullfile(data_dir, '*.dat'));
    end
    num_files = length(file_list);

    % Initianate results matrix (E_total and E_peak)
    Energy_Results = zeros(num_files, 1);
    Peak_Results = zeros(num_files, 1);
    TOI_Used = zeros(num_files, 2); % [start_sample, end_sample]
    
    sampling_period_s = 1 / fs;
    W_samples = 500; % TOI window, pts, around 26 us
tic
    for i = 1:num_files
        file_path = fullfile(data_folder, file_list(i).name);        
        try
            
            if MATfile
               data_struct = load(file_path);
               % !!! Make sure data format is NumSamples x NumChannels) 
               field_names = fieldnames(data_struct);
               RF_matrix = data_struct.(field_names{1}); 
            else
               f = fopen(file_path,'r');
               RF_matrix = fread(f,'int16');
               RF_matrix = reshape(RF_matrix,[],256);
               fclose(f);
            end
    
            % Find RF peak
            roi_data_full = RF_matrix(:, roi_indices(1)); 
            [max_val_cols, max_idx_cols] = max(abs(roi_data_full));
            [max_peak, i_channel_roi] = max(max_val_cols); % V_Peak
            Peak_Results(i) = max_peak;
            % Find where is RF peak
            N_peak = max_idx_cols(i_channel_roi); 
            
            % Plot RF data
%            figure;plot(RF_matrix(:,roi_indices(1)))

            % Define TOI
            N_buffer_pre = 100; % Get 100 pts buffer to capture peak
            start_sample = max(1, N_peak - N_buffer_pre); 
            end_sample = min(length(roi_data_full), start_sample + W_samples - 1);
            TOI_Used(i, :) = [start_sample, end_sample]; % Save TOI
            
            % Energy calculation (spatial + temporal integration)
            % Extract data within ROI & TOI 
            roi_data_subset = RF_matrix(start_sample:end_sample, roi_indices); 
            V_squared = roi_data_subset.^2; % Calculate V^2
            energy_per_channel = sum(V_squared, 1) *1;% sampling_period_s; % V^2*s Energy calculation (time-domain integration)
            Total_Acoustic_Energy = sum(energy_per_channel);% Total vacuum radiation energy (space integration)
            Energy_Results(i) = Total_Acoustic_Energy;
            
            fprintf('Processed file %d/%d: %s. E_Total: %.3e V^2*s, V_Peak: %.3f V, TOI: %d-%d\n', ...
                    i, num_files, file_list(i).name, Total_Acoustic_Energy, max_peak, start_sample, end_sample);
        
        catch ME   % If process failed, give NaN
            warning('Error processing file %s: %s', file_list(i).name, ME.message);
            Energy_Results(i) = NaN;
            Peak_Results(i) = NaN;
        end
    end
    
    disp('PCM E_total and E_peak based on Spatio-Temporal Integration');
% end
toc
% How to use the function:
% data_dir = 'C:\MyCavitationData'; 
% [E_results, P_results, TOI_log] = process_cavitation_data_focused(data_dir, roi_indices_spatial, fs, num_samples_total);
% disp(E_results);
% disp(P_results);

%% Plot results
figure;plot(Energy_Results);title('Total Engery along pN')
figure;plot(Peak_Results);title('Peak Engery along pN')