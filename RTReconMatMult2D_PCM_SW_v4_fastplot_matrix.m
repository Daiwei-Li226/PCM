function RTReconMatMult2D_PCM_SW_v4_fastplot_matrix(ReceiveData)
% External function for realtime recon using matrix multiplication
% YT - 2022
% Add Phase correlation and save it to .mat file, DL,RY - 2025
% Change the way to show MAP to get fast

% persistent 
global frameCount GPU sensarr gmul_mat mult_mat area1 temperature_C delay_t travelingTime reconTime startIdx PCM_img max_frame temp_img XYMAP XZMAP CoordinatesSavedir map_buffer hPlot hImg1 hImg2 firstMAP 

tic

T = temperature_C;
sos = 1402.4+5.01*T-0.055*T^2+0.00022*T^3;
sos = sos*1e-3;

% Save an empty coordiantes file for frame #1 
if isempty(frameCount)
a = [0,0,0,0];
global a
CoordinatesSavedir = 'C:\Davia\';   %'C:\Users\PI-Lab\Desktop\PyDobot\';
CoordinatesSavename = [CoordinatesSavedir,'coordinates.mat'];
save(CoordinatesSavename,'a')
else 
end

% global laser_choice 
RFData = ReceiveData(:);
RFData = double(reshape(RFData,[round(size(RFData,1)/256),256]));
Receive = evalin('base','Receive');

% determine startIdx
fs = Receive(1).decimSampleRate;
reconStep = 950;
travelingTime = 40/1.5*fs;
%  peak = find(RFData(:,120) > 250, 1);
 [peak,Index] = max(RFData(:,120));
%  Index = find(RFData(:,120) > 4000, 1);
if startIdx > travelingTime
startIdx = round(Index - travelingTime);  %-100;300 steps 
else
    startIdx =1;
end
% if isempty(peak)
%    startIdx = 7500;
% disp(['NO Peak!'])
% end

step = 5; % 0.064 us per step
reconTime = 0:step:reconStep;

rfdata = RFData(startIdx:startIdx+1000,:,:);   % only 1001 points for rfdata

% tic
% create rf class for rfdata
pulse_index = round(delay_t*fs);
rf_class = RFclass2D_PCM(rfdata,pulse_index,fs,sos); 

if isempty(frameCount)
mult_mat = IndexMatrix2D_PCM(sensarr,area1,rf_class,"matrix"); 
    if GPU
       gmul_mat = mult_mat.Send2GPU;
    end
end
% elapsedTime = toc; 
% disp(['Matrix Calculation time: ', num2str(elapsedTime), ' s'])

% Recon
% tic
for iframe = 1:1
    for t = 1:length(reconTime)
            ishift = step*(t-1)+startIdx;
            rfdata_recon = squeeze(rfdata(:,:,iframe));
            rfdata_recon = circshift(rfdata_recon,[-ishift,0]);
            rf_class = RFclass2D(rfdata_recon,0,fs,sos);

            temp_rec = DAS_mmult(rf_class,gmul_mat);
            temp_rec = reshape(temp_rec,length(area1.x_arr),...
                length(area1.y_arr),length(area1.z_arr));
            PCM_img(:,:,:,t) = gather(temp_rec);

    end
end

% elapsedTime = toc; 
% disp(['Find StartIdx + Recon time: ', num2str(elapsedTime), ' s'])

% Show MAP
% tic;

if isempty(frameCount)
fig1 = figure(2); clf; set(fig1, 'Position', [100, 500, 400, 300]);
fig2 = figure(3); clf; set(fig2, 'Position', [600, 500, 400, 300]);
fig3 = figure(4); clf; set(fig3, 'Position', [1100, 500, 400, 300]);
y = zeros(1, size(reconTime,2));   
img_data1 = zeros(size(area1.x_arr,1), size(area1.y_arr,1));  
img_data2 = zeros(size(area1.x_arr,1), size(area1.z_arr,1));  
figure(2); hPlot = plot(y);
figure(3);   % hImg1 = imshow(img_data1, [],'Colormap', hot, 'InitialMagnification', 'fit');axis on;xlabel('x');ylabel('y');title('Top MAP (XY)');
hImg1 = imagesc(area1.x_arr, area1.y_arr, img_data1);  axis image;set(gca, 'YDir', 'normal'); % set y axis upwards
colormap(hot);xlabel('x(mm)'); ylabel('y(mm)');title('Top MAP (XY)');
figure(4);  % hImg2 = imshow(img_data2, [],'Colormap', hot, 'InitialMagnification', 'fit');axis on; xlabel('z');ylabel('x');title('Side MAP (XZ)');
hImg2 = imagesc(area1.z_arr, area1.x_arr, img_data2);  axis image;set(gca, 'YDir', 'normal');
colormap(hot);xlabel('z(mm)'); ylabel('x(mm)');title('Side MAP (XZ)');

frameCount =1;
else 
end

% Find most converge time point
PCM_img = abs(PCM_img);
for t = 1:length(reconTime)
    temp_amp(t) = max(PCM_img(:,:,:,t),[],'all');
end
[~,max_frame] = max(temp_amp);

% Get MAP   
temp_img = abs(PCM_img(:,:,:,max_frame))./max(PCM_img(:));
XYMAP = squeeze(max(temp_img,[],3));
XZMAP = squeeze(max(temp_img,[],2));

% Save map_buffer as rolling (only 2 frames each time)
if size(map_buffer,3) >1
    map_buffer(:,:,1) = map_buffer(:,:,2);    % Old MAP
    map_buffer(:,:,2) = XZMAP;                % New MAP
else
    map_buffer = cat(3,map_buffer,XZMAP);
    if frameCount ==1 
        firstMAP = XZMAP;
    end
end
% Phase Correlation - Generate dx , dy and save to coordinates.mat file
bb = (frameCount-1)./10;
if mod(bb,1)==0
    if frameCount > 1
[dy,dx] = phaseCorr(map_buffer(:,:,1),map_buffer(:,:,2));
a = [0,0,0,0];
a(1) = dx;a(2) = dy; a(4)=frameCount; 
CoordinatesSavename = [CoordinatesSavedir,'coordinates.mat'];
save(CoordinatesSavename,'a')
disp(['saved new coordinates:', num2str(a)])
    end 
end

elapsedTime = toc; 
disp(['Get coordinates time: ', num2str(elapsedTime), ' s'])

% Show images
set(hPlot,'YData',temp_amp);
set(hImg1, 'CData', XYMAP);%XYMAP is Side View
set(hImg2, 'CData', XZMAP);%XZMAP is Top View

drawnow;  % Ensure MATLAB updates the figures in real-time

% elapsedTime = toc; 
% disp(['Only Display time: ', num2str(elapsedTime), ' s'])



    function [dy,dx] = phaseCorr(ref,mov)
        ref_fft = fft2(ref);
        mov_fft = fft2(mov);
        cor = ifft2(((ref_fft.*conj(mov_fft))./abs(ref_fft.*conj(mov_fft))));
        [~,idx] = max(cor(:));
        [dy,dx] = ind2sub(size(cor),idx);
        dy = dy-1;
        dx = dx-1;
        if dy>size(ref,1)/2
            dy = dy - size(ref,1);
        end
        if dx> size(ref,2)/2
            dx = dx - size(ref,2);
        end
    end



frameCount = frameCount + 1




end