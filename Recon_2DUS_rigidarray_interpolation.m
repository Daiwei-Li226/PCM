%% Vectorized 2D US recon
% Edit by Nanchao
% 2023/09/19
% Simplify the recon script
% Updated by Nanchao
% 2023/09/26
% Fit for the softarray
%% Recon_PA_2D
clear, clc
g = gpuDevice(1);
reset(g);
fileName = 't5_z100_x5__US.mat';
dirName = 'C:\Users\daiwe\Dropbox\CavitationMapping_2024_Davia_Chaorui\20250807_5MHzV2_PCMsensitivtity\2_US\';
load([dirName fileName]);

%%
Res = 0.16; %0.2
cAngle = 80*pi/180;  %
x_range = [-10,10];
y_range = [-10,10];
z_range = [95,105];
x_img = x_range(1):Res:x_range(2);
y_img = y_range(1):Res:y_range(2);
z_img = z_range(1):Res:z_range(2);
%% reading sensor location

% Freq = 1;
% 
% Trans2561.name = 'qcr2561';
% Trans2561.units = 'wavelengths';
% Trans2561 =  computeTrans_256_2D_Martix(Trans2561);    % L7-4 transducer is 'known' transducer so we can use computeTrans.
% 
% Trans2562.name = 'qcr2562';
% Trans2562.units = 'wavelengths'; 
% Trans2562 = computeTrans_256_2D_Martix(Trans2562);
% 
% waveLength = (1.5e3/1000)/Freq;
% Trans.ElementPos = [Trans2561.ElementPos; Trans2562.ElementPos];
% Trans.ElementPosMm = Trans.ElementPos * waveLength;
% 
% x_trans = Trans.ElementPosMm(:,1);
% y_trans = Trans.ElementPosMm(:,2);
% z_trans = Trans.ElementPosMm(:,3);

RecoFreq = 5.6; 
temperature_c= 21;
sos= (1402.4+5.01*temperature_c-0.055*temperature_c^2+0.00022*temperature_c^3);

Trans2561.name = 'qcr2561';
Trans2561.units = 'wavelengths';
Trans2561 =  computeTrans_256_2D_Martix(Trans2561,sos);    % L7-4 transducer is 'known' transducer so we can use computeTrans.

Trans2562.name = 'qcr2562';
Trans2562.units = 'wavelengths'; 
Trans2562 = computeTrans_256_2D_Martix(Trans2562,sos);

Trans.frequency=Trans2561.frequency;
waveLength = (sos/1000)/RecoFreq;
Trans.ElementPosMm = [Trans2561.ElementPosMm; Trans2562.ElementPosMm];
Trans.ElementPosWL = Trans.ElementPosMm /waveLength;

x_trans = Trans.ElementPosMm(:,1);
y_trans = Trans.ElementPosMm(:,2);
z_trans = Trans.ElementPosMm(:,3); 

% x_trans = Trans.ElementPos(:,1);
% y_trans = Trans.ElementPos(:,2);
% z_trans = Trans.ElementPos(:,3);

% % % % % % % % % % % % % % % % % % % % % % % % % load('C:\Users\daiwe\Dropbox\CavitationMapping_2024_Davia_Chaorui\Trans_and_Normals_fitted.mat');
figure, scatter(x_trans,y_trans);

    % % x_trans = x_trans(:);y_trans = y_trans(:);z_trans = z_trans(:);
    shapeflag = 1;
% catch ME
%     if true
% %         Trans_location_Name = dir([dirName 'Base_real_position.mat']);
%         Trans_location_Name = 'Base_real_position.mat';
%         load([dirName Trans_location_Name]);
%         x_trans = reorder_x;
%         y_trans = reorder_y;
%         z_trans = -(reorder_z - max(reorder_z));
%     end
% end

% get RF data
clear RFData
Trans.numelements = 256;
RFData = RcvData{1};
RFData = RFData(1:Receive(end).endSample,:,:);
% RFData = mean(RFData,3);
trans_num = 256;
% RFData(:,101:128) = [];


% [nx,ny,nz] = surfnorm(reshape(x_trans,[10 10]),reshape(y_trans,[10 10]),reshape(z_trans,[10 10]));
% nx = nx(:); ny = ny(:); nz = nz(:);
%%  Interpolation if needed 

% Average RF data
RFData = mean(RFData,3);

int = false;
int_dis = 4;    % 7 will be too slow

if int 
[num_samples, num_elem, num_frames] = size(RFData);
x_orig = x_trans(:);
y_orig = y_trans(:);
z_orig = z_trans(:);

% === do interpolation pair once (ielem, ipair) ===
interpolation_pairs = [];

for ielem = 1:num_elem
    for ipair = ielem+1:num_elem  % symmetrical
        dis = sqrt((x_orig(ielem) - x_orig(ipair))^2 + ...
                   (y_orig(ielem) - y_orig(ipair))^2 + ...
                   (z_orig(ielem) - z_orig(ipair))^2);
        if dis < int_dis
            interpolation_pairs = [interpolation_pairs; ielem, ipair];
        end
    end
end

% === Transducer position after interpolation ===
x_new = x_orig;
y_new = y_orig;
z_new = z_orig;

for i = 1:size(interpolation_pairs,1)
    ielem = interpolation_pairs(i,1);
    ipair = interpolation_pairs(i,2);

    x_new = [x_new; mean([x_orig(ielem), x_orig(ipair)])];
    y_new = [y_new; mean([y_orig(ielem), y_orig(ipair)])];
    z_new = [z_new; mean([z_orig(ielem), z_orig(ipair)])];
end

new_num_elem = num_elem + size(interpolation_pairs,1);
RFData_interp = zeros(num_samples, new_num_elem, num_frames, 'like', RFData);  

% === Interpolation for all frames ===
for iFrame = 1:num_frames
    RF_slice = double(RFData(:,:,iFrame)); 

    % save original dataset
    RF_interp = RF_slice;

    % Interpolation
    for i = 1:size(interpolation_pairs,1)
        ielem = interpolation_pairs(i,1);
        ipair = interpolation_pairs(i,2);

        new_RF = mean([RF_slice(:,ielem), RF_slice(:,ipair)], 2);
        RF_interp(:, end+1) = new_RF;
    end

    RFData_interp(:,:,iFrame) = RF_interp;
    iFrame
end

% Renew transducer position
x_trans = x_new;
y_trans = y_new;
z_trans = z_new;
else
end
%% 
figure, hold on
for ielem = 1:length(x_trans)
    scatter(x_trans(ielem),y_trans(ielem));
    text(x_trans(ielem),y_trans(ielem),num2str(ielem));
end
figure
scatter3(x_trans,y_trans,-z_trans);
% hold on
% scatter3(x_out,y_out,(z_out - max(z_out)));
xlabel('x (mm)');ylabel('x (mm)');zlabel('x (mm)');
% hold on
% for ielem = 1:numel(x_trans)
%     text(x_trans(ielem),y_trans(ielem),num2str(ielem));
% end

mat = [x_trans(:), y_trans(:), z_trans(:)];
ptCloud = pointCloud(mat);
normals = pcnormals(ptCloud,3);
% ptCloud.Normal = normals;
nx = normals(:,1);
ny = normals(:,2);
nz = normals(:,3);



%%  Recon
%
    clear transSenAll rRCV
    Trans.numelements = size(x_trans,1);
    % init
    T = mediumTemp;                                                  
    sos = (1402.4+5.01*T-0.055*T^2+0.00022*T^3);
    fs = Receive(1).decimSampleRate*1e6;
    elemWidth = 2; % mm
    elemWidthwl = elemWidth/(sos/1e3/Trans.frequency);
    wl = sos/Trans.frequency/1e3;
    if ~exist('transSenAll','var')
        transSenAll = calc_angular_sens(x_trans,y_trans,z_trans,x_img,...
            y_img,z_img,elemWidthwl,Trans,nx,ny,nz);
        transSenAll = gpuArray(single(transSenAll));
    end
    if ~exist('rRCV','var')
        rRCV = calc_deley_mat(x_trans,y_trans,z_trans,x_img,...
            y_img,z_img,Trans);
    end
    TXNum = size(TX,2);
    RFshift = TW.peak/(Trans.frequency)*Receive(1).decimSampleRate;
%      RFshift = 2.4720/(Trans.frequency)*Receive(1).decimSampleRate;      % Sav TW.peak   
    % RFshift = -5;
    % RFshift = -10;
    RFshift = RFshift + 4e3/sos/(1/Receive(1).decimSampleRate) + 2e3/sos/(1/Receive(1).decimSampleRate);
  
%%%%%%%%%%%%%%%%%
if int
  RFData  = RFData_interp;
end

        % recon
        framenum = size(RFData,3);
        us_rec = single(zeros(length(x_img),length(y_img),length(z_img),framenum));
        tic

       % Parameters
        num_coefficients = 256; % Number of coefficients to generate
        mean_val = 1.0;         % Mean sensitivity value
        std_dev = 0.05;         % Standard deviation (fluctuation level)
        seed = 42;              % Seed for reproducibility
        
        % Generate and visualize the coefficients
        weight = generateSensitivityCoefficients(num_coefficients, mean_val, std_dev, seed);
        

        for itx = 1:TXNum
            itx
            alineStart = Receive(itx).startSample;
            alineEnd = Receive(itx).endSample;
            iRF = weight(itx) * single(RFData(alineStart:alineEnd,:,:));
            iRF = imtranslate(iRF,[0 -RFshift],'FillValues',0);

            iRF = gpuArray(single(hilbert(iRF)));
            txelem = find(TX(itx).Apod>0);
%             txelem = itx;
            idxAll = gpuArray(single(round((repmat(rRCV(:,:,:,txelem),[1 1 1 Trans.numelements])...
                + rRCV - 2*P.startDepth*wl)/1e3/sos*fs)));
            for iframe = 1:framenum
                us_bmode = gpuArray(single(zeros(size(us_rec(:,:,:,1)))));
                for ielem = 1:Trans.numelements
                    [us_bmode,us_bmode1] = GPUReconLoop(iRF,idxAll,iframe,ielem,us_bmode,transSenAll);
                end
                us_rec(:,:,:,iframe) = us_rec(:,:,:,iframe) + single(gather(us_bmode));
            end
        end
        toc

        CR = 30;
%         CR = 5;
        % display
        xymap = squeeze(max(abs(mean(us_rec(:,:,:,:),4)),[],3))';
        xzmap = squeeze(max(abs(mean(us_rec(:,:,:,:),4)),[],1))';
        yzmap = squeeze(max(abs(mean(us_rec(:,:,:,:),4)),[],2))';
%         xymap = squeeze(max(mean(us_rec(:,:,:,:),4)),[],3)';
%         xzmap = squeeze(max(mean(us_rec(:,:,:,:),4)),[],1)';
%         yzmap = squeeze(max(mean(us_rec(:,:,:,:),4)),[],2)';
        fig = figure();
         subplot(1,3,1),imagesc(x_img,y_img,20*log10(xymap/2.2245e+04))       % max(xymap(:))
%         subplot(1,3,1),imagesc(x_img,y_img,20*log10(xymap/6.2736e+04))
        colormap gray,caxis([-CR 0]) 
        pbaspect([length(x_img)/length(y_img) 1 1]),xlabel('x [mm]'),ylabel('y [mm]')
        subplot(1,3,2),imagesc(x_img,z_img,20*log10(xzmap/max(xzmap(:))))
%         subplot(1,3,2),imagesc(x_img,z_img,20*log10(xzmap/6.2736e+04))
        colormap gray,caxis([-CR 0]) 
        pbaspect([length(y_img)/length(z_img) 1 1]),xlabel('y [mm]'),ylabel('z [mm]')
        subplot(1,3,3),imagesc(y_img,z_img,20*log10(yzmap/max(yzmap(:))))
%         subplot(1,3,3),imagesc(y_img,z_img,20*log10(yzmap/6.2736e+04))
        colormap gray,caxis([-CR 0]) 
        pbaspect([length(x_img)/length(z_img) 1 1]),xlabel('x [mm]'),ylabel('z [mm]')
        drawnow

        saving_path = [dirName,'Results',filesep];
        mkdir(saving_path)
        
        if int == true
            fileName=erase(fileName,'.mat');
            exportgraphics(fig,[saving_path,fileName,'_recon.tif']); 
                  save([saving_path,fileName,'_recon_4_intern_c.mat'],'us_rec');
        else
            exportgraphics(fig,[saving_path,fileName,'.tif']); 
        end
        
        
        
        
        
        
        
        
        %%
% % [ius_svd] = SVD_filter(us_rec,20);
% 
% ius_rec = us_rec;
% cutoff = 20;
% 
% 
% fprintf('Clutter filtering ... \n');
% tcluter = tic;
% [row,col,page,ensemble] = size(ius_rec);
% S = reshape(ius_rec,[row*col*page,ensemble]);
%     % direct SVD
%     [U,D,V] = svd(S,0);
%     % Dtemp = log10(diag(D));
%     Dtemp = 10*log10(diag(D)/max(D(:)+1e-10));   % need to prevent the case of D = 0
%     % mean doppler frequency
%     fs = 1/(36*150e-6);   % sampling rate
%     Nv = size(V,1);
%     for k = 1:Nv
%         r(k) = sum(conj(V(1:Nv-1,k)).*V(2:Nv,k)); % calculate the self correlation
%     end
%     phi = angle(r);
%     phi = abs(phi)/pi*fs;
% 
%     figure(),subplot(2,1,1)
%     plot(Dtemp), xlabel('Singular value order'), ylabel('Mag(dB)')
%     subplot(2,1,2)
%     plot(phi), xlabel('Singular value order'), ylabel('Frequency(Hz)')
% 
%     auto_cor = corr(abs(U));
%     figure()
%     imagesc(corr(abs(U)));
%     pbaspect([1 1 1])
%     colormap('turbo');
%     colorbar
% 
%     % cor_dif = auto_cor(2:end,:) - auto_cor(1:end-1,:);
%     % cor_dif_sum = sum(cor_dif,2);
%     % figure, plot(cor_dif_sum);
%     
%     cutoff_low = cutoff; cutoff_high = 300;
%     D0 = D;
%     D0(:,1:cutoff_low) = 0;
%     D0(:,cutoff_high:end) = 0;
%     S = U*D0*V;
%     ius_svd = single(abs(reshape(S,[row,col,page,ensemble])));
% 
% toc(tcluter)
% 
% 

% end
% %%
%         xymap = squeeze(max(abs(mean(ius_svd(:,:,:,:),4)),[],3))';
%         xzmap = squeeze(max(abs(mean(ius_svd(:,:,:,:),4)),[],1))';
%         yzmap = squeeze(max(abs(mean(ius_svd(:,:,:,:),4)),[],2))';
% %         xymap = squeeze(max(mean(us_rec(:,:,:,:),4)),[],3)';
% %         xzmap = squeeze(max(mean(us_rec(:,:,:,:),4)),[],1)';
% %         yzmap = squeeze(max(mean(us_rec(:,:,:,:),4)),[],2)';
%         fig = figure();
%         subplot(1,3,1),imagesc(x_img,y_img,20*log10(xymap/max(xymap(:))))
%         colormap gray,caxis([-CR 0]) 
%         pbaspect([length(x_img)/length(y_img) 1 1]),xlabel('x [mm]'),ylabel('y [mm]')
%         subplot(1,3,2),imagesc(x_img,z_img,20*log10(xzmap/max(xzmap(:))))
%         colormap gray,caxis([-CR 0]) 
%         pbaspect([length(y_img)/length(z_img) 1 1]),xlabel('y [mm]'),ylabel('z [mm]')
%         subplot(1,3,3),imagesc(y_img,z_img,20*log10(yzmap/max(yzmap(:))))
%         colormap gray,caxis([-CR 0]) 
%         pbaspect([length(x_img)/length(z_img) 1 1]),xlabel('x [mm]'),ylabel('z [mm]')
%         drawnow
%         exportgraphics(fig,[saving_path,fileName,'after-SVD.tif']); 

function rRCV = calc_deley_mat(x_trans,y_trans,z_trans,x_img,y_img,z_img,Trans)
x_img0 = repmat(x_img',[1 length(y_img) length(z_img)]);
y_img0 = repmat(y_img,[length(x_img) 1 length(z_img)]);
z_img0 = repmat(z_img',[1 length(x_img) length(y_img)]);
z_img0 = permute(z_img0,[2 3 1]);
% if ~exist('rRCV','var')
    rRCV = zeros(length(x_img),length(y_img),length(z_img),Trans.numelements);
    % transSenAll = zeros(size(rRCV));
    elemAngleAll = zeros(1,Trans.numelements);
    tdend = 0;
    wb = waitbar(0,'Calculating RCV dist. map ...');
    for ielem = 1:Trans.numelements
        tdstart = tic;
        elemAngleAll(ielem) = abs(atan(sqrt(x_trans(ielem).^2 + y_trans(ielem).^2)/z_trans(ielem)));
        rRCV(:,:,:,ielem) = sqrt((x_img0 - x_trans(ielem)).^2 + ...
            (y_img0 - y_trans(ielem)).^2 + ...
            (z_img0 - z_trans(ielem)).^2);
        tdend = tdend + toc(tdstart);
        tdavg = tdend/ielem;
        tdrem = (Trans.numelements - ielem)*tdavg;
        waitbar(ielem/Trans.numelements,wb,...
            [sprintf('%12.1f',tdrem) ' sec remaining calculating RCV idx']);
    end
    close(wb)
% end
end

function transSenAll = calc_angular_sens(x_trans,y_trans,z_trans,x_img,y_img,z_img,elemWidthwl,Trans,nx,ny,nz)
x_img0 = repmat(x_img',[1 length(y_img) length(z_img)]);
y_img0 = repmat(y_img,[length(x_img) 1 length(z_img)]);
z_img0 = repmat(z_img',[1 length(x_img) length(y_img)]);
z_img0 = permute(z_img0,[2 3 1]);

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
        dotUV = nx(ielem)*(x_img1 - x_trans(ielem)) + ...
            ny(ielem)*(y_img1 - y_trans(ielem)) + ...
            nz(ielem)*(z_img1 - z_trans(ielem));


        normU = norm([nx(ielem) ny(ielem) nz(ielem)]);
        normV = sqrt((x_img1 - x_trans(ielem)).^2 + ...
            (y_img1 - y_trans(ielem)).^2 + ...
            (z_img1 - z_trans(ielem)).^2);
        temp = abs(dotUV./(normU*normV));
        temp(temp>1) = 1;
        %     temp(temp<-1) = -1;
        theta0 = acos(temp) + 0.0001;
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
function [us_bmode,us_bmode1] = GPUReconLoop(rfdataAll,idxAll,iframe,iline,us_bmode,transSenAll)
if iline == 1, us_bmode = 0*us_bmode; end
rfdata = rfdataAll(:,iline,iframe);
us_bmode1 = interp1(rfdata,idxAll(:,:,:,iline));
us_bmode = us_bmode + us_bmode1.*transSenAll(:,:,:,iline);
% us_bmode = us_bmode + us_bmode1;
end

function coefficients = generateSensitivityCoefficients(num_coefficients, mean_val, std_dev, seed)
    % Generate an array of coefficients to mimic transducer sensitivity fluctuations
    % 
    % Parameters:
    % num_coefficients - Number of coefficients to generate (int)
    % mean_val - Mean value of the sensitivity (float)
    % std_dev - Standard deviation for fluctuations (float)
    % seed - Random seed for reproducibility (int)
    %
    % Returns:
    % coefficients - Array of sensitivity coefficients

    if nargin < 4
        rng('shuffle'); % Use current time as the seed if not provided
    else
        rng(seed); % Set the seed for reproducibility
    end

    % Generate normally distributed sensitivity coefficients
    coefficients = mean_val + std_dev .* randn(1, num_coefficients);

    % Plot the coefficients to visualize the fluctuations
    figure;
    plot(coefficients, 'o-', 'LineWidth', 1.5, 'MarkerSize', 4);
    title('Fluctuation of Transducer Sensitivity Coefficients');
    xlabel('Coefficient Index');
    ylabel('Sensitivity Coefficient');
    grid on;
end