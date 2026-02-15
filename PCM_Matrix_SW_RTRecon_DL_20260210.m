clear all
g = gpuDevice(1);
reset(g);
addpath(genpath('C:\Davia\Davia_codeOnly\'))
%% Use C drive or SSD to store data

%% === Set commonly modified parameters ========================================
% Saving options
targetdir = 'C:\Davia\20260210\';     
targetfolder = 't1';  

TEST_MODE = 1; %
SIM_MODE = 0;
P(1).startDepth =0;          % Acquisition start depth in wavelengths
P(1).endDepth = 2000;        % Acquisition end depth    
P(1).startAcqDepth = 0;      % 300-500 us for first peak, 600-900 us for second peak
P(1).numFrames = 10;

TXflag = 0; % tx=1
mediumTemp = 21; % temperature in Celsius
decimSampleNum = 1;
BPfilter = 0;
% % % caculation
% Trans.frequency=4;
% T=21;
% sos = 1402.4+5.01*T-0.055*T^2+0.00022*T^3;
% waveLength = (sos/1000)/Trans.frequency;
% acqtime = 2*1300*waveLength*1e-3/sos*15.625*1e6;

% Refresh coordinates.mat to [0,0,0,0]
a = [0,0,0,0];
% CoordinatesSavedir = 'C:\Davia\20251207_Bmode_DB\';   %'C:\Users\PI-Lab\Desktop\PyDobot\';
% CoordinatesSavename = [CoordinatesSavedir,'coordinates.mat'];
% save(CoordinatesSavename,'a')
%% Continuous saving directory and corresponding flags
global fdir;
fdir = [targetdir,targetfolder,'\'];

if TEST_MODE
    RT_RECON = 1;     % in test mode, real-time recon is enable
    P(1).numFrames = 1;    % in test mode, number of frames are small for quicker VSX preprocessing
    P(2).numFrames = 1;
    RT_SAVING = 0;
else
    RT_RECON = 0;     % real-time recon or not
    RT_SAVING = 1;
end

if RT_SAVING
    if ~exist(fdir,'dir')
        mkdir(fdir)
    else
        error('FOLDER ALREADY EXIST. DATA WILL BE OVERWRITTEN.')
    end
end

%% delay time setting
sos = 1402.4+5.01*mediumTemp-0.055*mediumTemp^2+0.00022*mediumTemp^3;
flash2Qdelay=0;  %laser control for quantel laser

na = 1; % Set na = number of flash angles for 2D.
ne = 0; % ne = number of acquisitions in PA ensemble for coherent addition in I/Q buffer.

dtheta2D = 0;
oneway = 1;     % (logical) oneway=1 turns off the transmitters by setting TX.Apod to zero for all transmitters

PA_Angle = 0;   % angle of transmit plane wave for use in testing the PA mode in simulation
PA_PRF = 1;   % PA PRF in Hz. To be set in accordance with laser rep rate, when not using the input trigger mode.
if oneway==1
    disp(' *** PCM mode: Using one-way F-only reconstruction ***')
else
    disp(' *** Ultrasound Transmit mode: Using conventional T/R reconstruction ***')
end

%% Specify system parameters(Vantage 128)
Resource.Parameters.numTransmit = 256;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 256;    % number of receive channels.
Resource.Parameters.speedOfSound = sos;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = SIM_MODE;
Resource.Parameters.Connector = [1,2];  % use both of the connectors

%% Specify Trans structure array--Customized 2D array
% para_path='C:\Users\daiwe\Dropbox\CavitationMapping_2024_Davia_Chaorui\PCM_CODE_0619-0625\US_Para\';
% % % % % % para_path = 'E:\Davia_codeOnly\US_Para\';   % Vera PC
para_path = 'D:\Davia\20250716_new5MHz\';  % Matrix array  -  Vera PC
load([para_path,'x_trans.mat']);
load([para_path,'z_trans.mat']);
load([para_path,'y_trans.mat']);
% % % % % % % load([para_path,'ele_trans.mat']);
% % % % % % % load([para_path,'azi_trans.mat']);

Trans.id = -1;
Trans.name = '2D Array';
Trans.units = 'mm';
Trans.spacingMm = 1;
Trans.frequency = 5.8;
Trans.wavelength = (sos/1000)/Trans.frequency;

Trans.numelements = 256;
Trans.elementWidth = 1;
Trans.type = 2;
Trans.impedance = 1e3;
Trans.connType = 1;
Trans.ElementPos = zeros(Trans.numelements,5);
Trans.ElementPos(:,1) = x_trans;
Trans.ElementPos(:,2) = y_trans;
Trans.ElementPos(:,3) = z_trans;
% % % % % % Trans.ElementPos(:,4) = azi_trans; %angle
% % % % % % Trans.ElementPos(:,5) = ele_trans;  %angle
Trans.maxHighVoltage = 20;
Trans.radiusMm = 40;
Trans.lensCorrection = 0;

%Intermediate Variables
if strcmp(Trans.units,'mm')
    Trans.ElementPosMm = Trans.ElementPos;
    Trans.ElementPosWL = Trans.ElementPos./Trans.wavelength;
    Trans.spacing = Trans.spacingMm./Trans.wavelength;
    Trans.radius = Trans.radiusMm./Trans.wavelength;
else
    Trans.ElementPosMm = Trans.ElementPos.*Trans.wavelength;
    Trans.ElementPosWL = Trans.ElementPos;
end

eleWidthWl = Trans.elementWidth ./ Trans.wavelength;
Theta = (-pi/2:pi/100:pi/2);
Theta(51) = 0.0000001;
Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./...
    (eleWidthWl*pi*sin(Theta))));  % set a reasonable high voltage limit.
%% Setup real time recon matrix
global sensarr area1 delay_t GPU temperature_C shiftVal map_buffer
% sensarr: position of the elements 
% areal: imaging area
map_buffer = ones(59,53);  %y_range,z_range
GPU = 1;
ResFactor = 3.25; %6.25
pixel_size = 1/ResFactor;
x_range = [-9 9];
y_range = [-9 9];% [-5 5]
z_range = [42 58];  %[26 34]
x0 = x_trans; xintval = 1/ResFactor;
y0 = y_trans; yintval = 1/ResFactor;
z0 = z_trans; zintval = 1/ResFactor;

delay_t = 0;
sensarr = SensorArray2D_PCM(y0,x0,z0);
sensarr.x0 = sensarr.x0';
sensarr.y0 = sensarr.y0';
sensarr.z0 = sensarr.z0';
area1 = ImageArea2D_PCM(y_range(2), y_range(1), yintval, ...
    x_range(2), x_range(1), xintval, ...
    z_range(2), z_range(1), zintval);
temperature_C = mediumTemp;
shiftVal = 150*Trans.frequency*4;  % 2400 ?

%% Specify PData structure arrays.
% - 2D PData structure
PData(1).PDelta(1) = 1.0;
PData(1).PDelta(3) = 0.5;
PData(1).Size(1,1) = ceil((P(1).endDepth-P(1).startDepth)/PData(1).PDelta(3)); % rows
PData(1).Size(1,2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1)); % cols
PData(1).Size(1,3) = 1;             % single image page
PData(1).Origin = [-Trans.spacing*63.5,0,P(1).startDepth]; % x,y,z of uppr lft crnr.

%% Specify Media object and point displacement function
pt1;
Media.function = 'movePoints';

%% Specify Resources.
% - RcvBuffer(1) is for both 2D and PA acquisitions.
fs = 4*Trans.frequency;
ds = (1/(fs*10^6))*sos;
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = ceil(2*(P(1).endDepth-P(1).startDepth)*Trans.wavelength*10^-3 / ds / 128)*128 *(na + ne)*P(1).numFrames;    %16000*2
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 1;           % 20 frames allocated for RF acqusitions.
% InterBuffer(1) is for 2D reconstructions.
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;          % one intermediate frame needed for 2D.
% InterBuffer(2) is for PA reconstructions.
Resource.InterBuffer(2).datatype = 'complex';
Resource.InterBuffer(2).numFrames = 1;          % one intermediate frame needed for PA.
% ImageBuffer(1) is for 2D image.
Resource.ImageBuffer(1).datatype = 'double';    % image buffer for 2D
Resource.ImageBuffer(1).numFrames = P.numFrames;
% ImageBuffer(2) is for PA image.
Resource.ImageBuffer(2).datatype = 'double';    % image buffer for PA
Resource.ImageBuffer(2).numFrames = P.numFrames;
% DisplayWindow is for 2D combined with PA
% Resource.DisplayWindow(1).Title = mfilename;
% Resource.DisplayWindow(1).pdelta = 0.4;
% ScrnSize = get(0,'ScreenSize');
% DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
% DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
% Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
%     DwWidth, DwHeight];
% Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)]; % 2D imaging is in the X,Z plane
% Resource.DisplayWindow(1).Type = 'Verasonics';
% Resource.DisplayWindow(1).numFrames = 20;
% Resource.DisplayWindow(1).AxesUnits = 'mm';
% Resource.DisplayWindow(1).Colormap = grayscaleCFImap;
% Resource.DisplayWindow(1).splitPalette = 1;
% Resource.DisplayWindow(1).Orientation = 'xy';

%% Specify Transmit waveforms structure
% - 2D transmit waveform
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,0.67,2,1];
% - PA transmit waveform
TW(2).type = 'parametric';
TW(2).Parameters = [Trans.frequency,0.67,6,1];

%% Specify Transmit beams structure
TX = repmat(struct('waveform', 1, ...
    'Origin', [0.0,0.0,0.0], ...
    'focus', 0.0, ...
    'Steer', [0.0,0.0], ...
    'Apod', zeros(1,Trans.numelements), ...
    'Delay', zeros(1,Trans.numelements)), 1, na+ne); % na TXs for 2D + 1 for PA
% - Set event specific TX attributes.
if fix(na/2) == na/2       % if na even
    startAngle = (-(fix(na/2) - 1) - 0.5)*dtheta2D;
else
    startAngle = -fix(na/2)*dtheta2D;
end
for n = 1:na   % na transmit events for 2D
    TX(n).Steer = [(startAngle+(n-1)*dtheta2D),0.0];
    TX(n).Delay = computeTXDelays(TX(n));
end
% -- only one TX struct needed for PA
TX(na+1).waveform = 2;
if oneway
    TX(na+1).Apod = zeros(1,Trans.numelements);     % THIS COMMAND TURNS OFF ALL TRANSMITTERS AND INVOKES THE RECEIVE-ONLY BEAMFORMER
else
    TX(na+1).Apod =  ones(1,Trans.numelements);     % This is the conventional T/R condition and invokes the default beamformer
end
TX(na+1).Steer = [PA_Angle,0.0];            % only relevant when transmitters are active
TX(na+1).Delay = computeTXDelays(TX(na+1)); % only relevant when transmitters are active
if TXflag == 0
    TX(1).Apod = zeros(1,Trans.numelements);
end
%% Specify TPC structures.
TPC(1).name = '2D';
TPC(1).maxHighVoltage = 30;
TPC(1).hv = 3;
% This allows one to use different transmit profile for PA ... only relevant if transmitters are active
% --- currently TPC(2) is not used ---
TPC(2).name = 'PA';
TPC(2).maxHighVoltage = 30;

%% Analog front end gain settings.
RcvProfile(1).LnaGain = 18;     % 12, 18, or 24 dB  (18=default)
RcvProfile(1).condition = 'immediate';

RcvProfile(2).LnaGain = 18;
RcvProfile(2).condition = 'immediate';

if BPfilter == 1
    BPF = [];
else
    BPF1 = zeros(1,20);
    BPF=[BPF1,1];
end
%% Specify Receive structure arrays.
%   We need to acquire all the 2D and PA data within a single RcvBuffer frame.  This allows
%   the transfer-to-host DMA after each frame to transfer a large amount of data, improving throughput.
% - We need 2*na Receives for a 2D frame and ne Receives for a PA frame.
maxAcqLngth2D = sqrt(P(1).endDepth^2 + (sqrt(Trans.numelements)*Trans.spacing)^2) - P(1).startDepth;
wl4sPer128 = 128/(4*2);  % wavelengths in a 128 sample block for 4 smpls per wave round trip.
% wl2sPer128 = 128/(2*2);  % wavelengths in a 128 sample block for 2 smpls per wave round trip.
if TEST_MODE
    Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
    'startDepth', P(1).startDepth, ...'endDepth', P(1).startDepth + wl4sPer128*ceil(maxAcqLngth2D/wl4sPer128), ...
    'endDepth', P(1).endDepth,...
    'TGC', 1, ...           % TGC(1) is tied to the GUI sliders
    'bufnum', 1, ...
    'framenum', 1, ...
    'acqNum', 1, ...
    'sampleMode','NS200BW', ...%
    'InputFilter', BPF, ...
    'mode', 0, ...
    'callMediaFunc', 0), 1, (na+ne)*P(1).numFrames);
else 
    Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
    'startDepth', P(1).startAcqDepth, ...'endDepth', P(1).startDepth + wl4sPer128*ceil(maxAcqLngth2D/wl4sPer128), ...
    'endDepth', P(1).endDepth,...
    'TGC', 1, ...           % TGC(1) is tied to the GUI sliders
    'bufnum', 1, ...
    'framenum', 1, ...
    'acqNum', 1, ...
    'sampleMode','NS200BW', ...%
    'InputFilter', BPF, ...
    'mode', 0, ...
    'callMediaFunc', 0), 1, (na+ne)*P(1).numFrames);
end
% - Set event specific Receive attributes.
for i = 1:P(1).numFrames
    k = (na + ne)*(i-1); % k keeps track of Receive index increment per frame.
    % - Set attributes for each frame.
    Receive(k+1).callMediaFunc = 1; % move points before doing ensemble of different angle plane waves
    % acquisitions for 2D
    for j = 1:na
        Receive(j+k).framenum = 1;
        Receive(j+k).acqNum = i;
    end
    % PA acquisitions
%     for j = (na+1):(na+ne)
%         Receive(j+k).framenum = i;
%         Receive(j+k).acqNum = j;        % PA acqNums continue after 2D
%         Receive(j+k).startDepth = P(2).startDepth;
%         Receive(j+k).endDepth = P(2).startDepth + wl4sPer128*ceil(maxAcqLngthPA/wl4sPer128);
%         Receive(j+k).TGC = 2;           % TGC(1) is tied to the GUI sliders
%         if j==na+1, Receive(j+k).callMediaFunc = 1; end % move points between 2D and PA to see difference in simulation
%     end
end

%% Specify TGC Waveform structures.
% - 2D TGC
% TGC.CntrlPts = [0,141,275,404,510,603,702,782];
TGC.CntrlPts = [100 100 100 100 100 100 100 100];
TGC(1).rangeMax = P(1).endDepth;
TGC(1).Waveform = computeTGCWaveform(TGC(1));

%% Specify Recon structure arrays.
% - We need two Recon structures, one for 2D, one for PA. These will be referenced in the same
%   event, so that they will use the same (most recent) acquisition frame.
Recon = repmat(struct('senscutoff', 0.7, ...
    'pdatanum', 1, ...
    'rcvBufFrame', -1, ...
    'IntBufDest', [1,1], ...
    'ImgBufDest', [1,-1], ...
    'RINums', zeros(1,1)), 1, 1);
% - Set Recon values for 2D frame.
Recon(1).RINums(1,1:na) = 1:na;  % na ReconInfos needed for na angles
k = na + 1;
% - Set Recon values for PA ensemble.
% Recon(2).pdatanum = 2;
% Recon(2).IntBufDest = [2,1];
% Recon(2).ImgBufDest = [2,-1];
% Recon(2).RINums(1,1:ne) = k:(k+ne-1);   % 'ne' ReconInfos needed for PA ensemble.

%% Define ReconInfo structures.
% - For 2D, we need na ReconInfo structures for na steering angles.
% - For PA, we need ne ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...    % 4=accumulate IQ data.
    'txnum', 1, ...
    'rcvnum', 1, ...
    'regionnum', 1), 1, na + ne);
% - Set specific ReconInfo attributes.
%   - ReconInfos for 2D frame.
ReconInfo(1).mode = 'replaceIQ';          % 3=replace IQ data (expect to use mode 5 on last acquisition)
for j = 1:na
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = j;
end
if na>1
    ReconInfo(na).mode = 'accumIQ_replaceIntensity';     % 5=Reconstruct IQ data, add values to InterBuffer and compute magnitude, replacing data in ImageBuffer.
else
    ReconInfo(na).mode = 'replaceIntensity';     % 0=replace IQ data, detect, and replace Intensity data in ImageBuffer. (single acquisition)
end

% %  - ReconInfos for PA ensemble.
% k = na;
% for j = 1:ne
%     if j==1, ReconInfo(k+j).mode = 'replaceIQ'; end
%     ReconInfo(k+j).txnum = na + 1;
%     ReconInfo(k+j).rcvnum = na + j;
% end
% if ne>1
%     ReconInfo(na+ne).mode = 'accumIQ_replaceIntensity'; % 5=accum and detect
% else
%     ReconInfo(na+ne).mode = 'replaceIntensity'; % 0=replace IQ data, detect, and replace Intensity data;  1=Add the new reconstructed intensity data to the data in the ImageBuffer
% end

%% Specify Process structure arrays.
cpt = 22;       % define here so we can use in UIControl below
cpers = 80;     % define here so we can use in UIControl below

Process(1).classname = 'External';
Process(1).method = 'saveRcvDataRT_feedback_2D_PCM'; % 'saveRcvDataRealTime2D_PCM';
Process(1).Parameters = {'srcbuffer','receive',...
    'srcbufnum',1,...
    'srcframenum',-1,...
    'dstbuffer','none'};

Process(2).classname = 'External';
Process(2).method = 'RTReconMatMult2D_PCM_SW_v4_fastplot_matrix';
Process(2).Parameters = {'srcbuffer','receive',...
    'srcbufnum',1,...
    'srcframenum',-1,...
    'dstbuffer','none'};

%% Specify SeqControl structure arrays.
% -- Time between 2D flash angle acquisitions
SeqControl(1).command = 'timeToNextAcq';
SeqControl(1).argument = 400;  %UNIT is us, from the start of this DAQ to the start of next one
% -- Change to Profile 2 (PA)
SeqControl(2).command = 'setTPCProfile';
SeqControl(2).condition = 'next';
SeqControl(2).argument = 2;
% -- Time between 2D acquisition and PA ensemble. Set to allow time for profile change.
SeqControl(3).command = 'timeToNextAcq';
SeqControl(3).argument = 10; % time in usec  %Time between US and PA
% ACM_DAQ= maxAcqLngth2D*8/20  %[us]
% Total_ACM=SeqControl(3).argument;
% -- PRF for PA ensemble
SeqControl(4).command = 'timeToNextAcq';
SeqControl(4).argument = round(1/(PA_PRF*1e-06)); % us  (10 msecs for PA_PRF=100 Hz)
% -- Change to Profile 1 (2D)
SeqControl(5).command = 'setTPCProfile';
SeqControl(5).condition = 'next';
SeqControl(5).argument = 1;
% -- Time between PA and next 2D acquisition. Set to allow time for profile change.
SeqControl(6).command = 'timeToNextAcq';
SeqControl(6).argument = 7000; % time in usec
% -- Jump back to start.
SeqControl(7).command = 'jump';
SeqControl(7).argument = 1;
% set receive profile
SeqControl(8).command = 'setRcvProfile';
SeqControl(8).argument = 1;
SeqControl(9).command = 'setRcvProfile';
SeqControl(9).argument = 2;
% output trigger
SeqControl(10).command = 'triggerOut';
% SeqControl(10).argument =0; % 250*13.8;  %250 is 1us
% input trigger
SeqControl(11).command = 'triggerIn';
SeqControl(11).condition = 'Trigger_1_Rising'; % Trigger input 1, enable with rising edge
SeqControl(11).argument = 0; % 500 msec timeout delay
% (Timeout range is 1:255 in 250 msec steps; 0 means timeout disabled)
% noop delay between trigger in and start of acquisition
SeqControl(12).command = 'noop';
SeqControl(12).argument = fix(flash2Qdelay)*5; % noop counts are in 0.2 microsec increments
% sync command
SeqControl(13).command = 'sync';
SeqControl(13).argument = 20000000; % 30 sec timeout for software sequencer (default is 0.5 seconds)
% - The remainder of the SeqControl structures are defined dynamically in the sequence events.
%   The variable nsc keeps track of the next SeqControl structure index.

nsc = 14;  % next SeqControl number
% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 3;

%% Specify Event structure arrays.
n = 1;
Event(n).info = 'noop'; %wait the bubble to grow
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 13;
n = n+1;

for i = 1:Resource.RcvBuffer(1).numFrames
    % Wait for input trigger from flash lamp firing
%     Event(n).info = 'Wait for Trigger IN';
%     Event(n).tx = 0;
%     Event(n).rcv = 0;
%     Event(n).recon = 0;
%     Event(n).process = 0;
%     Event(n).seqControl = 11;
%     n = n+1;
    
    %set delay time before sending out interrogation pulse
%     Event(n).info = 'noop'; %wait the bubble to grow
%     Event(n).tx = 0;
%     Event(n).rcv = 0;
%     Event(n).recon = 0;
%     Event(n).process = 0;
%     Event(n).seqControl = 13;
%     n = n+1;
    
    % Acquire PCM frame
    for j = 1:na           %1:5   na=1
        Event(n).info = 'Acquire PCM';
        Event(n).tx = j;
        Event(n).rcv = (na+ne)*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = [11];   %
        n = n+1;
    end
end
%     if na~=1
%         Event(n-na).seqControl = 1;
%         Event(n-1).seqControl =[9,3,10];   %(1/3)% replace last 2D acquisition Event's seqControl (longer TTNA and new RCV profile)
%     else
%         Event(n-1).seqControl=[9,3,10];  %(1/3) [9,3]
%     end

    if RT_RECON            
        Event(n).info = 'Transfer Data';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = nsc;
        n = n+1;
        SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer (each needs a different value of nsc)
        nsc = nsc+1;
        
        Event(n).info = ['3D Reconstruct Frame ' num2str(i)];
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 2;
        Event(n).seqControl = 0;
        n = n+1;
    else
%         Event(n).info = 'Transfer To Host';
%         Event(n).tx = 0;         % use 1st TX structure.
%         Event(n).rcv = 0;    % use 1st Rcv structure of frame.
%         Event(n).recon = 0;      % no reconstruction.
%         Event(n).process = 0;    % no processing
%         Event(n).seqControl = nsc; % set wait time and transfer data
%         SeqControl(nsc).command = 'transferToHost';
%         nsc = nsc + 1;
%         n = n+1;     
        Event(n).info = 'Transfer To Host';
        Event(n).tx = 0;         % use 1st TX structure.
        Event(n).rcv = 0;    % use 1st Rcv structure of frame.
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = [nsc]; % set wait time and transfer data
        SeqControl(nsc).command = 'transferToHost';
        tf_nsc = nsc;
        nsc = nsc + 1;
        SeqControl(nsc).command = 'waitForTransferComplete';
        SeqControl(nsc).argument = tf_nsc;
        nsc = nsc + 1;
        SeqControl(nsc).command = 'markTransferProcessed';
        SeqControl(nsc).argument = tf_nsc;
        mtf = nsc;
        nsc = nsc+1;
        n = n+1;
        
        Event(n).info = 'sync'; %wait the bubble to grow
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 13;
        n = n+1;
    end
    
%     Event(n).info = 'Return to Matlab';
%     Event(n).tx = 0;
%     Event(n).rcv = 0;
%     Event(n).recon = 0;
%     Event(n).process = 0;
%     Event(n).seqControl = nsc;
%     SeqControl(nsc).command = 'returnToMatlab';
%     nsc = nsc + 1;
%     n = n+1;
    
    if RT_SAVING
        Event(n).info = 'Realtime saving of PA data';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 1;
        Event(n).seqControl = 0;
        n = n+1;
    end
    
%     Event(n).info = 'HW/SW Sync';
%     Event(n).tx = 0;         % use 1st TX structure.
%     Event(n).rcv = 0;      % use ith Rcv structure.
%     Event(n).recon = 0;      % no reconstruction.
%     Event(n).process = 0;    % no processing
%     Event(n).seqControl = nsc;
%     SeqControl(nsc).command = 'sync';
%     nsc = nsc + 1;
%     n = n+1;

% %
Event(n).info = 'Jump back';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 7;


%% User specified UI Control Elements
% - Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
    'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
    'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%-UI#1Callback');

% - Color Priority Threshold Slider
UI(2).Control = {'UserB2','Style','VsSlider','Label','Color Priority','SliderMinMaxVal',[0,255,cpt],...
    'SliderStep',[1/255,0.1],'ValueFormat','%3.0f'};
UI(2).Callback = text2cell('%-UI#2Callback');

% - Color Persistence Slider
UI(3).Control = {'UserB1','Style','VsSlider','Label','Color Persistence','SliderMinMaxVal',[0,100,cpers],...
    'SliderStep',[1/100,0.1],'ValueFormat','%3.0f'};
UI(3).Callback = text2cell('%-UI#3Callback');

% - Wavelength delay SOS 
UI(4).Control =  {'UserB4','Style','VsSlider','Label','RF Delay (ind)',...
                  'SliderMinMaxVal',[shiftVal-1500,shiftVal+1500,shiftVal],...
                  'SliderStep',[0.01,0.1],'ValueFormat','%.2f'};
UI(4).Callback = text2cell('%-UI#4Callback');
%% Save all the structures to a .mat file, and run VSX automatically
filename = ('Array2DPCM_ACM');   % define variable 'filename' to permit VSX to skip user query for matfile
save (['MatFiles/',filename])                 % save the structures to a matfile
% VSX                             % invoke VSX automatically when running this Setup script

return


%% **** Callback routines to be encoded by text2cell function. ****
% ---------------------------------------------------------------------
%-UI#1Callback - Sensitivity cutoff change
ReconL = evalin('base', 'Recon');
for i = 1:size(ReconL,2)
    ReconL(i).senscutoff = UIValue;
end
assignin('base','Recon',ReconL);
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'Recon'};
assignin('base','Control', Control);
return
%-UI#1Callback


%-UI#2Callback - Color Threshold change
% Set the value in the Process structure for use in cineloop playback.
Process = evalin('base','Process');
for k = 1:2:length(Process(2).Parameters)
    if strcmp(Process(2).Parameters{k},'threshold'), Process(2).Parameters{k+1} = UIValue; end
end
assignin('base','Process',Process);
% Set Control.Command to set Image.threshold.
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Process',2,'threshold',UIValue};
assignin('base','Control', Control);
%-UI#2Callback

%-UI#3Callback - Color Persistence change
% Set the value in the Process structure for use in cineloop playback.
Process = evalin('base','Process');
for k = 1:2:length(Process(2).Parameters)
    if strcmp(Process(2).Parameters{k},'persistLevel'), Process(2).Parameters{k+1} = UIValue; end
end
assignin('base','Process',Process);
% Set Control.Command to set Image.persistLevel.
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Process',2,'persistLevel',UIValue};
assignin('base','Control', Control);
%-UI#3Callback

%-UI#4Callback - Wavelength delay for RT recon change
assignin('base','delay_t',UIValue);
%-UI#4Callback