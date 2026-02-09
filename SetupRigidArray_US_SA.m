% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpL7_4FlashAngles.m - Example of plane wave imaing with
%                                       steering angle transmits
% Description:
%   Sequence programming file for L7-4 Linear array, using plane wave
%   transmits with multiple steering angles. All 128 transmit and receive
%   channels are active for each acquisition. Processing is asynchronous
%   with respect to acquisition.
%
% Last update:
% 12/13/2015 - modified for SW 3.0

clear all
close all
P.startDepth = 0;   % Acquisition depth in wavelengths
P.endDepth = 150;   % This should preferrably be a multiple of 128 samples.
% P.endDepth = 300;   % This should preferrably be a multiple of 128 samples.
mediumTemp = 21;
sos = round(1402.4 + 5.01*mediumTemp - 0.055*mediumTemp^2 + 0.00022*mediumTemp^3);
frame_num_rcv_yt = 10;
frame_num_img_yt = 2; 
transimpedance = 1000;
TEST_MODE = 1;

PW = 1;

% Compute storage matrix size
fc = 5.6;
lambda = sos/(fc*1e6)*1e3; % mm
AcqDepth = P.endDepth * lambda;
LateralWidth = 18.9*2;
DiagDist = sqrt(AcqDepth^2+LateralWidth^2);
DiagPnts = DiagDist/lambda;
% DiagPntsSample = round(8*DiagPnts*1.2/2); % used for Resource.RcvBuffer.rowsPerFrame
% % % % % % DiagPntsSample = 2560; % used for Resource.RcvBuffer.rowsPerFrame
DiagPntsSample = 128*8; % used for Resource.RcvBuffer.rowsPerFrame   need 3900 pts for 110mm


USVolt_yt = 20;
% TX_ELEMS = 1:256;
TX_ELEMS = [];
ind = 0;
for x_line = 1:3:16
    for y_line = 1:3:16
        ielem = (x_line-1)*16+y_line;
        ind = ind + 1;
        TX_ELEMS(ind) = ielem;
    end
end

if PW
    na = 1;
else
    na = length(TX_ELEMS);      % Set na = number of angles.
end

%% Define system parameters.
Resource.Parameters.numTransmit = 256;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 256;    % number of receive channels.
Resource.Parameters.speedOfSound = sos;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
Resource.Parameters.Connector = [1,2];
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

%% Specify Trans structure array.
para_path = 'C:\Davia\20251125_RT_feedback\';  % Matrix array  -  Vera PC
load([para_path,'x_trans.mat']);
load([para_path,'z_trans.mat']);
load([para_path,'y_trans.mat']);

Trans.id = -1;
Trans.name = '2D Array';
Trans.units = 'mm';
Trans.spacingMm = 2;  % pitch size
Trans.frequency = 5.6;
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

Trans.ElementPos(:,1) = Trans.ElementPos(:,1)+ abs(min(Trans.ElementPos(:,1)));
Trans.ElementPos(:,2) = -(Trans.ElementPos(:,2)+ abs(min(Trans.ElementPos(:,2))));
% % % % % % Trans.ElementPos(:,4) = azi_trans; %angle
% % % % % % Trans.ElementPos(:,5) = ele_trans;  %angle
Trans.maxHighVoltage = 40;
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
% Specify PData structure array.
PData(1).PDelta = [0.5, 0.5, 0.5];
PData(1).Size(3) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((Trans.spacing*(Trans.numelements-0)^0.5)/PData(1).PDelta(2));
PData(1).Size(1) = PData(1).Size(2);      % single image page
% PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)^0.5 /2,Trans.spacing*(Trans.numelements-1)^0.5 /2,P.startDepth]; % x,y,z of upper lft crnr.
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.
PData(1).Origin = [0,0,P.startDepth];

% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

%% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = DiagPntsSample*na*2; % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = frame_num_rcv_yt;    % 30 frames stored in RcvBuffer.
Resource.InterBuffer(1).numFrames = 1;   % one intermediate buffer needed.
Resource.ImageBuffer(1).numFrames = 1; %frame_num_img_yt;
z_plane = 100;
% XZ plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).Title = '3D Flash Image - XZ plane';
Resource.DisplayWindow(1).pdelta = 0.5;
Resource.DisplayWindow(1).Position = [0,480, ...
    ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta), ... % width
    ceil(PData(1).Size(3)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta)];    % height
Resource.DisplayWindow(1).Orientation = 'xz';
% Resource.DisplayWindow(1).ReferencePt = [0,0.0,0.0];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),PData(1).PDelta(2)*PData(1).Size(2)/2,PData(1).Origin(3)];
% Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0.0,-40];
Resource.DisplayWindow(1).Colormap = grayscaleCFImap;
Resource.DisplayWindow(1).splitPalette = 1;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).mode = '2d';
% 
% % YZ plane
Resource.DisplayWindow(2).Type = 'Verasonics';
Resource.DisplayWindow(2).Title = '3D Flash Image - YZ plane';
Resource.DisplayWindow(2).pdelta = Resource.DisplayWindow(1).pdelta;
Resource.DisplayWindow(2).Position = [660,480, ...
    ceil(PData(1).Size(1)*PData(1).PDelta(2)/Resource.DisplayWindow(2).pdelta), ... % width
    ceil(PData(1).Size(3)*PData(1).PDelta(3)/Resource.DisplayWindow(2).pdelta)];    % height
Resource.DisplayWindow(2).Orientation = 'yz';
Resource.DisplayWindow(2).ReferencePt = [PData(1).PDelta(1)*PData(1).Size(1)/2,-PData(1).Origin(2),PData(1).Origin(3)];
% Resource.DisplayWindow(2).ReferencePt = [0,-PData(1).Origin(2),-40];
Resource.DisplayWindow(2).Colormap = grayscaleCFImap;
Resource.DisplayWindow(2).splitPalette = 1;
Resource.DisplayWindow(2).AxesUnits = 'mm';
% Resource.DisplayWindow(2).mode = '2d';
%
% % XY plane
Resource.DisplayWindow(3).Type = 'Verasonics';
Resource.DisplayWindow(3).Title = '3D Flash Image - XY plane';
Resource.DisplayWindow(3).pdelta = Resource.DisplayWindow(1).pdelta;
Resource.DisplayWindow(3).Position = [0,40, ...
    ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(3).pdelta), ... % width
    ceil(PData(1).Size(1)*PData(1).PDelta(2)/Resource.DisplayWindow(3).pdelta)];    % height
Resource.DisplayWindow(3).Orientation = 'xy';
% Resource.DisplayWindow(3).ReferencePt = [PData(1).Origin(1),...  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     -PData(1).Origin(2),(P.endDepth-P.startDepth)/2];%PData.Region(end).Shape.oPAIntersect];
Resource.DisplayWindow(3).ReferencePt = [PData(1).Origin(1),...
    -PData(1).Origin(2),36];
%Resource.DisplayWindow(3).Colormap = grayscaleCFImap;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Resource.DisplayWindow(3).Colormap = gray(256);
Resource.DisplayWindow(3).splitPalette = 1;
Resource.DisplayWindow(3).AxesUnits = 'mm';
Resource.DisplayWindow(3).mode = '2d';
%% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,2,1];

if PW
    TX_APOD = ones(1,Trans.numelements);
else
    TX_APOD = zeros(length(TX_ELEMS),Trans.numelements);
    for itx = 1:na
        TX_APOD(itx,TX_ELEMS(itx)) = 1;
    end
end
% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', kaiser(Resource.Parameters.numTransmit,1)', ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Trans.numelements)), 1, na);
% - Set event specific TX attributes.
for n = 1:na   % na transmit events
    TX(n).Steer = [0.0,0.0];
    TX(n).Delay = computeTXDelays(TX(n));
    TX(n).Apod = TX_APOD(n,:);
end
%%
TPC(1).name = '2D';
TPC(1).maxHighVoltage = 40;
TPC(1).hv = USVolt_yt;
%% Specify TGC Waveform structure.
% TGC.CntrlPts = [0,141,275,404,510,603,702,782];
TGC.CntrlPts = [100,100,100,100,100,100,100,100];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays.
% - We need na Receives for every frame.
maxAcqLength = ceil(sqrt(P.endDepth^2 + (15*Trans.spacing)^2));
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength,...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, na*Resource.RcvBuffer(1).numFrames);

% - Set event specific Receive attributes for each frame.
for i = 1:Resource.RcvBuffer(1).numFrames
    Receive(na*(i-1)+1).callMediaFunc = 1;
    for j = 1:na
        Receive(na*(i-1)+j).framenum = i;
        Receive(na*(i-1)+j).acqNum = j;
    end
end

% Specify Recon structure arrays.
% - We need one Recon structures which will be used for each frame.
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums', 1:na);

% Define ReconInfo structures.
% We need na ReconInfo structures for na steering angles.
ReconInfo = repmat(struct('mode', 'accumIQ', ...  % default is to accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1), 1, na);
% - Set specific ReconInfo attributes.
if na>1
    ReconInfo(1).mode = 'replaceIQ'; % replace IQ data
    for j = 1:na  % For each row in the column
        ReconInfo(j).txnum = j;
        ReconInfo(j).rcvnum = j;
    end
    ReconInfo(na).mode = 'accumIQ_replaceIntensity'; % accum and detect
else
    ReconInfo(1).mode = 'replaceIntensity';
end

% Specify Process structure array.
pers = 0;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',1.0,...            % pgain is image processing gain
                         'reject',2,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',60,...
                         'mappingMethod','full',...
                         'display',0,...      % display image after processing
                         'displayWindow',1};
Process(2).classname = 'Image';
Process(2).method = 'imageDisplay';
Process(2).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',1.0,...            % pgain is image processing gain
                         'reject',2,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','full',...
                         'display',0,...      % display image after processing
                         'displayWindow',2};
Process(3).classname = 'Image';
Process(3).method = 'imageDisplay';
Process(3).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',100,...            % pgain is image processing gain
                         'reject',0,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','log',... 
                         'compressFactor',20,...
                         'mappingMethod','full',...
                         'display',0,...      % display image after processing
                         'displayWindow',3};
% Process(4): 自定义 XY 平面显示
Process(4).classname = 'External';
Process(4).method    = 'myImageDisplayFunction'; % 您的自定义函数名
Process(4).Parameters = {'srcbuffer', 'image', ... % 源数据是 ImageBuffer
                         'srcbufnum', 1, ...
                         'srcframenum', -1, ...    % 处理最近的帧
                         'dstbuffer', 'none'};     % 没有输出目的地
% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
% SeqControl(2).argument = 900;  % 160 usec
SeqControl(2).argument = 150;  % 160 usec
SeqControl(3).command = 'timeToNextAcq';  % time between frames
% SeqControl(3).argument = 40000 - (na-1)*200;  % 20 msec
% SeqControl(3).argument = 20000;  % 20 msec
SeqControl(3).argument = 150;  % 20 msec
SeqControl(4).command = 'returnToMatlab';
SeqControl(5).command = 'triggerIn';   % the sequencer will pause and wait for the trigger input signal immediately before
SeqControl(5).condition = 'Trigger_2_Rising';
SeqControl(6).command = 'sync';
SeqControl(6).argument = 30000000;
SeqControl(7).command = 'stop';
nsc = 8; % nsc is count of SeqControl objects

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:na                     
        Event(n).info = 'Full aperture.';
        Event(n).tx = j;
        Event(n).rcv = na*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = [2];
        n = n+1;
    end

    if ~TEST_MODE
        Event(n).info = 'sync';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 6;
        n = n+1;
    end

    Event(n).info = 'Transfer To Host';
    Event(n).tx = 0;         % use 1st TX structure.
    Event(n).rcv = 0;    % use 1st Rcv structure of frame.
    Event(n).recon = 0;      % no reconstruction.
    Event(n).process = 0;    % no processing
    Event(n).seqControl = nsc;
    SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
    nsc = nsc+1;
    n = n + 1;
    
    if TEST_MODE
%         Event(n).info = 'recon and process';
%         Event(n).tx = 0;
%         Event(n).rcv = 0;
%         Event(n).recon = 1;
%         Event(n).process = 1;
%         Event(n).seqControl = 0;
%         n = n + 1;
        
%         Event(n).info = 'recon and process';
%         Event(n).tx = 0;
%         Event(n).rcv = 0;
%         Event(n).recon = 1;
%         Event(n).process = 2;
%         Event(n).seqControl = 0;
%         n = n + 1;
% 
        Event(n).info = 'recon and process';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 1;
        Event(n).process = 3;
        Event(n).seqControl = 0;
        n = n + 1;

        % 替换原有 Event(n).process = 3 的事件，改为 Process(4)
Event(n).info    = 'Perform custom XY image display.';
Event(n).tx      = 0;
Event(n).rcv     = 0;
Event(n).recon   = 0;
Event(n).process = 4; % <-- 引用新的 Process(4)
Event(n).seqControl = 0;
        n = n + 1;
    end

%     if floor(i/5) == i/5 
        Event(n).info = 'return to matlab';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 4;
        n = n + 1;
%     end
end

if TEST_MODE
    Event(n).info = 'Jump back';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 1;
else
    Event(n).info = 'stop';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 7;
end

% User specified UI Control Elements
% - Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%SensCutoffCallback');

% - Range Change
MinMaxVal = [64,300,P.endDepth]; % default unit is wavelength
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        AxesUnit = 'mm';
        MinMaxVal = MinMaxVal * (Resource.Parameters.speedOfSound/1000/Trans.frequency);
    end
end
UI(2).Control = {'UserA1','Style','VsSlider','Label',['Range (',AxesUnit,')'],...
                 'SliderMinMaxVal',MinMaxVal,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
UI(2).Callback = text2cell('%RangeChangeCallback');

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 5;

%% Save all the structures to a .mat file.
save('MatFiles/Soft2DArray');
% VSX

return
%%
% **** Callback routines to be converted by text2cell function. ****
%SensCutoffCallback - Sensitivity cutoff change
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
%SensCutoffCallback

%RangeChangeCallback - Range change
simMode = evalin('base','Resource.Parameters.simulateMode');
% No range change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',evalin('base','P.endDepth'));
    return
end
Trans = evalin('base','Trans');
Resource = evalin('base','Resource');
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

P = evalin('base','P');
P.endDepth = UIValue;
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        P.endDepth = UIValue*scaleToWvl;
    end
end
assignin('base','P',P);

evalin('base','PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));');
evalin('base','PData(1).Region = computeRegions(PData(1));');
evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');
Receive = evalin('base', 'Receive');
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
for i = 1:size(Receive,2)
    Receive(i).endDepth = maxAcqLength;
end
assignin('base','Receive',Receive);
evalin('base','TGC.rangeMax = P.endDepth;');
evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','Receive','TGC','Recon'};
assignin('base','Control', Control);
assignin('base', 'action', 'displayChange');
return
%RangeChangeCallback





