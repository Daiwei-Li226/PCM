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

clear all
P.startDepth = 0;   % Acquisition depth in wavelengths
% P.endDepth = 128;   % This should preferrably be a multiple of 128 samples.
P.endDepth = 300;   % This should preferrably be a multiple of 128 samples.
mediumTemp = 21;
sos = round(1402.4 + 5.01*mediumTemp - 0.055*mediumTemp^2 + 0.00022*mediumTemp^3);
frame_num_rcv_yt = 10;
frame_num_img_yt = 10; 
transimpedance = 2000;  %1000

% Compute storage matrix size
fc = 4.5;
lambda = sos/(fc*1e6)*1e3; % mm
AcqDepth = P.endDepth * lambda;
LateralWidth = 18.9*2;
DiagDist = sqrt(AcqDepth^2+LateralWidth^2);
DiagPnts = DiagDist/lambda;
% DiagPntsSample = round(8*DiagPnts*1.2/2); % used for Resource.RcvBuffer.rowsPerFrame
DiagPntsSample = 4096; % used for Resource.RcvBuffer.rowsPerFrame
USVolt_yt = 20;
% TX_ELEMS = 1:256;
% TX_ELEMS = [];
ind = 0;
% for x_line = 1:3:16
%     for y_line = 1:3:16
for x_line = 1:1:16
    for y_line = 1:1:16
        ielem = (x_line-1)*16+y_line;
        ind = ind + 1;
        TX_ELEMS(ind) = ielem;
%         text(x_trans(ielem),y_trans(ielem),num2str(ielem));
    end
end
% na = 36;      % Set na = number of angles.
na = 256;      % Set na = number of angles.
% na = 64;      % Set na = number of angles.
if (na > 1), dtheta = (36*pi/180)/(na-1); P.startAngle = -36*pi/180/2; else dtheta = 0; P.startAngle=0; end % set dtheta to range over +/- 18 degrees.

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
Trans.name = 'Deformable2D_6x6';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
% Trans = computeTrans(Trans);  % L7-4 transducer is 'known' transducer so we can use computeTrans.
% Trans.maxHighVoltage = 50;  % set maximum high voltage limit for pulser supply.
Trans.id = -1;
Trans.units = 'mm';
Trans.lensCorrection = 0;
Trans.frequency = fc;
Trans.spacingMm =1;
Trans.elementWidth = 6.73;
Trans.numelements = 256;
Trans.connType = 1;
Trans.type = 2;
Trans.maxHighVoltage = 20;
Trans.radiusMm = 0;
Trans.impedance = transimpedance;
% x0 = linspace(-4.5,6,8);
% y0 = linspace(-10.5,12,16); 
% [x_trans,y_trans] = meshgrid(x0,y0);
% x_trans = reshape(x_trans,1,[]);
% y_trans = reshape(y_trans,1,[]);
% z_trans = zeros(size(x_trans));
para_path = 'C:\Davia\20251125_RT_feedback\';  % Matrix array  -  Vera PC
load([para_path,'x_trans.mat']);
load([para_path,'z_trans.mat']);
load([para_path,'y_trans.mat']);

Trans.ElementPos = zeros(Trans.numelements,5);
Trans.ElementPos(:,1) = x_trans;
Trans.ElementPos(:,2) = y_trans;
Trans.ElementPos(:,3) = z_trans;
Trans.ElementPos(:,4) = 0; %angle
Trans.ElementPos(:,5) = 0;  %angle

waveLength = (Resource.Parameters.speedOfSound/1000)/Trans.frequency;
if strcmp(Trans.units,'mm')
    Trans.ElementPosMm = Trans.ElementPos;
    Trans.ElementPosWL = Trans.ElementPos;
    Trans.ElementPosWL = Trans.ElementPos(:,1:3)./waveLength;
    Trans.spacing = Trans.spacingMm./waveLength;
    Trans.radius = Trans.radiusMm./waveLength;
else
    Trans.ElementPosMm = Trans.ElementPosMm;
    Trans.ElementPosMm = Trans.ElementPosMm(:,1:3)./waveLength;
    Trans.ElementPosWL = Trans.ElementPos;
end
eleWidthWl = Trans.elementWidth ./ waveLength;
Theta = linspace(-pi/4,pi/4,101);
Theta(51) = 0.0000001;
Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./...
    (eleWidthWl*pi*sin(Theta))));
% Specify PData structure array.
PData(1).PDelta = [0.5, 0, 0.5];
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.

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
Resource.ImageBuffer(1).numFrames = frame_num_img_yt;
z_plane = 2;
% XZ plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).Title = '3D Flash Image - XZ plane';
Resource.DisplayWindow(1).pdelta = 0.5;
Resource.DisplayWindow(1).Position = [0,480, ...
    ceil(PData(1).Size(2)*PData(1).PDelta(2)/Resource.DisplayWindow(1).pdelta), ... % width
    ceil(PData(1).Size(3)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta)];    % height
Resource.DisplayWindow(1).Orientation = 'xz';
% Resource.DisplayWindow(1).ReferencePt = [0,0.0,0.0];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0.0,PData(1).Origin(3)];
% Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0.0,-40];
Resource.DisplayWindow(1).Colormap = grayscaleCFImap;
Resource.DisplayWindow(1).splitPalette = 1;
Resource.DisplayWindow(1).AxesUnits = 'mm';
% Resource.DisplayWindow(1).mode = '2d';

% YZ plane
Resource.DisplayWindow(2).Type = 'Verasonics';
Resource.DisplayWindow(2).Title = '3D Flash Image - YZ plane';
Resource.DisplayWindow(2).pdelta = Resource.DisplayWindow(1).pdelta;
Resource.DisplayWindow(2).Position = [660,480, ...
    ceil(PData(1).Size(1)*PData(1).PDelta(1)/Resource.DisplayWindow(2).pdelta), ... % width
    ceil(PData(1).Size(3)*PData(1).PDelta(3)/Resource.DisplayWindow(2).pdelta)];    % height
Resource.DisplayWindow(2).Orientation = 'yz';
Resource.DisplayWindow(2).ReferencePt = [0,-PData(1).Origin(2),PData(1).Origin(3)];
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
    ceil(PData(1).Size(2)*PData(1).PDelta(2)/Resource.DisplayWindow(3).pdelta), ... % width
    ceil(PData(1).Size(1)*PData(1).PDelta(1)/Resource.DisplayWindow(3).pdelta)];    % height
Resource.DisplayWindow(3).Orientation = 'xy';
Resource.DisplayWindow(3).ReferencePt = [PData(1).Origin(1),...
    -PData(1).Origin(2),z_plane];%PData.Region(end).Shape.oPAIntersect];
Resource.DisplayWindow(3).Colormap = grayscaleCFImap;
Resource.DisplayWindow(3).splitPalette = 1;
Resource.DisplayWindow(3).AxesUnits = 'mm';

%% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,1,1]; % [Trans.frequency,.67,1,1]
TW(1).euqalize = 0;

TX_APOD = zeros(length(TX_ELEMS),Trans.numelements);
for itx = 1:na
    TX_APOD(itx,TX_ELEMS(itx)) = 1;
end
% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', kaiser(Resource.Parameters.numTransmit,1)', ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Trans.numelements)), 1, na);
% - Set event specific TX attributes.
if fix(na/2) == na/2       % if na even
    P.startAngle = (-(fix(na/2) - 1) - 0.5)*dtheta;
else
    P.startAngle = -fix(na/2)*dtheta;
end
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
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
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
pers = 20;
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
                         'compressFactor',40,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};

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
nsc = 5; % nsc is count of SeqControl objects

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:na                      % Acquire frame
        Event(n).info = 'Full aperture.';
        Event(n).tx = j;
        Event(n).rcv = na*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;
    end
    Event(n-1).seqControl = [3,nsc]; % modify last acquisition Event's seqControl
      SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
      nsc = nsc+1;

    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    if floor(i/5) == i/5     % Exit to Matlab every 5th frame
        Event(n).seqControl = 4;
    else
        Event(n).seqControl = 0;
    end
    n = n+1;
end

Event(n).info = 'Jump back';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 1;

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
