function Trans = computeTrans_256_2D_Martix(Trans,sos)

KnownTransducers = {'custom2561','0000'...
    'custom2562','0000'...
    'qcr2561','0000'...
    'qcr2562','0000'};

%Transducer paramter
frequency =5.6; % MHz
pitch = 2;     % mm
kerf =1;      % mm

if ~isstruct(Trans)  % if a structure is not provided as input, assume input is ID to translate into string.
    n = find(strcmpi(Trans, KnownTransducers));
    if isempty(n), Trans = 'Unknown';
    else Trans = KnownTransducers{n(1)-1};
    end
    return
end

if ~isfield(Trans,'name'), error('computeTrans: Trans.name must be provided in input structure.\n'); end
n = find(strcmpi(Trans.name, KnownTransducers));
if isempty(n), error('computeTrans: Trans.name not recognized as known transducer.'); end
speedOfSound = 1.540;  % default speed of sound in mm/usec
if evalin('base','exist(''Resource'',''var'')&&isfield(Resource,''Parameters'')')
    if evalin('base','isfield(Resource.Parameters,''speedOfSound'')')
        speedOfSound = evalin('base','Resource.Parameters.speedOfSound')/1000; % speed of sound in mm/usec
    end
end

    switch KnownTransducers{n}
          
         case 'qcr2561'
            Trans.frequency = frequency;
            Trans.type = 2;     % Array geometry is linear (x values only).
            Trans.id = hex2dec('0000');
            Trans.connType = 2; % HDI connector
            Trans.numelements = 128;
            Trans.elementWidth = pitch-kerf;   %mm
            Trans.spacingMm = pitch;           %mm
            waveLength=(sos/1000)/Trans.frequency;
            Trans.spacing   = Trans.spacingMm / waveLength;
            eleWidthWl = Trans.elementWidth /waveLength;
            Trans.ElementPos = zeros(Trans.numelements,5);
         

            x_temp = (([1:16]-8.5))';
            x = repmat(x_temp,8,1);
            y = [];
            for i = 1:8
                y_temp = repmat(((8.5-i)),16,1);
                y = [y;y_temp];
            end
            z = zeros(128,1);
            temp = [x,y,z];
            Trans.ElementPosMm(:,1:3) = temp * Trans.spacingMm;
            Trans.ElementPos(:,1:3) = temp * Trans.spacing;
            
            if ~isfield(Trans,'ElementSens')
                % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                Theta = (-pi/2:pi/100:pi/2);
                Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
            end
            Trans.impedance = 50; % using default value 
            
            if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 50; end
        case 'qcr2562'
            if ~isfield(Trans,'frequency'), Trans.frequency = frequency; end
            Trans.type = 2;     % Array geometry is linear (x values only).
            Trans.id = hex2dec('0000');
            Trans.connType = 2; % HDI connector
            Trans.numelements = 128;
            Trans.elementWidth = pitch-kerf;   %mm
            Trans.spacingMm = pitch;      %mm
            waveLength=(sos/1000)/Trans.frequency;
            Trans.spacing   = Trans.spacingMm / waveLength;
            eleWidthWl = Trans.elementWidth /waveLength; 
            Trans.ElementPos = zeros(Trans.numelements,5);

            x_temp = (([1:16]-8.5))';
            x = repmat(x_temp,8,1);
            y = [];
            for i = 1:8
                y_temp = -repmat(((i-0.5)),16,1);
                y = [y;y_temp];
            end
            z = zeros(128,1);
            temp = [x,y,z];
            Trans.ElementPosMm(:,1:3) = temp * Trans.spacingMm;
            Trans.ElementPos(:,1:3) = temp * Trans.spacing;
            Trans.frequency=frequency;

            if ~isfield(Trans,'ElementSens')
                % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                Theta = (-pi/2:pi/100:pi/2);
                Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
            end
            Trans.impedance = 50; % using default value 
            if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 50; end
    end
end