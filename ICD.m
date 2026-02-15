% 参数
fs = 19.2308e6;  c = 1500;         % Hz, m/s
fc = 5.6e6;  fbw = 0.75;
bp = [fc*(1-fbw/2), fc*(1+fbw/2)]; % [3.5 7.7] MHz
center = [5, 1, 24.5];             % mm
radius = 3;  step = 0.5;           % mm
win_s = 5e-7; hop_s = 5e-7;

% 1) 读取+裁剪 -> (样本, 通道) = (16896,256) → 取 10000:11100
% fid = fopen('_0020.dat','rb');
% A = fread(fid, [16896,256], 'double'); fclose(fid);
% A = A(10000:11100, 1:256);          % (1101 x 256)
% W = A.';                             % (256 x 1101)

f = fopen('_0020.dat','r');
iRFData = fread(f,'int16');
iRFData = reshape(iRFData,[32000,256]);
fclose(f);
RFData = iRFData(10100:10400,:);
 W = RFData .';   

% 2) 阵元坐标（mm→m）
xmm = load('x.mat'); fn = fieldnames(xmm); xmm = xmm.(fn{1})(:);
ymm = load('y.mat'); fn = fieldnames(ymm); ymm = ymm.(fn{1})(:);
zmm = load('z.mat'); fn = fieldnames(zmm); zmm = zmm.(fn{1})(:);
pos = [xmm(1:256), ymm(1:256), zmm(1:256)] * 1e-3;  % m

% 3) 带通滤波
d = designfilt('bandpassiir','FilterOrder',8, ...
    'HalfPowerFrequency1',bp(1),'HalfPowerFrequency2',bp(2), ...
    'SampleRate',fs);
Wf = filtfilt(d, W.').';             % (256 x 1101)

% 4) 体素与掩膜
[xg,yg,zg] = ndgrid(center(1)-radius:step:center(1)+radius, ...
                    center(2)-radius:step:center(2)+radius, ...
                    center(3)-radius:step:center(3)+radius);
vox = [xg(:), yg(:), zg(:)];
mask = vecnorm(vox - center, 2, 2) <= radius;
vox_m = vox(mask,:) * 1e-3;          % m
dV = (step*1e-3)^3;

% 5) 滑动时间窗
L = round(win_s*fs); H = round(hop_s*fs);
idx = 1:H:(size(Wf,2)-L+1);

% 6) DAS + ICD
ICD = zeros(numel(idx),1);

for k = 1:numel(idx)
    seg = Wf(:, idx(k):idx(k)+L-1);  % (C x L)
    % 延时（矢量化计算）
    d = vecnorm(pos - permute(vox_m,[3 2 1]), 2, 2);   % (C x V)
    tau = d / c;  idxs = round(-tau*fs);               % (C x V)
    % 对每个体素做对齐叠加
    pow = zeros(size(vox_m,1),1);
    for v = 1:size(vox_m,1)
        sh = idxs(:,v);
        s0 = max(0, -min(sh)); e0 = min(L, L - max(sh));
        if e0 - s0 <= 8, continue; end
        acc = zeros(1, e0-s0);
        for cch = 1:size(seg,1)
            acc = acc + seg(cch, (s0+sh(cch)+1):(e0+sh(cch)));
        end
        pow(v) = mean(acc.^2);
    end
    ICD(k) = sum(pow) * dV * win_s;
end
disp(ICD)
