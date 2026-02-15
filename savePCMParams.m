temperature_C = mediumTemp;
fxngen_delay = 100; %us
laser_PRF = 10;
filename = [targetfolder,'_setupParams.mat'];
% filename = [dataname,'_setupParams.mat'];
filedir = fdir;
checkMakeDir(filedir)
RFDataSize = [Resource.RcvBuffer.rowsPerFrame, Resource.RcvBuffer.colsPerFrame, Resource.RcvBuffer.numFrames];
save([filedir,filename],...
    'PData','RFDataSize','fxngen_delay','laser_PRF',...
    'Receive','TX','Trans','temperature_C');