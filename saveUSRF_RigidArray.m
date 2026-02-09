savedir = 'C:\Davia\20260204_Random\3_USphantom\';
% savedir = 'E:\Nanchao Wang\20240728\PulseEcho\Softarray_US\';
checkMakeDir(savedir)
target = 't1';
% target = 'T1-Matrix-RePos-Reconnect-100-GraphiteCross-1MHz-10frame-90v';
timeSA = '_';

savename = [savedir,target,timeSA,'_US.mat'];
tic
save(savename,'mediumTemp','RcvData','Trans','P','Receive','TW','TX',...
    'SeqControl','TX_ELEMS','USVolt_yt','-v7.3')
toc
