savedir = 'C:\Davia\20260213_random\1_USphantom\';
checkMakeDir(savedir)
target = 't1_USphantom';
% timeSA = '_';
% savename = [savedir,target,timeSA,'_US.mat'];
savename = [savedir,target,'_US.mat'];
tic
save(savename,'mediumTemp','RcvData','Trans','P','Receive','TW','TX',...
    'SeqControl','TX_ELEMS','USVolt_yt','-v7.3')
toc
