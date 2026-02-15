
us_rec_total = us_tatol_DMAS;

% xymap = squeeze(max(abs(mean(us_rec(:,:,:,:),4)),[],3))';
xymap = squeeze(max(abs(us_rec_total),[],3))';
xzmap = squeeze(max(abs(us_rec_total),[],1))';
yzmap = squeeze(max(abs(us_rec_total),[],2))';

fig = figure();
subplot(1,3,1),imagesc(x_img,y_img,20*log10(xymap))
colormap gray
pbaspect([length(x_img)/length(y_img) 1 1]),xlabel('x [mm]'),ylabel('y [mm]'), title('DMAS XY')
subplot(1,3,2),imagesc(x_img,z_img,20*log10(xzmap))
% subplot(1,3,2),imagesc(x_img,z_img,xzmap)
colormap gray
pbaspect([length(y_img)/length(z_img) 1 1]),xlabel('y [mm]'),ylabel('z [mm]'), title('DMAS YZ')
subplot(1,3,3),imagesc(y_img,z_img,20*log10(yzmap))
colormap gray
pbaspect([length(x_img)/length(z_img) 1 1]),xlabel('x [mm]'),ylabel('z [mm]'), title('DMAS XZ')
drawnow

%%
us_rec = mean(us_rec(:,:,:,:),4);
%%
Res = 0.12*2;
x_range = [-20,20];
y_range = [-20,20];
z_range = [40,130];
x_img = x_range(1):Res:x_range(2);
y_img = y_range(1):Res:y_range(2);
z_img = z_range(1):Res:z_range(2);


%%

saving_path = 'C:\Users\daiwe\OneDrive - Duke University\02_PI LAB\02_3DPCM_Matrixarray\2025_IEEETMI\figure\fig7\';
save([saving_path,'us_phantom_DMAS_recon.mat'],'us_total_DMAS','x_img','y_img','z_img','Res');
