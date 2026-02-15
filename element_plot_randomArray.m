clear all; clc; %close all;

load ("C:\Users\daiwe\Dropbox\CavitationMapping_2024_Davia_Chaorui\20260204_RandomArray_ele\run001_eval0000774_centers.mat")

x_trans = elem_centers_mm(:,1);
y_trans = elem_centers_mm(:,2);
z_trans = zeros(256,1);


M = reshape(x_trans, 16, 16);  
N = reshape(y_trans, 16, 16);   

M =rot90(M,-1);  
N =rot90(N,-1);

x_trans = reshape(M, [], 1);
y_trans = reshape(N, [], 1);

figure, hold on
for ielem = 1:length(x_trans)
    scatter(x_trans(ielem),y_trans(ielem));
    text(x_trans(ielem),y_trans(ielem),num2str(ielem));
end
