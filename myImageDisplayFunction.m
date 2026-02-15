function myImageDisplayFunction(ImgData)
    % 存储图形句柄
    persistent hImg1 hImg2 hImg3
    % 存储物理参数，避免重复 evalin
    persistent x_mm y_mm z_mm initialized

    % --- 1. 仅在第一次运行时获取物理参数 ---
    if isempty(initialized)
        PData = evalin('base', 'PData(1)');
        Trans = evalin('base', 'Trans');
        wl = Trans.wavelength;
        
        y_wl = PData.Origin(2) + (0:PData.Size(1)-1) * PData.PDelta(2);
        x_wl = PData.Origin(1) + (0:PData.Size(2)-1) * PData.PDelta(1);
        z_wl = PData.Origin(3) + (0:PData.Size(3)-1) * PData.PDelta(3);
        
        y_mm = y_wl * wl;
        x_mm = x_wl * wl;
        z_mm = z_wl * wl;
        initialized = true;
    end

    % --- 2. 提取数据并投影 ---
    A = abs(ImgData(:,:,:,1)); % 先取绝对值
    
    % 执行 MAP (这里可以用 max(A, [], dim) 替代 max(abs(A)))
    XY = rot90(log(squeeze(max(A, [], 3)) + eps), 3);
    XZ = rot90(log(squeeze(max(A, [], 2)) + eps), 3);
    YZ = rot90(log(squeeze(max(A, [], 1)) + eps), 3);

    % --- 3. 初始化或更新显示 ---
    if isempty(hImg1) || ~ishandle(hImg1)
        % --- XY Plane ---
        f1 = figure('Name', 'XY Plane');
        ax1 = axes(f1);
        hImg1 = image(ax1, y_mm, x_mm, XY, 'CDataMapping', 'scaled');
        colormap(ax1, 'gray'); axis(ax1, 'image');
        xlabel(ax1, 'Y [mm]'); ylabel(ax1, 'X [mm]');
        
        % --- XZ Plane ---
        f2 = figure('Name', 'XZ Plane');
        ax2 = axes(f2);
        hImg2 = image(ax2, y_mm, z_mm, XZ, 'CDataMapping', 'scaled');
        colormap(ax2, 'gray'); axis(ax2, 'image');
        xlabel(ax2, 'X [mm]'); ylabel(ax2, 'Z [mm]');
        
        % --- YZ Plane ---
        f3 = figure('Name', 'YZ Plane');
        ax3 = axes(f3);
        hImg3 = image(ax3, x_mm, z_mm, YZ, 'CDataMapping', 'scaled');
        colormap(ax3, 'gray'); axis(ax3, 'image');
        xlabel(ax3, 'Y [mm]'); ylabel(ax3, 'Z [mm]');
    else
        % --- 关键提速：直接更新数据属性 ---
        set(hImg1, 'CData', XY);
        set(hImg2, 'CData', XZ);
        set(hImg3, 'CData', YZ);
    end

    drawnow limitrate
end





% function myImageDisplayFunction (ImgData)
%     persistent myHandle1
%     persistent myHandle2
%     persistent myHandle3
% 
% PData = evalin('base', 'PData(1)');
% Trans = evalin('base', 'Trans');
% wl = Trans.wavelength; % mm
% 
%  A =ImgData(:,:,:,1);
%  XY = squeeze(max(abs(A),[],3));
%  XZ = squeeze(max(abs(A),[],2));
%  YZ = squeeze(max(abs(A),[],1));
% 
% y_wl = PData.Origin(2) + (0:PData.Size(1)-1) * PData.PDelta(2);
% x_wl = PData.Origin(1) + (0:PData.Size(2)-1) * PData.PDelta(1);
% z_wl = PData.Origin(3) + (0:PData.Size(3)-1) * PData.PDelta(3);
% y_mm = y_wl * wl;
% x_mm = x_wl * wl;
% z_mm = z_wl * wl;
% 
% XY_disp = XY;
% XZ_disp = XZ;
% YZ_disp = YZ;
% 
% XY_disp = rot90(log(XY + eps), 3);
% XZ_disp = rot90(log(XZ + eps), 3);
% YZ_disp = rot90(log(YZ + eps), 3);
% 
%     
%     % --- 3. 初始化和显示 ---
%     if isempty(myHandle1) || ~ishandle(myHandle1.figure)
%         myHandle1.figure = figure('Name', 'Custom XY Plane Display');
%         myHandle1.axes = axes(myHandle1.figure);
%         xlabel(myHandle1.axes, 'Y [mm]');
%         ylabel(myHandle1.axes, 'X [mm]');
%     end
%         if isempty(myHandle2) || ~ishandle(myHandle2.figure)
%         myHandle2.figure = figure('Name', 'Custom XZ Plane Display');
%         myHandle2.axes = axes(myHandle2.figure);
%         xlabel(myHandle2.axes, 'X [mm]');
%         ylabel(myHandle2.axes, 'Z [mm]');
%         end
%             if isempty(myHandle3) || ~ishandle(myHandle3.figure)
%         myHandle3.figure = figure('Name', 'Custom YZ Plane Display');
%         myHandle3.axes = axes(myHandle3.figure);
%         xlabel(myHandle3.axes, 'Y [mm]');
%         ylabel(myHandle3.axes, 'Z [mm]');
%     end
%     
%     imagesc(myHandle1.axes, y_mm, x_mm, XY_disp);
%     colormap(myHandle1.axes, "gray");
%     axis(myHandle1.axes, 'image'); % 保持 1:1 比例，防止图像变形
%     
%     % XZ Display: 横轴 Z, 纵轴 Y
%     imagesc(myHandle2.axes, y_mm,z_mm,  XZ_disp);
%     colormap(myHandle2.axes, "gray");
%     axis(myHandle2.axes, 'image');
%     
%     % YZ Display: 横轴 Z, 纵轴 X
%     imagesc(myHandle3.axes,  x_mm, z_mm,YZ_disp);
%     colormap(myHandle3.axes, "gray");
%     axis(myHandle3.axes, 'image');
% 
%     drawnow limitrate
%     
% end