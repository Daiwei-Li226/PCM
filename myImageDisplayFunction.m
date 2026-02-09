function myImageDisplayFunction (ImgData)
    persistent myHandle1
    persistent myHandle2
    persistent myHandle3

 A =ImgData(:,:,:,1);
 XY = squeeze(max(abs(A),[],3));
 XZ = squeeze(max(abs(A),[],2));
 YZ = squeeze(max(abs(A),[],1));
 XY = rot90(XY, 3);
 XZ = rot90(XZ, 3);
 YZ = rot90(YZ, 3);
    
    % --- 3. 初始化和显示 ---
    if isempty(myHandle1) || ~ishandle(myHandle1.figure)
        myHandle1.figure = figure('Name', 'Custom XY Plane Display');
        myHandle1.axes = axes(myHandle1.figure);
    end
        if isempty(myHandle2) || ~ishandle(myHandle2.figure)
        myHandle2.figure = figure('Name', 'Custom XZ Plane Display');
        myHandle2.axes = axes(myHandle2.figure);
        end
            if isempty(myHandle3) || ~ishandle(myHandle3.figure)
        myHandle3.figure = figure('Name', 'Custom YZ Plane Display');
        myHandle3.axes = axes(myHandle3.figure);
    end
    
    imagesc(myHandle1.axes, log(XY));
    colormap(myHandle1.axes, "parula");
    drawnow 
        imagesc(myHandle2.axes, log(XZ));
    colormap(myHandle2.axes, "parula");
    drawnow 
        imagesc(myHandle3.axes, log(YZ));
    colormap(myHandle3.axes, "parula");
    drawnow 
    
end