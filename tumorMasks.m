function [maskE, maskTh] = tumorMasks(v, rect, sli, slf)
    dimC = size(imcrop(v(:,:,1),rect));
    
    roi = zeros(dimC(1), dimC(2), slf-sli+1, "uint8");
    lenS = slf-sli+1;
    for s = 1:lenS
        % figure
        I = imadjust(imcrop(v(:,:,s+sli-1),rect));
        % imshow(I)
        % figure
        % imshow(imadjust(I))
        roi(:,:,s) = I;
        
        %figure
        %imshow(roi(:,:,s))
    end
    
    %% Using Edge detection
    
    close all
    
    edges = zeros(dimC(1), dimC(2), slf-sli+1);
    for s = 1:lenS
        edges(:, :, s) = edge(roi(:,:,s) , 'log');
    end
    
    
    se90 = strel('line',3,90);
    se0 = strel('line',3,0);
    seD = strel('diamond',1);
    
    mask = zeros(dimC(1), dimC(2), slf-sli+1);
    
    for s = 1:lenS
        idil = imdilate(edges(:,:,s),[se90 se0]);
        ifill = imfill(idil);
        mf = imerode(ifill,seD);
        mf = imerode(mf,seD);
    
        mask(:,:,s) = mf;
    
        %figure
        %imshow(mask(:,:,s))
    end 
    
    % for s = 1:lenS
    %     figure
    %     imshow(labeloverlay(roi(:,:,s),mask(:,:,s)))
    % end
    
    %% Using just threshold
    % Image adjustment is important here, in edges about nothing changes
    
    close all
    
    threshold = 130;
    disp("Threshold used in TH masking: "+string(threshold));
    
    I = roi(:,:,10);
    figure
    subplot(1,2,1)
    imshow(I)
    subplot(1,2,2)
    histogram(I)
    
    seD = strel('diamond',1);
    
    thMask = zeros(dimC(1), dimC(2), slf-sli+1, "double");
    for s = 1:lenS
        I = roi(:,:,s);
        for i = 1:dimC(1)
            for j = 1:dimC(2)
                thMask(i,j, s) = I(i,j) >= threshold;    
            end
        end
        
        % figure
        % imshow(255*I)
    end
    
    thMask = 255*thMask;
    for s = 1:lenS
        thMask(:,:,s) = imerode(imfill(thMask(:,:,s)), seD);
        % figure
        % imshow(labeloverlay(roi(:,:,s),thMask(:,:,s)))
    end
    
    
    % for s = sli:slf
    %     figure
    %     imshow(logical(roi(:,:,s)))
    % end
    %%
    maskE = mask;
    maskTh = thMask;

    %%
    % figure
    % montage(roi, colormap("gray"));
    % hold on
    % montage(mask, [0 0 1; 0 1 0]);
    % drawnow;
end

