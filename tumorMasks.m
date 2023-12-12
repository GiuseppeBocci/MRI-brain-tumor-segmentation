function [maskE, maskTh] = tumorMasks(vol, rect, sli, slf, ngraybits, varargin)
    dimC = size(imcrop(vol(:,:,1),rect));
    
    roi = zeros(dimC(1), dimC(2), slf-sli+1, "uint8");
    lenS = slf-sli+1;
    for s = 1:lenS

        I = imcrop(vol(:,:,s+sli-1),rect);

        % figure
        % subplot(1,3,1)
        % imshow(I)
        % subplot(1,3,3)
        % bins = 2^ngraybits;
        % histogram(I,bins)
        % h = histcounts(I, bins);
        % if s>lenS*0.1 && s<lenS*0.9 
        %     I = imadjust(I);
        % 
        %     subplot(1,3,2)
        %     imshow(I)
        % end

        roi(:,:,s) = I;
        
        %figure
        %imshow(roi(:,:,s))
    end

    figure
    montage(roi)
    title("ROI")
    
    %% Using Edge detection
    
    edges = zeros(dimC(1), dimC(2), slf-sli+1);
    for s = 1:lenS
        edges(:, :, s) = edge(roi(:,:,s) , 'log');
        % figure
        % imshow(edges(:,:,s))
    end
    
    
    se90 = strel('line',3.4,90);
    se0 = strel('line',3.4,0);
    seD = strel('diamond',1);
    
    mask = zeros(dimC(1), dimC(2), slf-sli+1);
    
    for s = 1:lenS
        idil = imdilate(edges(:,:,s),[se90 se0]);
        
        ifill = imfill(idil);
        % figure
        % imshow(ifill)
        mf = imerode(ifill,seD);
        mf = imerode(mf,seD);
    
        mask(:,:,s) = mf;
    
        % figure
        % imshow(mask(:,:,s))
    end 
    
    % for s = 1:lenS
    %     figure
    %     imshow(labeloverlay(roi(:,:,s),mask(:,:,s)))
    % end
    
    %% Using threshold
    % Image adjustment is important here, in edges about nothing changes
    % For automatic threshold is better without imadjust to consider
    % less of the fibrotic part of the brain
    % For manual threshold imadjust is useful in same cases. Manual
    % threshold suggested with imadjust is 150

    if length(varargin) == 1 
       threshold = cell2mat(varargin(1))/(2^ngraybits-1);
    else
        disp("Threshold automatically set");
        threshold = graythresh(roi);
    end
    disp("Threshold used in TH masking: "+string(threshold));
    
    thMask = imbinarize(roi, threshold);
    
    % thMask = 255*thMask; % imfill needs a numerical matrix
    % thMask= imfill(thMask); improves manual threshold

    % for s = 1:lenS
    %     % figure
    %     % imshow(labeloverlay(roi(:,:,s),thMask(:,:,s)))
    % end

    %%
    maskE = mask;
    maskTh = thMask;

end

