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
    
    % roi = histeq(roi);
    % figure
    % montage(roi)
    % title("ROI")

    % for s=1:lenS
    %     roi(:,:,s) = imadjust(roi(:,:,s));
    % end
    % figure
    % montage(roi)
    % title("ROI")
    % figure
    % for s = 1:lenS
    %     subplot(6, 7, s)
    %     imhist(roi(:,:,s), 256)
    % end
    % sgtitle('Histogram ROI original')
    
    %%
    se90 = strel('line',3,90);
    se0 = strel('line',3,0);
    seD = strel('diamond',1);
    h = [0 -1 0; -1 4 -1; 0 -1 0];
    edges2 = zeros(dimC(1), dimC(2), slf-sli+1);
    roi2 = roi/2;
    for s = 1:lenS
        
        for i=1:dimC(1)
            for j=1:dimC(2)
                if roi2(i,j,s) > 75
                    roi2(i,j,s) = roi2(i,j,s)+127;
                    if(roi2(i,j,s) > 230)
                        roi2(i,j,s) = 0;
                    end
                else
                    %roi2(i,j,s) = 0;
                end
            end
        end
     
    

        

        
        
    end
    
    
    
    figure
    for s = 1:lenS
        subplot(6,7,s)
        imshow(roi2(:,:,s))
    end
    
    roi3 = roi2;
    dim_med = 3;
    for s = 1:lenS
        roi3(:,:,s) = medfilt2(roi3(:,:,s), [dim_med dim_med]);
    end
    figure
    for s = 1:lenS
        subplot(6,7,s)
        imshow(roi3(:,:,s))
    end
    sgtitle("Solo avg");


    figure
    for s = 1:lenS
        roi2(:,:,s) = imfill(roi2(:,:,s));
        subplot(6,7,s)
        imshow(roi2(:,:,s))
    end


    dim_med = 3;
    for s = 1:lenS
        roi2(:,:,s) = medfilt2(roi2(:,:,s), [dim_med dim_med]);
    end

    figure
    for s = 1:lenS
        subplot(6,7,s)
        imshow(roi2(:,:,s))
    end
    figure
    for s=1:lenS
        edges2(:,:,s) = edge(roi2(:,:,s) , 'log');
        
        edges2(:,:,s) = imdilate(edges2(:,:,s),[se90 se0]);
        subplot(6,7,s)
        imshow(edges2(:,:,s))
        
        ifill = imfill(edges2(:,:,s));
        %mf = imerode(ifill,seD);
        edges2(:,:,s) = imerode(ifill,seD);

        
    end
   



    %% Using Edge detection
    
    edges = zeros(dimC(1), dimC(2), slf-sli+1);
    for s = 1:lenS
        edges(:, :, s) = edge(roi(:,:,s) , 'canny'); %'log');
        % figure
        % imshow(edges(:,:,s))
    end
    
    
    se90 = strel('line',3.4,90);
    se0 = strel('line',3.4,0);
    seD = strel('diamond',1);
    
    mask = zeros(dimC(1), dimC(2), slf-sli+1);
    
    for s = 1:lenS
        idil = imdilate(edges2(:,:,s),[se90 se0]);
        
        ifill = imfill(idil);
        % figure
        % imshow(ifill)
        mf = imerode(ifill,seD);
        mf = imerode(mf,seD);
    
        mask(:,:,s) = mf;
    
        % figure
        % imshow(mask(:,:,s))
    end 
    figure
    for s = 1:lenS
        subplot(6,7,s)
        imshow(labeloverlay(roi2(:,:,s),mask(:,:,s)))
    end
    %% Using threshold
    % Image adjustment is important here, in edges about nothing changes
    % For automatic threshold is better without imadjust to consider
    % less of the fibrotic part of the brain
    % For manual threshold imadjust is useful in same cases. Manual
    % threshold suggested with imadjust is 150
    roi = roi2; %TODO
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
    figure
    for s = 1:lenS
        subplot(6,7,s)
        imshow(thMask(:,:,s))
    end
    maskE = mask;
    maskTh = thMask;

end

