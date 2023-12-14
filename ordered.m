clear
close all
clc

%Import data: axial view
MRId = load('MRIdata.mat');
v_ax = MRId.vol;

figure
montage(v_ax)
title("MRI axial")

%Sagittal view
trasversalRot = flip(v_ax, 3);
trasversalRot = flip(trasversalRot, 1);
v_sag = permute(trasversalRot, [3 2 1]);

figure
montage(v_sag);
title('MRI sagittal')

%% Slide 135

%Select a ROI
rect = [137.5100   20.5100   42.9800   29.9800];
dimC = size(imcrop(v_sag(:,:,1),rect));
slice = 135;

roi_slice = imcrop(v_sag(:,:,slice),rect);
figure
subplot(1,3,1)
imshow(roi_slice)
str = sprintf('ROI of slice %.0f', slice);
title(str)

%Edge detection
edge_slice = edge(roi_slice, 'sobel');

subplot(1,3,2)
imshow(edge_slice)
str = sprintf('Edge of slice %.0f', slice);
title(str)

%Binary
ifill_slice = imfill(edge_slice, 'holes');

subplot(1,3,3)
imshow(ifill_slice)
str = sprintf('Binary of slice %.0f', slice);
title(str)

%Area
area_slice = bwarea(ifill_slice);

%% Sagittal volume

%Select a ROI
slf_sag = 150;
sli_sag = 112;
lenS = slf_sag-sli_sag+1;
sub_row = 6;
sub_col = 7;

for s = 1:lenS
    I = imcrop(v_sag(:,:,s+sli_sag-1),rect);
    roi(:,:,s) = I;
end

figure
montage(roi)
title("ROI original")

figure
for s = 1:lenS
    subplot(sub_row, sub_col, s)
    imhist(roi(:,:,s), 256)
end
sgtitle('Histogram ROI original')

%% Increase the contrast with average filter
dim_avg = 3;
avg_filt = (1/dim_avg^2).*ones(dim_avg, dim_avg);
roi_contrast = imfilter(roi, avg_filt, 'conv');

figure
montage(roi_contrast)
title('ROI after applying average filter')

figure
for s = 1:lenS
    subplot(sub_row, sub_col, s)
    imhist(roi_contrast(:,:,s), 256)
end
sgtitle('Histogram ROI after applying average filter')

%% Increase the contrast with median filter
%med remove salt and pepper noise better than avg
dim_med = 3;
for s = 1:lenS
    roi_contrast(:,:,s) = medfilt2(roi(:,:,s), [dim_med dim_med]);
end

figure
montage(roi_contrast)
title('ROI after applying median filter')

figure
for s = 1:lenS
    subplot(sub_row, sub_col, s)
    imhist(roi_contrast(:,:,s), 256)
end
sgtitle('Histogram ROI after applying median filter')

%% Increase the contrast adjusting image intensity
roi = im2double(roi);
max_roi = max(roi(:));
min_roi = min(roi(:));
LOW_in = min_roi;
HIGH_in = max_roi;
LOW_out = 0;
HIGH_out = 1;
gamma = 1.7; %increasing gamma, the image results more black

for s = 1:lenS
    roi_contrast(:,:,s) = imadjust(roi(:,:,s), [LOW_in HIGH_in], [LOW_out HIGH_out], gamma);
end

figure
montage(roi_contrast)
str = sprintf('ROI with contrast %.1f', gamma);
title(str)

figure
for s = 1:lenS
    subplot(sub_row, sub_col, s)
    imhist(roi_contrast(:,:,s), 256)
end
str = sprintf('Histogram ROI with contrast %.1f', gamma);
sgtitle(str)

%% Select a lesion using the threshold
disp("Threshold automatically set");
threshold = graythresh(roi_contrast);

disp("Threshold used in TH masking: "+string(threshold));

thMask = imbinarize(roi_contrast, threshold);
seD = strel('diamond',1);
thMask = imerode(thMask,seD);

figure
montage(thMask)
sgtitle('Binary image')


figure
for s = 1:lenS
    subplot(sub_row, sub_col, s)
    imshow(labeloverlay(roi_contrast(:,:,s), thMask(:,:,s)))
end
sgtitle('Check edges')

%volumeViewer(thMask)

%% Select a lesion using Edge detection
se90 = strel('line',3.4,90);
se0 = strel('line',3.4,0);
seD = strel('diamond',1);
%mask = zeros(dimC(1), dimC(2), lenS);
%edges = zeros(dimC(1), dimC(2), lenS);

figure
for s = 1:lenS
    edges(:, :, s) = edge(roi_contrast(:,:,s), 'sobel'); %log is better, but it doesn't find all
    subplot(sub_row, sub_col, s)
    imshow(edges(:,:,s))
end
sgtitle('Edges detection')

figure
for s = 1:lenS
    idil = imdilate(edges(:,:,s),[se90 se0]);
    ifill = imfill(idil);
    mf = imerode(ifill,seD);
    mf = imerode(mf,seD);
    mask(:,:,s) = mf;

    subplot(sub_row, sub_col, s)
    imshow(mask(:,:,s))
end
sgtitle('Binary image')

figure
for s = 1:lenS
    subplot(sub_row, sub_col, s)
    imshow(labeloverlay(roi_contrast(:,:,s),mask(:,:,s)))
end
sgtitle('Check edges')

%volumeViewer(mask)

%% Area of the tumor
area= zeros(lenS,1);
for s = 1:lenS
    area(s,1) = bwarea(thMask(:,:,s));
end

vol_tumor = sum(area);

%% Axial view
clear roi
clear roi_contrast
clear thMask

%Select ROI
rect = [133.5100  103.5100   51.9800   44.9800];
slf_ax = 90;
sli_ax = 65;
lenS = slf_ax-sli_ax+1;
sub_row = 6;
sub_col = 5;

for s = 1:lenS
    I = imcrop(v_ax(:,:,s+sli_ax-1),rect);
    roi(:,:,s) = I;
end

figure
montage(roi)
title("ROI original")

figure
for s = 1:lenS
    subplot(sub_row, sub_col, s)
    imhist(roi(:,:,s), 256)
end
sgtitle('Histogram ROI original')

%% Increase the contrast with median filter
dim_med = 3;
for s = 1:lenS
    roi_contrast(:,:,s) = medfilt2(roi(:,:,s), [dim_med dim_med]);
end

figure
montage(roi_contrast)
title('ROI after applying median filter')

figure
for s = 1:lenS
    subplot(sub_row, sub_col, s)
    imhist(roi_contrast(:,:,s), 256)
end
sgtitle('Histogram ROI after applying median filter')

%% Increase the contrast adjusting image intensity
roi = im2double(roi);
max_roi = max(roi(:));
min_roi = min(roi(:));
LOW_in = min_roi;
HIGH_in = max_roi;
LOW_out = 0;
HIGH_out = 1;
gamma = 1.5; %increasing gamma, the image results more black

for s = 1:lenS
    roi_contrast(:,:,s) = imadjust(roi(:,:,s), [LOW_in HIGH_in], [LOW_out HIGH_out], gamma);
end

figure
montage(roi_contrast)
str = sprintf('ROI with contrast %.1f', gamma);
title(str)

figure
for s = 1:lenS
    subplot(sub_row, sub_col, s)
    imhist(roi_contrast(:,:,s), 256)
end
str = sprintf('Histogram ROI with contrast %.1f', gamma);
sgtitle(str)

%% Select a lesion using the threshold
disp("Threshold automatically set");
threshold = graythresh(roi_contrast);

disp("Threshold used in TH masking: "+string(threshold));

thMask = imbinarize(roi_contrast, threshold);
seD = strel('diamond',1);
thMask = imerode(thMask,seD);

figure
montage(thMask)
sgtitle('Binary image')

figure
for s = 1:lenS
    subplot(sub_row, sub_col, s)
    imshow(labeloverlay(roi_contrast(:,:,s), thMask(:,:,s)))
end
sgtitle('Check edges')

%volumeViewer(thMask)

%% Select a lesion using Edge detection
se90 = strel('line',3.4,90);
se0 = strel('line',3.4,0);
seD = strel('diamond',1);

figure
for s = 1:lenS
    edges(:, :, s) = edge(roi_contrast(:,:,s), 'sobel');
    subplot(sub_row, sub_col, s)
    imshow(edges(:,:,s))
end
sgtitle('Edges detection')

figure
for s = 1:lenS
    idil = imdilate(edges(:,:,s),[se90 se0]);
    ifill = imfill(idil);
    mf = imerode(ifill,seD);
    mf = imerode(mf,seD);
    mask(:,:,s) = mf;

    subplot(sub_row, sub_col, s)
    imshow(mask(:,:,s))
end
sgtitle('Binary image')

figure
for s = 1:lenS
    subplot(sub_row, sub_col, s)
    imshow(labeloverlay(roi_contrast(:,:,s),mask(:,:,s)))
end
sgtitle('Check edges')

%volumeViewer(mask)
