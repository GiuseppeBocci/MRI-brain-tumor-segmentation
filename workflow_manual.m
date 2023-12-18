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

%% Select a ROI
rect = [137.5100   20.5100   42.9800   29.9800];
dimC = size(imcrop(v_sag(:,:,1),rect));
slice = 256-135; %Because we flipped the slices

roi_slice = imcrop(v_sag(:,:,slice), rect);
figure
subplot(2,2,1)
imshow(roi_slice)
str = sprintf('ROI of slice %.0f', slice);
title(str)

%% Automatic segmentation of slice 135
%Edge detection
edge_slice = edge(roi_slice, 'canny');

subplot(2,2,2)
imshow(edge_slice)
str = sprintf('Edge of slice %.0f', slice);
title(str)

%Binary
ifill_slice = imfill(edge_slice, 'holes');

subplot(2,2,3)
imshow(ifill_slice)
str = sprintf('Binary of slice %.0f', slice);
title(str)

%Area
area_automatic = bwarea(ifill_slice);

%% Manually segmentation of slice 135
imageSegmenter(roi_slice);

%% Plot of manually segmentation
maskedImage = logical(maskedImage);

subplot(2,2,4)
imshow(maskedImage)
str = sprintf('Manual binary of slice %.0f', slice);
title(str)

%Area
area_manual = bwarea(maskedImage);

%% Segmentation performances
for j = 1:dimC(1)
    for i = 1: dimC(2)
        TP(j,i) = ((ifill_slice(j,i)==1) & (maskedImage(j,i)==1));
        TN(j,i) = ((ifill_slice(j,i)==0) & (maskedImage(j,i)==0));
        FP(j,i) = ((ifill_slice(j,i)==1) & (maskedImage(j,i)==0));
        FN(j,i) = ((ifill_slice(j,i)==0) & (maskedImage(j,i)==1));
    end
end

TP = sum(sum(TP));
TN = sum(sum(TN));
FP = sum(sum(FP));
FN = sum(sum(FN));

%Sensitivity
sen = TP/(TP+FN);

%Specificity
spec = TN/(TN+FP);

%Dice coefficient
dice = 2*TP/(area_automatic+area_manual);

