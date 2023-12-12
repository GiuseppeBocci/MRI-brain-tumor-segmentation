close all
clear all

MRId = load('MRIdata.mat');

%slice = 77 ; % (90 esculsa mail tuomere termina tra 89e90) 89 - 64 (64 non si vede molto potrebbe essere esculsa)
v = MRId.vol;
%volumeViewer(v);

% for s = 1:90 %MRId.dim(3)
%     figure
%     imshow(v(:,:,s))
% end

figure
montage(v)
title("MRI transversal")

%% ROI T
close all
%I = v(:,:,slice);
%[cropped, rect] = imcrop(I);
rect = [133.5100  103.5100   51.9800   44.9800];

[mask, thMask] = tumorMasks(v, rect, 65, 90, 8);

%% VOLUMES T
close all

pixelVol = MRId.pixdim(1)*MRId.pixdim(2)*MRId.pixdim(3) % mm^3

[areaE, areaT] = areasFromMasks(mask, thMask);

volE = pixelVol*sum(areaE)
volT = pixelVol*sum(areaT)
volumeViewer(mask)
volumeViewer(thMask)

%% SAGITTAL
trasversalRot = flip(v, 3); % order of sicels
% from rx to sx and not the inverse. No computational relevant
trasversalRot = flip(trasversalRot, 1);
sagittal = permute(trasversalRot, [3 2 1]);
figure
montage(sagittal)
title("MRI sagittal")
% slices 112(113) - 152(150)

%% 
close all
volumeViewer close

%I = sagittal(:,:,135);
%[cropped, rect] = imcrop(I);
rect = [137.5100   20.5100   42.9800   29.9800];

[mask, thMask] = tumorMasks(sagittal, rect, 113, 150, 8); % threshold);

[areaE, areaT] = areasFromMasks(mask, thMask);

volE = pixelVol*sum(areaE)
volT = pixelVol*sum(areaT)
volumeViewer(mask)
volumeViewer(thMask)

%% CORONAL
trasversalRot = flip(v, 3); % order of sicels
% from front to back and not the inverse. No computational relevant
trasversalRot = flip(trasversalRot, 2);
coronal = permute(trasversalRot, [3 1 2]);
figure
montage(coronal)
title("MRI coronal")
% slices 81 - 116(115)

%% 
close all
volumeViewer close

% I = coronal(:,:,100);
% [cropped, rect] = imcrop(I);
rect = [104.5100   20.5100   42.9800   31.9800];

[mask, thMask] = tumorMasks(coronal, rect, 81, 116, 8); % threshold);

[areaE, areaT] = areasFromMasks(mask, thMask);

volE = pixelVol*sum(areaE)
volT = pixelVol*sum(areaT)
volumeViewer(mask)
volumeViewer(thMask)