close all
clear all

MRId = load('MRIdata.mat');

slice = 77 ; % (90 esculsa mail tuomere termina tra 89e90) 89 - 64 (64 non si vede molto potrebbe essere esculsa)
v = MRId.vol;
%volumeViewer(v);

% for s = 1:90 %MRId.dim(3)
%     figure
%     imshow(v(:,:,s))
% end

%% ROI
close all
%I = v(:,:,slice);
%[cropped, rect] = imcrop(I);
rect = [133.5100  103.5100   51.9800   44.9800];

[mask, thMask] = tumorMasks(v, rect, 65, 90);

%% VOLUMES
close all

pixelVol = MRId.pixdim(1)*MRId.pixdim(2)*MRId.pixdim(3) % mm^3

[areaE, areaT] = areasFromMasks(mask, thMask);

volE = pixelVol*sum(areaE)
volT = pixelVol*sum(areaT)
volumeViewer(mask)
volumeViewer(thMask)

%%
volumeViewer close
