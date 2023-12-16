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

pixelVol = MRId.pixdim(1)*MRId.pixdim(2)*MRId.pixdim(3); % mm^3
garylevelsbits = 8;

%%
%Gaussian noise
v = imnoise(MRId.vol, 'gaussian', 0, 1e-3);
%v = imnoise(MRId.vol, 'gaussian', 0, 0.1);
%v = imnoise(MRId.vol, 'gaussian', 0.5, 1e-3); 
% v = imnoise(MRId.vol, 'gaussian', -0.5, 1e-3); 
% v = 255-v;

%Salt and pepper
v = imnoise(v, 'salt & pepper', 0.05);
%v = imnoise(MRId.vol, 'salt & pepper', 0.2);
% SP3 = imnoise(v_ax, 'salt & pepper', 0.35);

%% ROI T
close all
volumeViewer close

figure
montage(v)
title("MRI transversal")

%I = v(:,:,slice);
%[cropped, rect] = imcrop(I);
rect = [133.5100  103.5100   51.9800   44.9800];

[mask, thMask] = tumorMasks(v, rect, 65, 90, garylevelsbits);

%% VOLUMES T
close all

[areaE, areaT] = areasFromMasks(mask, thMask);

volE = pixelVol*sum(areaE)
volT = pixelVol*sum(areaT)
volumeViewer(mask)
volumeViewer(thMask)

%% SAGITTAL
% No computational rilevant
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

[mask, thMask] = tumorMasks(sagittal, rect, 113, 150, garylevelsbits); % threshold);

[areaE, areaT] = areasFromMasks(mask, thMask);

volE = pixelVol*sum(areaE)
volT = pixelVol*sum(areaT)
%volumeViewer(mask)
%volumeViewer(thMask)

%% CORONAL
% No computational rilevant
trasversalRot = flip(v, 3); % order of sicels
% from front to back and not the inverse. No computational relevant
trasversalRot = flip(trasversalRot, 2);
% sx and dx not inverted. No computational relevant
trasversalRot = flip(trasversalRot, 1);
coronal = permute(trasversalRot, [3 1 2]);
figure
montage(coronal)
title("MRI coronal")
% slices 81 - 116(115)

%% 
close all
volumeViewer close

%I = coronal(:,:,100);
%[cropped, rect] = imcrop(I);
rect = [109.5100   21.5100   42.9800   29.9800];

[mask, thMask] = tumorMasks(coronal, rect, 81, 116, garylevelsbits); % threshold);

[areaE, areaT] = areasFromMasks(mask, thMask);

volE = pixelVol*sum(areaE)
volT = pixelVol*sum(areaT)
%volumeViewer(mask)
%volumeViewer(thMask)