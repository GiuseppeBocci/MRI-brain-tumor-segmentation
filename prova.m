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

lmasks = size(mask);
lenS = lmasks(3);
pixelVol = MRId.pixdim(1)*MRId.pixdim(2)*MRId.pixdim(3) % mm^3
areaE = zeros(lenS,1);
areaT = zeros(lenS,1);
for s = 1:lenS
    for i = 1:lmasks(1)
        for j = 1:lmasks(2)
            if(mask(i,j,s) ~= 0)
                areaE(s,1) = areaE(s,1) + 1;
            end
            if(thMask(i,j,s) ~= 0)
                areaT(s,1) = areaT(s,1) + 1;
            end
        end
    end
end

volE = pixelVol*sum(areaE)
volT = pixelVol*sum(areaT)
volumeViewer(mask)
volumeViewer(thMask)

%%
volumeViewer close

%%
numSlices = round(MRId.pixdim(3)*size(v,3)/MRId.pixdim(1));
DTransverseIsotropic = imresize3(v,[256 256 numSlices]);
DSagittal = permute(DTransverseIsotropic,[1 3 2]);
DSagittalRotated = imrotate3(DSagittal,90,[0 0 1],"cubic");
figure
montage(DSagittalRotated,colormap("gray"))
title("Sagittal Slices")
volumeViewer(DSagittalRotated)
figure 
imshow(DSagittalRotated(:,:,135))
%%
% figure
% imshow(I)
% h = 2*[0 -1 0; -1 4 -1; 0 -1 0];
% fI = imfilter(I,h);
% figure
% imshow(fI)


%%
close all
I = imcrop(v(:,:,65),rect);
[BWs, threshold] = edge(I,'sobel');
figure
imshow(BWs)
% [BWs, threshold] = edge(I,'log');
% figure
% imshow(BWs)
% [BWs, threshold] = edge(I,'canny');
% figure
% imshow(BWs)


%%
close all
a = load('MRIdata.mat')
D = a

disp( size(D) )       % it can be used by "montage" and "immovie" functions

v = D.vol;
figure
I = v(:,:,73);
imshow(v(:,:,73))
figure
Y = fft2(I,2^nextpow2(a.dim(1)),2^nextpow2(a.dim(2)));
imagesc(abs(fftshift(Y)));

volumeViewer(v)
%% high pass
l = [-1 -1 -1; -1 8 -1; -1 -1 -1];
flI = imfilter(I,l);
figure
imshow(flI)

%% high pass
h = 2*[0 -1 0; -1 4 -1; 0 -1 0];
fI = imfilter(I,h);
figure
imshow(fI)

%%
[centers, radii, metric] = imfindcircles(fI,[2 50]);
centers
viscircles(centers, radii,'EdgeColor','b');

%%
% https://it.mathworks.com/help/images/detecting-a-cell-using-image-segmentation.html
[~,threshold] = edge(v(:,:,73),'sobel');
fudgeFactor = 0.5; % introduced to allow a margin of error in unknown quantities
BWs = edge(v(:,:,73),'sobel',threshold * fudgeFactor);
figure
imshow(BWs)
%%
figure
Y = fft2(BWs,2^nextpow2(a.dim(1)),2^nextpow2(a.dim(2)));
imagesc(abs(fftshift(Y)));
h = 20*[0 -1 0; -1 4 -1; 0 -1 0];
fBWs = imfilter(BWs,h)
figure
imshow(fBWs)

%%
figure
J = imadjust(I);
imshow(J)
figure
histogram(I)
figure
histogram(J)

[~,threshold] = edge(v(:,:,73),'sobel');
fudgeFactor = 1;
BWs = edge(v(:,:,73),'sobel',threshold * fudgeFactor);
figure
imshow(BWs)
%visboundaries(v(:,:,73))
%volumeViewer(v)
%%
se90 = strel('line',3,90);
se0 = strel('line',3,0);
BWsdil = imdilate(BWs,[se90 se0]);
imshow(BWsdil)
title('Dilated Gradient Mask')
%%
BWdfill = imfill(BWsdil,'holes');
imshow(BWdfill)
title('Binary Image with Filled Holes')
%%
BWnobord = imclearborder(BWdfill,4);
imshow(BWnobord)
title('Cleared Border Image')
%%
seD = strel('diamond',1);
BWfinal = imerode(BWnobord,seD);
BWfinal = imerode(BWfinal,seD);
imshow(BWfinal)
title('Segmented Image');
%%
imshow(labeloverlay(I,BWfinal))
title('Mask Over Original Image')

