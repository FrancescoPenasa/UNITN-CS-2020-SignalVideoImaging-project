filename = "../Data/volume-covid19-A-0000.nii.gz";
V = niftiread(filename);

V = im2single(V);
% volumeViewer(V);

XY = V(:,:,170);
XZ = squeeze(V(:,256,:));

figure
imshow(XY,[],'Border','tight');

figure
imshow(XZ,[],'Border','tight');

% THRESHOLD SELECTION: 
%  - run imageSegmenter(XY)
%  - click on Threshold
%  - chose manual threshold
%  - move the threshold until you obtain a good segmentation
%  - use that threshold in the next line of code
BW = XY > 4.902000e-01;

BW = imcomplement(BW);
BW = imclearborder(BW);
BW = imfill(BW, 'holes');
radius = 3;
decomposition = 0;
se = strel('disk',radius,decomposition);
BW = imerode(BW, se);
maskedImageXY = XY;
maskedImageXY(~BW) = 0;
imshow(maskedImageXY)

BW = imbinarize(XZ);
BW = imcomplement(BW);
BW = imclearborder(BW);
BW = imfill(BW,'holes');
radius = 13;
decomposition = 0;
se = strel('disk',radius,decomposition);
BW = imerode(BW, se);
maskedImageXZ = XZ;
maskedImageXZ(~BW) = 0;
imshow(maskedImageXZ)

mask = false(size(V));
mask(:,:,170) = maskedImageXY;
mask(:,256,:) = mask(:,256,:)|reshape(maskedImageXZ,[512,1,341]);

V = histeq(V);

BW = activecontour(V,mask,100,'Chan-Vese');

segmentedImage = V.*single(BW);

volumeViewer(segmentedImage);
