filename = "../Data/volume-covid19-A-0000.nii.gz";
V = niftiread(filename);

% Convert image to single precision
V = im2single(V);

% Create an intensity matrix for the region growing, starting from the
% original image
intensity_matrix = imadjustn(V); intensity_matrix = im2uint8(intensity_matrix);

% Select two slices in the XY and XZ planes
xy_slice = 170;
xz_slice = 280;

XY = V(:,:,xy_slice);
XZ = squeeze(V(:,xz_slice,:));


%% PLANE XY
% Create black and wite image by selecting a threshold 
BW = XY > 4.902000e-01;

figure; imshow(BW); title("BW")

% Isolate the lungs from the XY image
BW = imcomplement(BW);
BW = imclearborder(BW);
BW = imfill(BW, 'holes');
radius = 5;
decomposition = 0;
se = strel('disk',radius,decomposition);
BW = imerode(BW, se);
figure; imshow(BW);

intensity_XY = intensity_matrix(:,:,xy_slice);
[m,n] = size(BW);
lungs_XY = zeros(m,n); lungs_XY = im2uint8(lungs_XY);

% select the the a piece of the XY image only containing the lungs
for i=1:m
    for j=1:n
        if BW(i,j) == 1
           lungs_XY(i,j) = intensity_XY(i,j);
        end
    end
end

% Extract the edges of the lungs
edges_XY = edge(lungs_XY);
figure; imshow(edges_XY);


%% PLANE XY
% (Same procedure as for the XY plane)
BW = imbinarize(XZ);
BW = imcomplement(BW);
BW = imclearborder(BW);
BW = imfill(BW,'holes');
radius = 5;
decomposition = 0;
se = strel('disk',radius,decomposition);
BW = imerode(BW, se);
figure; imshow(BW);

intensity_XZ = intensity_matrix(:,xz_slice,:);
[m,n] = size(BW);
seeds_XZ = zeros(m,n); seeds_XZ = im2uint8(seeds_XZ);

for i=1:m
    for j=1:n
        if BW(i,j) == 1
           seeds_XZ(i,j) = intensity_XZ(i,j);
        end
    end
end

edges_XZ = edge(seeds_XZ);
figure; imshow(edges_XZ);

%% Create seeds matrix

% Create the final seeds matrix merging the lungs edges in the XY and XZ
% planes
seeds = zeros(size(V)); seeds = im2uint8(seeds);
seeds(:,:,xy_slice) = edges_XY;
seeds(:,xz_slice,:) = seeds(:,xz_slice,:)|reshape(edges_XZ,[512,1,341]);

volumeViewer(seeds);


%% Region Growing

% Create a mask selecting the region of the image to be considerated for
% the region growing (a matrix full of ones means to execute the algorithm on the entire image)
mask = true(size(intensity_matrix));

% Run the Region Growing algorithm
tic;
output = regionGrowing3d(seeds, intensity_matrix, mask, 2, "6n", 7500000);
toc;

volumeViewer(output);


