filename = "../Data/volume-covid19-A-0000.nii.gz";
V = niftiread(filename);
whos V
V = im2single(V);
[m,n,z] = size(V);
J = imadjust(V(:,:,160));
figure;
imshow(J); title('Original image')

D=imadjustn(V);
figure;
montage(D,'Indices', 150:169);title('Original image') % multiple plot

% c = V(:,:,6); %n is the slice number that you want to visualize.
% imshow(c,[])

XY = V(:,:,160);
XZ = squeeze(V(256,:,:));
figure
imshow(XY,[],'Border','tight');
imshow(XZ,[],'Border','tight');


% Gaussian filter ?
sigma=2;
volSmooth = imgaussfilt3(D, sigma);
figure;
montage(volSmooth,'Indices', 150:169); title('Gaussian filtered image volume')



figure(1)
volshow(V)
BW = edge3(V,'approxcanny',0.4);

figure(2)
volshow(BW)

% Filter and clear border
D=imadjustn(V);
%Gaussian
sigma=6;
volSmooth = imgaussfilt3(D, sigma);
% Box
volSmooth=imboxfilt3(volSmooth,[5 5 3]);
%figure;
%montage(volSmooth,'Indices', 150:169); title('Gaussian filtered image volume')
%edge
BW = edge3(volSmooth,'approxcanny',0.4);
figure
imshow(BW(:,:,50)); title('Edges');
volumeViewer(BW);
BWc2 = imclearborder(BW,8);
volumeViewer(BWc2)

% remove edges
% cut image XY
XY=zeros(341,251);
for k=1:z
XY = imcrop(BWc2(:,:,z),[160 90 250 340]); %cut figure [xmin ymin width height]
XY = imclearborder(XY,18);
XY = bwpropfilt(XY,'EulerNumber',[0 0]);
XY = imresize(XY,[512 512]);
BW(:,:,z)=XY;
end 


% cut image YZ
YZ=zeros(331,276);
for j=1:n
YZ = imcrop(squeeze(BWc2(:,j,:)),[50 95 275 330]); %cut figure [xmin ymin width height]
YZ = imclearborder(YZ,18);
YZ = bwpropfilt(YZ,'EulerNumber',[0 0]);
YZ= imresize(YZ,[512 341]);
BW(:,j,:)=YZ;
end

% cut image XZ
XZ=zeros(241,261);
for i=1:m
XZ = imcrop(squeeze(BWc2(i,:,:)),[60 160 260 240]); %cut figure [xmin ymin width height]
XZ = imclearborder(XZ,18);
XZ = bwpropfilt(XZ,'EulerNumber',[0 0]);
XZ = imresize(XZ,[512 341]);
BW(i,:,:)=XZ;
end

BW=imclearborder(BW,8);
volumeViewer(BW)
J=medfilt3(BW);
volumeViewer(J)






[m,n,z]=size(V); 

%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEEDING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l=1;
i=1;
j=1;
k=1;
% Zone map
ZM = zeros(m,n,z); 
% Counter of the space between two different boundaries
cs = 1;

% Searching among all the pixel of the image (m X n X z) 
 
while (i<m)
    j=1;
    while(j<n)
        k=1;
        while(k<z)
        % Verify not to exceed the dimension of the image summing cs
            while ((k+cs)<z-1)
            % If the pixel is part of an edge, break
                if (BW(i,j,k) == 1)
                    break
                end
            % Enlarge cs until I find another edge, then label the zone
            % with 'l' in the center of the zone (in ZM)
                if ((BW(i,j,k+cs)==1)|| ((k+cs)==(z-1)))
                    w = round(((k+cs)+k)/2);
                    ZM(i,j,w)=l;  
                    l=l+1;
                    break 
                else
                    cs =cs+1;
                end 
                k = k+cs+1;
            end 
            k = k+1;
            cs =1;
        end
        % Increment the value for searching (manually, it's not a for
        % cycle)
        j = j+1;
    end
    % Increment the value for searching (manually, it's not a for
    % cycle)
    i = i+1;
end

% volumeViewer(ZM)
% Create a matrix of edge starting from ZM
CM  = zeros (m,n,z);
DN  = im2uint8(V); %figure
DNM = zeros (m,n,z);
DNM = im2uint8(DNM); %immagini con puntini 
DDIF= 255*ones(m,n,z);

% Matrix to check if a seed has been processed
PR  = zeros (m,n,z);

% Matrix that checks if a pixel has been processed in a certain cycle
% (infer uniform enlargment in all directions)
CC  = zeros (m,n,z);

% Matrix that checks the decimation process
CD = zeros(m,n,z);

%%%%%%%%%%%%%%%%%%%%%%%
% REMOVING SEED IN EXCESS
%%%%%%%%%%%%%%%%%%%%%%%%

% Removing the seed in vertical direction
for i=1:m
    for j=1:n
        for k=1:z
            if (ZM(i,j,k)>0 && CD(i,j,k)==0)
                c=1;
                di=ZM(i,j,k);
                while (((i+c)<m) && (ZM(i+c,j,k)>0))
                    c=c+1;
                    ZM(i+c-1,j,k) = 0;
                end
                if (c>1)
                    u = round(i+(c/2));
                    ZM(u,j,k)= di;
                    CM(u,j,k)=1;
                    CD(u,j,k)=1;
                    DNM (u,j,k)= DN (i,j,k);
                    ZM(i,j,k) =0;
                else
                 DNM(i,j,k)=DN(i,j,k);
                 CM(i,j,k)=1;
                end
            end
        end
    end
end

CD = zeros(m,n,z);


% Examinating the NE-SW direction
for i=1:m
    for j=1:n
        for k=1:z
            if (ZM(i,j,k)>0 && CD(i,j,k)==0)
                c=1;
                di=ZM(i,j,k); 
            while ((i+c<m) && (j-c>1) && (k-c>1) && ZM(i+c,j-c,k-c)>0)
                    c=c+1;
                    ZM(i+c-1,j-c+1,k-c+1) = 0;
                    CM(i+c-1,j-c+1,k-c+1) = 0;
            end
                if (c>1)
                    u = round(i+(c/2));
                    v = round(j-(c/2));
                    w = round(k-(c/2));
                    ZM(u,v,w)= di;
                    CM(u,v,w)=1;
                    CD(u,v,w)=1;
                    DNM (u,v,w)= DN (i,j,k);
                    ZM(i,j,k) =0;
                    CM(i,j,k)=0;
                end
            end
        end
    end
end

CD = zeros(m,n,z);

% Examinating the NW-SE direction
for i=1:m
    for j=1:n
        for k=1:z
            if (ZM(i,j,k)>0 && CD(i,j,k)==0)
                c=1;
                di=ZM(i,j,k); 
                while ((i+c<m) && (j+c<n) && (k+c<z) && ZM(i+c,j+c,k+c)>0)
                    c=c+1;
                    ZM(i+c-1,j+c-1,k+c-1) = 0;
                    CM(i+c-1,j+c-1,k+c-1) = 0;
                end
                if (c>1)
                    u = round(i+(c/2));
                    v = round(j+(c/2));
                    w = round(k+(c/2));
                    ZM(u,v,w)= di;
                    CM(u,v,w)=1;
                    CD(u,v,w)=1;
                    DNM (u,v,w)= DN (i,j,k);
                    ZM(i,j,k) =0;
                    CM(i,j,k)=0;
                end
            end
        end
    end
end

% figure
% volumeViewer(CM)
% figure;volumeViewer(DN)

%Removing seeds 2px-away from existing one
for i=1:m
    for j=1:n
        for k=1:z
            if (CM(i,j,k)==1 && i<(m-2) && j<(n-2) && k<(z-2) && i>2 && j>2 && k>2)
                CM(i-2,j-2,k-2)=0;
                CM(i-2,j-1,k-2)=0;
                CM(i-2,j,k-2)=0;
                CM(i-2,j+1,k-2)=0;
                CM(i-2,j+2,k-2)=0;
                CM(i-1,j-2,k-2)=0;
                CM(i-1,j-1,k-2)=0;
                CM(i-1,j,k-2)=0;
                CM(i-1,j+1,k-2)=0;
                CM(i-1,j+2,k-2)=0;
                CM(i,j-2,k-2)=0;
                CM(i,j-1,k-2)=0;
                CM(i,j+1,k-2)=0;
                CM(i,j+2,k-2)=0;
                CM(i+1,j-2,k-2)=0;
                CM(i+1,j-1,k-2)=0;
                CM(i+1,j,k-2)=0;
                CM(i+1,j+1,k-2)=0;
                CM(i+1,j+2,k-2)=0;
                CM(i+2,j-2,k-2)=0;
                CM(i+2,j-1,k-2)=0;
                CM(i+2,j,k-2)=0;
                CM(i+2,j+1,k-2)=0;
                CM(i+2,j+2,k-2)=0;
                
                CM(i-2,j-2,k-1)=0;
                CM(i-2,j-1,k-1)=0;
                CM(i-2,j,k-1)=0;
                CM(i-2,j+1,k-1)=0;
                CM(i-2,j+2,k-1)=0;
                CM(i-1,j-2,k-1)=0;
                CM(i-1,j-1,k-1)=0;
                CM(i-1,j,k-1)=0;
                CM(i-1,j+1,k-1)=0;
                CM(i-1,j+2,k-1)=0;
                CM(i,j-2,k-1)=0;
                CM(i,j-1,k-1)=0;
                CM(i,j+1,k-1)=0;
                CM(i,j+2,k-1)=0;
                CM(i+1,j-2,k-1)=0;
                CM(i+1,j-1,k-1)=0;
                CM(i+1,j,k-1)=0;
                CM(i+1,j+1,k-1)=0;
                CM(i+1,j+2,k-1)=0;
                CM(i+2,j-2,k-1)=0;
                CM(i+2,j-1,k-1)=0;
                CM(i+2,j,k-1)=0;
                CM(i+2,j+1,k-1)=0;
                CM(i+2,j+2,k-1)=0;
                
                CM(i-2,j-2,k)=0;
                CM(i-2,j-1,k)=0;
                CM(i-2,j,k)=0;
                CM(i-2,j+1,k)=0;
                CM(i-2,j+2,k)=0;
                CM(i-1,j-2,k)=0;
                CM(i-1,j-1,k)=0;
                CM(i-1,j,k)=0;
                CM(i-1,j+1,k)=0;
                CM(i-1,j+2,k)=0;
                CM(i,j-2,k)=0;
                CM(i,j-1,k)=0;
                CM(i,j+1,k)=0;
                CM(i,j+2,k)=0;
                CM(i+1,j-2,k)=0;
                CM(i+1,j-1,k)=0;
                CM(i+1,j,k)=0;
                CM(i+1,j+1,k)=0;
                CM(i+1,j+2,k)=0;
                CM(i+2,j-2,k)=0;
                CM(i+2,j-1,k)=0;
                CM(i+2,j,k)=0;
                CM(i+2,j+1,k)=0;
                CM(i+2,j+2,k)=0;
                
                CM(i-2,j-2,k+1)=0;
                CM(i-2,j-1,k+1)=0;
                CM(i-2,j,k+1)=0;
                CM(i-2,j+1,k+1)=0;
                CM(i-2,j+2,k+1)=0;
                CM(i-1,j-2,k+1)=0;
                CM(i-1,j-1,k+1)=0;
                CM(i-1,j,k+1)=0;
                CM(i-1,j+1,k+1)=0;
                CM(i-1,j+2,k+1)=0;
                CM(i,j-2,k+1)=0;
                CM(i,j-1,k+1)=0;
                CM(i,j+1,k+1)=0;
                CM(i,j+2,k+1)=0;
                CM(i+1,j-2,k+1)=0;
                CM(i+1,j-1,k+1)=0;
                CM(i+1,j,k+1)=0;
                CM(i+1,j+1,k+1)=0;
                CM(i+1,j+2,k+1)=0;
                CM(i+2,j-2,k+1)=0;
                CM(i+2,j-1,k+1)=0;
                CM(i+2,j,k+1)=0;
                CM(i+2,j+1,k+1)=0;
                CM(i+2,j+2,k+1)=0;
                
                CM(i-2,j-2,k+2)=0;
                CM(i-2,j-1,k+2)=0;
                CM(i-2,j,k+2)=0;
                CM(i-2,j+1,k+2)=0;
                CM(i-2,j+2,k+2)=0;
                CM(i-1,j-2,k+2)=0;
                CM(i-1,j-1,k+2)=0;
                CM(i-1,j,k+2)=0;
                CM(i-1,j+1,k+2)=0;
                CM(i-1,j+2,k+2)=0;
                CM(i,j-2,k+2)=0;
                CM(i,j-1,k+2)=0;
                CM(i,j+1,k+2)=0;
                CM(i,j+2,k+2)=0;
                CM(i+1,j-2,k+2)=0;
                CM(i+1,j-1,k+2)=0;
                CM(i+1,j,k+2)=0;
                CM(i+1,j+1,k+2)=0;
                CM(i+1,j+2,k+2)=0;
                CM(i+2,j-2,k+2)=0;
                CM(i+2,j-1,k+2)=0;
                CM(i+2,j,k+2)=0;
                CM(i+2,j+1,k+2)=0;
                CM(i+2,j+2,k+2)=0;
                
                
                
                
                

            end
        end
    end
end
%Removing seeds c-px-away from existing one
c=1;
b=1;
for i=1:m
    for j=1:n
        for k=1:z
            for c=1:10
                for b=1:10
                    l=max(c,b);
                    if (CM(i,j,k)==1 && i<(m-l) && j<(n-l) && k<(z-l) && i>l && j>l && k>l)
                        CM(i-c,j-c,k-c)=0;
                        CM(i-c,j-b,k-c)=0;
                        CM(i-c,j,k-c)=0;
                        CM(i-c,j+b,k-c)=0;
                        CM(i-c,j+c,k-c)=0;
                        CM(i-b,j-c,k-c)=0;
                        CM(i-b,j-b,k-c)=0;
                        CM(i-b,j,k-c)=0;
                        CM(i-b,j+b,k-c)=0;
                        CM(i-b,j+c,k-c)=0;
                        CM(i,j-c,k-c)=0;
                        CM(i,j-b,k-c)=0;
                        CM(i,j+b,k-c)=0;
                        CM(i,j+c,k-c)=0;
                        CM(i+b,j-c,k-c)=0;
                        CM(i+b,j-b,k-c)=0;
                        CM(i+b,j,k-c)=0;
                        CM(i+b,j+b,k-c)=0;
                        CM(i+b,j+c,k-c)=0;
                        CM(i+c,j-c,k-c)=0;
                        CM(i+c,j-b,k-c)=0;
                        CM(i+c,j,k-c)=0;
                        CM(i+c,j+b,k-c)=0;
                        CM(i+c,j+c,k-c)=0;

                        CM(i-c,j-c,k-b)=0;
                        CM(i-c,j-b,k-b)=0;
                        CM(i-c,j,k-b)=0;
                        CM(i-c,j+b,k-b)=0;
                        CM(i-c,j+c,k-b)=0;
                        CM(i-b,j-c,k-b)=0;
                        CM(i-b,j-b,k-b)=0;
                        CM(i-b,j,k-b)=0;
                        CM(i-b,j+b,k-b)=0;
                        CM(i-b,j+c,k-b)=0;
                        CM(i,j-c,k-b)=0;
                        CM(i,j-b,k-b)=0;
                        CM(i,j+b,k-b)=0;
                        CM(i,j+c,k-b)=0;
                        CM(i+b,j-c,k-b)=0;
                        CM(i+b,j-b,k-b)=0;
                        CM(i+b,j,k-b)=0;
                        CM(i+b,j+b,k-b)=0;
                        CM(i+b,j+c,k-b)=0;
                        CM(i+b,j-c,k-b)=0;
                        CM(i+c,j-b,k-b)=0;
                        CM(i+c,j,k-b)=0;
                        CM(i+c,j+b,k-b)=0;
                        CM(i+c,j+c,k-b)=0;

                        CM(i-c,j-c,k)=0;
                        CM(i-c,j-b,k)=0;
                        CM(i-c,j,k)=0;
                        CM(i-c,j+b,k)=0;
                        CM(i-c,j+c,k)=0;
                        CM(i-b,j-c,k)=0;
                        CM(i-b,j-b,k)=0;
                        CM(i-b,j,k)=0;
                        CM(i-b,j+b,k)=0;
                        CM(i-b,j+c,k)=0;
                        CM(i,j-c,k)=0;
                        CM(i,j-b,k)=0;
                        CM(i,j+b,k)=0;
                        CM(i,j+c,k)=0;
                        CM(i+b,j-c,k)=0;
                        CM(i+b,j-c,k)=0;
                        CM(i+b,j,k)=0;
                        CM(i+b,j+b,k)=0;
                        CM(i+b,j+c,k)=0;
                        CM(i+c,j-c,k)=0;
                        CM(i+c,j-b,k)=0;
                        CM(i+c,j,k)=0;
                        CM(i+c,j+b,k)=0;
                        CM(i+c,j+c,k)=0;

                        CM(i-c,j-c,k+b)=0;
                        CM(i-c,j-b,k+b)=0;
                        CM(i-c,j,k+b)=0;
                        CM(i-c,j+b,k+b)=0;
                        CM(i-c,j+c,k+b)=0;
                        CM(i-b,j-c,k+b)=0;
                        CM(i-b,j-b,k+b)=0;
                        CM(i-b,j,k+b)=0;
                        CM(i-b,j+b,k+b)=0;
                        CM(i-b,j+c,k+b)=0;
                        CM(i,j-c,k+b)=0;
                        CM(i,j-b,k+b)=0;
                        CM(i,j+b,k+b)=0;
                        CM(i,j+c,k+b)=0;
                        CM(i+b,j-c,k+b)=0;
                        CM(i+b,j-b,k+b)=0;
                        CM(i+b,j,k+b)=0;
                        CM(i+b,j+b,k+b)=0;
                        CM(i+b,j+c,k+b)=0;
                        CM(i+c,j-c,k+b)=0;
                        CM(i+c,j-b,k+b)=0;
                        CM(i+c,j,k+b)=0;
                        CM(i+c,j+b,k+b)=0;
                        CM(i+c,j+c,k+b)=0;

                        CM(i-c,j-c,k+c)=0;
                        CM(i-c,j-b,k+c)=0;
                        CM(i-c,j,k+c)=0;
                        CM(i-c,j+b,k+c)=0;
                        CM(i-c,j+c,k+c)=0;
                        CM(i-b,j-c,k+c)=0;
                        CM(i-b,j-b,k+c)=0;
                        CM(i-b,j,k+c)=0;
                        CM(i-b,j+b,k+c)=0;
                        CM(i-b,j+c,k+c)=0;
                        CM(i,j-c,k+c)=0;
                        CM(i,j-c,k+c)=0;
                        CM(i,j+b,k+c)=0;
                        CM(i,j+c,k+c)=0;
                        CM(i+b,j-c,k+c)=0;
                        CM(i+b,j-b,k+c)=0;
                        CM(i+b,j,k+c)=0;
                        CM(i+b,j+b,k+c)=0;
                        CM(i+b,j+c,k+c)=0;
                        CM(i+c,j-c,k+c)=0;
                        CM(i+c,j-b,k+c)=0;
                        CM(i+c,j,k+c)=0;
                        CM(i+c,j+b,k+c)=0;
                        CM(i+c,j+c,k+c)=0;
                    end
                end 
            end
        end
    end
end


figure(4)
volshow(CM)

%%

DN_flat = DN(:);
[GC,GR] = groupcounts(DN_flat);

[r,c,v] = ind2sub(size(CM),find(CM == 1));
sum(CM(:) == 1) %check the number of seeds


len=length(r);
masks=ones(512,512,341);

outputMask=zeros(512,512,341);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Growing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
neighbordMode="26n";
threshold=1;
i=1;
totc=0;
while totc<30
    i=1;
    disp(totc);
    len=length(r);
    disp(len);
    while len>0 %not 
        newItem=[r(i),c(i),v(i)];
        len=len-1;
        %disp(len);
        if neighbordMode == "26n"
             neighbors = {newItem(1)-1, newItem(2)-1, newItem(3)-1;   newItem(1)-1, newItem(2)-1, newItem(3);   newItem(1)-1, newItem(2)-1, newItem(3)+1;
                                 newItem(1)-1, newItem(2), newItem(3)-1;     newItem(1)-1, newItem(2), newItem(3);     newItem(1)-1, newItem(2), newItem(3)+1;
                                 newItem(1)-1, newItem(2)+1, newItem(3)-1;   newItem(1)-1, newItem(2)+1, newItem(3);   newItem(1)-1, newItem(2)+1, newItem(3)+1;
                                 newItem(1), newItem(2)-1, newItem(3)-1;     newItem(1), newItem(2)-1, newItem(3);     newItem(1), newItem(2)-1, newItem(3)+1;
                                 newItem(1), newItem(2), newItem(3)-1;       newItem(1), newItem(2), newItem(3)+1;     newItem(1), newItem(2)+1, newItem(3)-1;
                                 newItem(1), newItem(2)+1, newItem(3);       newItem(1), newItem(2)+1, newItem(3)+1;   newItem(1)+1, newItem(2)-1, newItem(3)-1;
                                 newItem(1)+1, newItem(2)-1, newItem(3);     newItem(1)+1, newItem(2)-1, newItem(3)+1; newItem(1)+1, newItem(2), newItem(3)-1;
                                 newItem(1)+1, newItem(2), newItem(3);       newItem(1)+1, newItem(2), newItem(3)+1;   newItem(1)+1, newItem(2)+1, newItem(3)-1;
                                 newItem(1)+1, newItem(2)+1, newItem(3);     newItem(1)+1, newItem(2)+1, newItem(3)+1};

             for n=1:length(neighbors)
                 x=neighbors(n,1); x=x{:};
                 y=neighbors(n,2); y=y{:};
                 z=neighbors(n,3); z=z{:};
                if (x <= 512 && y <= 512 && z <= 341 && x > 0 && y > 0 && z > 0 && CM(r(i),c(i),v(i)) == 1)
                    intensity_neig = DN(x,y,z); %intensity neighbors
                    intensity_point = DN(r(i),c(i),v(i)); %intensity seed
                    outputMask(r(i),c(i),v(i)) = intensity_point;
                    if  (intensity_neig <= intensity_point+threshold && intensity_neig > intensity_point-threshold && outputMask(x,y,z) == 0)
                        outputMask(x,y,z) = 1;
                        CM(x,y,z) = 1;            
                        outputMask(r(i),c(i),v(i)) = intensity_point;
                        %DNM(x,y,z) = DNM(r(i),c(i),v(i)); %assign at the neighbour the same intensity as the seeds
    %                     r = [r;x]; % add element at the end of the r array
    %                     c = [c;y];
    %                     v = [v;z];
    %                     len=len+1;
                        disp(len);
                    end
                end
             end
        end
        if (neighbordMode == "6n")
            neighbors ={newItem(1), newItem(2), newItem(3)-1;
                    newItem(1), newItem(2), newItem(3)+1;
                    newItem(1), newItem(2)-1, newItem(3);
                    newItem(1), newItem(2)+1, newItem(3);
                    newItem(1)-1, newItem(2), newItem(3);
                    newItem(1)+1, newItem(2), newItem(3)};
             for n=1:length(neighbors)
                 x=neighbors(n,1); x=x{:};
                 y=neighbors(n,2); y=y{:};
                 z=neighbors(n,3); z=z{:};
                if (x <= 512 && y <= 512 && z <= 341 && x > 0 && y > 0 && z > 0 && masks(x,y,z) == 1)
                    intensity = DN(x,y,z);
                    if  (intensity <= intensity+threshold && intensity > intensity-threshold && outputMask(x,y,z) == 0)
                        outputMask(z,y,x) = 1;
    %                     r = [r;x]; % add element at the end of the r array
    %                     c = [c;y];
    %                     v = [v;z];
    %                     len=len+1;
                    end
                end
             end

        end 
        i=i+1;
    end
  totc=totc+1;  
end     
volumeViewer(outputMatrix);



% figure(8)
% volshow(DNM);
% 
% figure(9)
% imshow(imadjust(DNM(:,:,160)));

%% REGION GROWING (Prova) 

masks=ones(512,512,341);
output = rg_function(CM, DN, masks, 1, "6n");
