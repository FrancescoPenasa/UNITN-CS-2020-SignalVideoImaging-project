clear

V = niftiread('../Data/volume-covid19-A-0000.nii.gz');
V = imadjustn(V); V = im2uint8(V);
V = V(91:420,161:400,71:310);

volumeViewer(V);

[m,n,z]=size(V); 

%% Filtering and edge detection

%Gaussian
sigma=6;
volSmooth = imgaussfilt3(V, sigma);

% Box filter
volSmooth=imboxfilt3(volSmooth,[5 5 3]);

%edge
BW = edge3(volSmooth,'approxcanny',0.4);
BW = imclearborder(BW,6);

% Display edges
volumeViewer(BW)

%% SEEDING %%%%%

i=1;
j=1;
k=1;

% Counter of the space between two different boundaries
cs = 1;

seeds = zeros(m,n,z); 

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
                % Enlarge cs until I find another edge, then create a seed
                % in the middle
                if ((BW(i,j,k+cs)==1) || ((k+cs)==(z-1)))
                    w = round(((k+cs)+k)/2);
                    seeds(i,j,w)=1;  
                    break 
                else
                    cs =cs+1;
                end 
                k = k+cs+1;
            end 
            k = k+1;
            cs =1;
        end
        j = j+1;
    end
    i = i+1;
end

% Matrix that checks the decimation process
visited = zeros(m,n,z);

%% REMOVING SEED IN EXCESS %%%%%

% Removing the seed in vertical direction
for i=1:m
    for j=1:n
        for k=1:z
            if (seeds(i,j,k) == 1 && visited(i,j,k)==0)
                c=1;
                while (((i+c)<m) && (seeds(i+c,j,k) == 1))
                    c=c+1;
                    seeds(i+c-1,j,k) = 0;
                end
                if (c>1)
                    u = round(i+(c/2));
                    seeds(u,j,k)= 1;
                    seeds(i,j,k) =0;
                    visited(u,j,k)=1;
                end
            end
        end
    end
end

% Reset visited matrix
visited = zeros(m,n,z);

% Examinating the NE-SW direction
for i=1:m
    for j=1:n
        for k=1:z
            if (seeds(i,j,k) == 1 && visited(i,j,k)==0)
                c=1;
                while ((i+c<m) && (j-c>1) && (k-c>1) && seeds(i+c,j-c,k-c)==1)
                    seeds(i+c-1,j-c+1,k-c+1) = 0;
                    c=c+1;
                end
                if (c>1)
                    u = round(i+(c/2));
                    v = round(j-(c/2));
                    w = round(k-(c/2));

                    seeds(u,v,w) = 1;
                    seeds(i,j,k) = 0;
                    visited(u,v,w) = 1;
                end
            end
        end
    end
end

% Reset visited matrix
visited = zeros(m,n,z);

% Examinating the NW-SE direction
for i=1:m
    for j=1:n
        for k=1:z
            if (seeds(i,j,k) == 1 && visited(i,j,k)==0)
                c=1;
                while ((i+c<m) && (j+c<n) && (k+c<z) && seeds(i+c,j+c,k+c) == 1)
                    seeds(i+c-1,j+c-1,k+c-1) = 0;
                    c=c+1;
                end
                if (c>1)
                    u = round(i+(c/2));
                    v = round(j+(c/2));
                    w = round(k+(c/2));
                    
                    seeds(u,v,w) = 1;
                    seeds(i,j,k) = 0;
                    visited(u,v,w) = 1;
                end
            end
        end
    end
end


%Removing seeds 10-px-away from existing one
n_px = 10;
for i=1:m
    for j=1:n
        for k=1:z
            for c=1:n_px
                for b=1:n_px
                    l=max(c,b);
                    if (seeds(i,j,k)==1 && i<(m-l) && j<(n-l) && k<(z-l) && i>l && j>l && k>l)
                        seeds(i-c,j-c,k-c)=0;
                        seeds(i-c,j-b,k-c)=0;
                        seeds(i-c,j,k-c)=0;
                        seeds(i-c,j+b,k-c)=0;
                        seeds(i-c,j+c,k-c)=0;
                        seeds(i-b,j-c,k-c)=0;
                        seeds(i-b,j-b,k-c)=0;
                        seeds(i-b,j,k-c)=0;
                        seeds(i-b,j+b,k-c)=0;
                        seeds(i-b,j+c,k-c)=0;
                        seeds(i,j-c,k-c)=0;
                        seeds(i,j-b,k-c)=0;
                        seeds(i,j+b,k-c)=0;
                        seeds(i,j+c,k-c)=0;
                        seeds(i+b,j-c,k-c)=0;
                        seeds(i+b,j-b,k-c)=0;
                        seeds(i+b,j,k-c)=0;
                        seeds(i+b,j+b,k-c)=0;
                        seeds(i+b,j+c,k-c)=0;
                        seeds(i+c,j-c,k-c)=0;
                        seeds(i+c,j-b,k-c)=0;
                        seeds(i+c,j,k-c)=0;
                        seeds(i+c,j+b,k-c)=0;
                        seeds(i+c,j+c,k-c)=0;

                        seeds(i-c,j-c,k-b)=0;
                        seeds(i-c,j-b,k-b)=0;
                        seeds(i-c,j,k-b)=0;
                        seeds(i-c,j+b,k-b)=0;
                        seeds(i-c,j+c,k-b)=0;
                        seeds(i-b,j-c,k-b)=0;
                        seeds(i-b,j-b,k-b)=0;
                        seeds(i-b,j,k-b)=0;
                        seeds(i-b,j+b,k-b)=0;
                        seeds(i-b,j+c,k-b)=0;
                        seeds(i,j-c,k-b)=0;
                        seeds(i,j-b,k-b)=0;
                        seeds(i,j+b,k-b)=0;
                        seeds(i,j+c,k-b)=0;
                        seeds(i+b,j-c,k-b)=0;
                        seeds(i+b,j-b,k-b)=0;
                        seeds(i+b,j,k-b)=0;
                        seeds(i+b,j+b,k-b)=0;
                        seeds(i+b,j+c,k-b)=0;
                        seeds(i+b,j-c,k-b)=0;
                        seeds(i+c,j-b,k-b)=0;
                        seeds(i+c,j,k-b)=0;
                        seeds(i+c,j+b,k-b)=0;
                        seeds(i+c,j+c,k-b)=0;

                        seeds(i-c,j-c,k)=0;
                        seeds(i-c,j-b,k)=0;
                        seeds(i-c,j,k)=0;
                        seeds(i-c,j+b,k)=0;
                        seeds(i-c,j+c,k)=0;
                        seeds(i-b,j-c,k)=0;
                        seeds(i-b,j-b,k)=0;
                        seeds(i-b,j,k)=0;
                        seeds(i-b,j+b,k)=0;
                        seeds(i-b,j+c,k)=0;
                        seeds(i,j-c,k)=0;
                        seeds(i,j-b,k)=0;
                        seeds(i,j+b,k)=0;
                        seeds(i,j+c,k)=0;
                        seeds(i+b,j-c,k)=0;
                        seeds(i+b,j-c,k)=0;
                        seeds(i+b,j,k)=0;
                        seeds(i+b,j+b,k)=0;
                        seeds(i+b,j+c,k)=0;
                        seeds(i+c,j-c,k)=0;
                        seeds(i+c,j-b,k)=0;
                        seeds(i+c,j,k)=0;
                        seeds(i+c,j+b,k)=0;
                        seeds(i+c,j+c,k)=0;

                        seeds(i-c,j-c,k+b)=0;
                        seeds(i-c,j-b,k+b)=0;
                        seeds(i-c,j,k+b)=0;
                        seeds(i-c,j+b,k+b)=0;
                        seeds(i-c,j+c,k+b)=0;
                        seeds(i-b,j-c,k+b)=0;
                        seeds(i-b,j-b,k+b)=0;
                        seeds(i-b,j,k+b)=0;
                        seeds(i-b,j+b,k+b)=0;
                        seeds(i-b,j+c,k+b)=0;
                        seeds(i,j-c,k+b)=0;
                        seeds(i,j-b,k+b)=0;
                        seeds(i,j+b,k+b)=0;
                        seeds(i,j+c,k+b)=0;
                        seeds(i+b,j-c,k+b)=0;
                        seeds(i+b,j-b,k+b)=0;
                        seeds(i+b,j,k+b)=0;
                        seeds(i+b,j+b,k+b)=0;
                        seeds(i+b,j+c,k+b)=0;
                        seeds(i+c,j-c,k+b)=0;
                        seeds(i+c,j-b,k+b)=0;
                        seeds(i+c,j,k+b)=0;
                        seeds(i+c,j+b,k+b)=0;
                        seeds(i+c,j+c,k+b)=0;

                        seeds(i-c,j-c,k+c)=0;
                        seeds(i-c,j-b,k+c)=0;
                        seeds(i-c,j,k+c)=0;
                        seeds(i-c,j+b,k+c)=0;
                        seeds(i-c,j+c,k+c)=0;
                        seeds(i-b,j-c,k+c)=0;
                        seeds(i-b,j-b,k+c)=0;
                        seeds(i-b,j,k+c)=0;
                        seeds(i-b,j+b,k+c)=0;
                        seeds(i-b,j+c,k+c)=0;
                        seeds(i,j-c,k+c)=0;
                        seeds(i,j-c,k+c)=0;
                        seeds(i,j+b,k+c)=0;
                        seeds(i,j+c,k+c)=0;
                        seeds(i+b,j-c,k+c)=0;
                        seeds(i+b,j-b,k+c)=0;
                        seeds(i+b,j,k+c)=0;
                        seeds(i+b,j+b,k+c)=0;
                        seeds(i+b,j+c,k+c)=0;
                        seeds(i+c,j-c,k+c)=0;
                        seeds(i+c,j-b,k+c)=0;
                        seeds(i+c,j,k+c)=0;
                        seeds(i+c,j+b,k+c)=0;
                        seeds(i+c,j+c,k+c)=0;
                    end
                end 
            end
        end
    end
end

% Display seeds
volumeViewer(seeds);

%% Region Growing

% Create intensity matrix from V
intensity_matrix  = im2uint8(V);

% Create a mask selecting the region of the image to be considerated for
% the region growing (a matrix full of ones means to execute the algorithm on the entire image)
mask=ones(m,n,z);

% Start Region Growing
format shortg
clock
output = regionGrowing3d(seeds, intensity_matrix, mask, 3, "6n");
format shortg
clock

% Display result
volumeViewer(output)

