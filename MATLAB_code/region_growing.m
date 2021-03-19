filename = "volume-covid19-A-0000.nii";
V = niftiread(filename);
whos V
[m,n,z] = size(V);

% c = V(:,:,6); %n is the slice number that you want to visualize.
% imshow(c,[])

figure(1)
volshow(V)
BW = edge3(V,'approxcanny',0.4);
figure(2)
volshow(BW)

% Zone map
ZM = zeros(m,n,z); 
% Counter of the space between two different boundaries
cs = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEEDING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l=1;
i=1;
j=1;
k=1;

% Searching among all the pixel of the image (m X n X z)
while (i<m)
    j=1;
    while(j<n)
        k=1;
        while(k<z)
        % Verify not to exceed the dimension of the image summing cs
            while ((k+cs)<z)
            
            % If the pixel is part of an edge, break
                if (BW(i,j,k) == 1)
                    break
                end
            
            % Enlarge cs until I find another edge, then label the zone
            % with 'l' in the center of the zone (in ZM)
                if ((BW(i,j,k+cs)==1)|| ((k+cs)==(z-1)))
                    u = round(((k+cs)+k)/2);
                    ZM(i,u,k)=l;  
                    l=l+1;
                    break 
                else
                    cs =cs+1;
                end               
            end 
            k = k+cs+1;
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

figure(3)
volshow(ZM)