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
                    w = round(((k+cs)+k)/2);
                    ZM(i,j,w)=l;  
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


% Create a matrix of edge starting from ZM
CM  = zeros (m,n,z);
DN  = im2uint8(V);
DNM = zeros (m,n,z);
DNM = im2uint8(DNM);
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

volshow(CM)
% Removing seeds 2px-away from existing one
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


volshow(CM)

dnt =2;
totc =0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Growing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while totc < 300
 % 4-Connected growing
     for i = 1:m
        for j = 1:n
            for k=1:z
                if (PR(i,j,k) == 0 && CM(i,j,k)==1 && CC(i,j,k)==0)
                    if (j>1)
                        dsx =abs(DNM(i,j,k)-DN(i,j-1,k));
                    end
                    if (j<n)
                        ddx =abs(DNM(i,j,k)-DN(i,j+1,k));
                    end
                    if (i>1)
                        dup =abs(DNM(i,j,k)-DN(i-1,j,k));
                    end
                    if (i<m)
                        ddown=abs(DNM(i,j,k)-DN(i+1,j,k));
                    end
                    if (k>1)
                        dinside=abs(DNM(i,j,k)-DN(i,j,k-1));
                    end
                    if (k<z)
                        doutside=abs(DNM(i,j,k)-DN(i,j,k+1));
                    end
                    % Growing to SX
                    if (j>1   && dsx<dnt && (DDIF(i,j-1,k)>dsx))

                        ZM(i,j-1,k) = ZM (i,j,k);
                        CM(i,j-1,k)= 1;
                        DNM(i,j-1,k) = DNM (i,j,k);
                        CC(i,j-1,k) = 1;
                        DDIF(i,j-1,k) = dsx;


                    end
                     % Growing to DX
                    if (j<n  && (ddx)<dnt && (DDIF(i,j+1,k)>ddx))

                        ZM(i,j+1,k) = ZM (i,j,k);
                        CM(i,j+1,k) =1;
                        DNM(i,j+1,k) = DNM (i,j,k);
                        CC(i,j+1,k) = 1;
                        DDIF(i,j+1,k)= ddx;


                    end

                     % Growing to N
                    if (i>1   && (dup)<dnt &&  (DDIF(i-1,j,k)>dup))

                        ZM(i-1,j,k) = ZM (i,j,k);
                        CM(i-1,j,k) =1;
                        DNM(i-1,j,k) = DNM (i,j,k);
                        CC(i-1,j,k) = 1;
                        DDIF(i-1,j,k) = dup;

                    end

                    % Growing to S
                    if (i<m && (ddown)<dnt && (DDIF(i+1,j,k)>ddown))
                        ZM(i+1,j,k) = ZM (i,j,k);
                        CM(i+1,j,k) =1;
                        DNM(i+1,j,k) = DNM (i,j,k);
                        CC(i+1,j,k) = 1;
                        DDIF(i+1,j,k) = ddown;

                    end
                    
                    % Growing inside
                    if (k>1   && (dinside)<dnt &&  (DDIF(i,j,k-1)>dinside))

                        ZM(i,j,k-1) = ZM (i,j,k);
                        CM(i,j,k-1) =1;
                        DNM(i,j,k-1) = DNM (i,j,k);
                        CC(i,j,k-1) = 1;
                        DDIF(i,j,k-1) = dinside;

                    end

                    
                    % Growing outside
                    if (k<z && (doutside)<dnt && (DDIF(i,j,k+1)>doutside))
                        ZM(i,j,k+1) = ZM (i,j,k);
                        CM(i,j,k+1) =1;
                        DNM(i,j,k+1) = DNM (i,j,k);
                        CC(i,j,k+1) = 1;
                        DDIF(i,j,k+1) = doutside;

                    end

                   CM(i,j,k)=0; 
                   PR (i,j,k) =1;
                end
            end
         end
         CC= zeros(m,n,z);
      totc = totc+1;

     end
end

imshow(DNM(:,:,100));
volshow(DNM);

imshow(CM(:,:,200));

% Not Processed pixels: using 4-connected growing some pixels are left
% behind==> fill them.

totc =0;
while totc<50
    for i=1:m
        for j=1:n
            for k=1:z
                if(ZM(i,j,k)==0)
                    cs=1;
                    cw=1;
                    while(cw<100)

                    if((j-cs)>0 &&(ZM(i,j-cs,k)>0) )
                        ZM(i,j)=ZM(i,j-cs,k);
                        DNM(i,j)=DNM(i,j-cs,k);
                        PR(i,j)=1;
                        break
                    end
                    if((j+cs)<n &&(ZM(i,j+cs,k)>0))
                        ZM(i,j,k)=ZM(i,j+cs,k);
                        DNM(i,j,k)=DNM(i,j+cs,k);
                        PR(i,j,k)=1;
                        break
                    end
                    if( (i-cs)>0 &&(ZM(i-cs,j,k)>0))
                        ZM(i,j,k)=ZM(i-cs,j,k);
                        DNM(i,j,k)=DNM(i-cs,j,k);
                        PR(i,j,k)=1;
                        break
                    end
                    if( (i+cs<m) && (ZM(i+cs,j,k)>0))
                        ZM(i,j,k)=ZM(i+cs,j,k);
                        DNM(i,j,k)=DNM(i+cs,j,k);
                        PR(i,j,k)=1;
                        break
                    end
                    if( (k-cs)>0 &&(ZM(i,j,k-cs)>0))
                        ZM(i,j,k)=ZM(i,j,k-cs);
                        DNM(i,j,k)=DNM(i,j,k-cs);
                        PR(i,j,k)=1;
                        break
                    end
                    if( (k+cs<z) && (ZM(i,j,k+cs)>0))
                        ZM(i,j,k)=ZM(i,j,k+cs);
                        DNM(i,j,k)=DNM(i,j,k+cs);
                        PR(i,j,k)=1;
                        break
                    end
                    cw = cw+1;
                    cs=cs+1;
                    end
                end
            end
        end
       totc= totc+1;
    end
end



volshow(DNM);
imshow(DNM(:,:,160));
