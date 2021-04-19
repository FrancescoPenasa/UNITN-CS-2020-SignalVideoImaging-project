
function outputMask = rg_function(seeds, intensity_matrix, mask, threshold, mode)
    [m, n, k] = size(intensity_matrix);
    outputMask = zeros(m,n,k);
    
    [seeds_x,seeds_y,seeds_z] = ind2sub(size(seeds),find(seeds == 1));
    
    [seeds_x,seeds_y,seeds_z]
    
    seeds_q = CQueue();
    for i = 1:length(seeds_x)
        seeds_q.push([seeds_x(i), seeds_y(i), seeds_z(i)])
    end
    
    while ~seeds_q.isempty()
        point = seeds_q.pop();
        neighbors = get_neighbors(point, mode);
        
        intensity_point = intensity_matrix(point(1), point(2), point(3)); %intensity seed
        outputMask(point(1), point(2), point(3)) = intensity_point;
        
        for n = 1:length(neighbors)
            
            x=neighbors(n,1);
            y=neighbors(n,2);
            z=neighbors(n,3);
             
            if (x <= 512 && y <= 512 && z <= 341 && x > 0 && y > 0 && z > 0 && mask(x,y,z) == 1)
                
                intensity_neig = intensity_matrix(x,y,z); %intensity neighbors
                                
                if  (intensity_neig <= (intensity_point + threshold) && intensity_neig >= (intensity_point - threshold) && outputMask(x,y,z) == 0)
                    outputMask(x,y,z) = idivide((intensity_neig + intensity_point), 2);
                    seeds_q.push();            
                end
            end
        end
    end
end



function neighbors = get_neighbors(point, mode)
    if mode == "26n"
    neighbors = {point(1)-1, point(2)-1, point(3)-1;   point(1)-1, point(2)-1, point(3);   point(1)-1, point(2)-1, point(3)+1;
                 point(1)-1, point(2), point(3)-1;     point(1)-1, point(2), point(3);     point(1)-1, point(2), point(3)+1;
                 point(1)-1, point(2)+1, point(3)-1;   point(1)-1, point(2)+1, point(3);   point(1)-1, point(2)+1, point(3)+1;
                 point(1), point(2)-1, point(3)-1;     point(1), point(2)-1, point(3);     point(1), point(2)-1, point(3)+1;
                 point(1), point(2), point(3)-1;       point(1), point(2), point(3)+1;     point(1), point(2)+1, point(3)-1;
                 point(1), point(2)+1, point(3);       point(1), point(2)+1, point(3)+1;   point(1)+1, point(2)-1, point(3)-1;
                 point(1)+1, point(2)-1, point(3);     point(1)+1, point(2)-1, point(3)+1; point(1)+1, point(2), point(3)-1;
                 point(1)+1, point(2), point(3);       point(1)+1, point(2), point(3)+1;   point(1)+1, point(2)+1, point(3)-1;
                 point(1)+1, point(2)+1, point(3);     point(1)+1, point(2)+1, point(3)+1};
    end
    if (neighbordMode == "6n")
            neighbors ={point(1), point(2), point(3)-1;
                    point(1), point(2), point(3)+1;
                    point(1), point(2)-1, point(3);
                    point(1), point(2)+1, point(3);
                    point(1)-1, point(2), point(3);
                    point(1)+1, point(2), point(3)};
    end
end
