
function outputMask = regionGrowing3d(seeds, intensity_matrix, mask, threshold, mode, iterations)
    with_iter = true;
    if nargin == 5
        with_iter = false;
    end
    
    [m, n, k] = size(intensity_matrix);
    outputMask = zeros(m,n,k);
    outputMask=im2uint8(outputMask);
    
    % Find the coordinates of all the initial seeds
    [seeds_x,seeds_y,seeds_z] = ind2sub(size(seeds),find(seeds == 1));
    
    % Create a Queue data structure for processing the seeds
    seeds_q = CQueue();
    
    % Insert seeds into the QUEUE and assing to the output matrix the
    % intensity of the seeds
    for i = 1:length(seeds_x)
        if (mask(seeds_x(i), seeds_y(i), seeds_z(i)) == 1)
            seeds_q.push([seeds_x(i), seeds_y(i), seeds_z(i)])
            outputMask(seeds_x(i), seeds_y(i), seeds_z(i)) = intensity_matrix(seeds_x(i), seeds_y(i), seeds_z(i));
        end
    end
    
    % Initialize waitbar to look at the execution progression
    if with_iter
        c = 1;
        f = waitbar(0,"Execution of 3D Region Growing algorithm");
    end
    
    % Iterate untill the Queue is empty or the maximum number of iterations
    % is reached
    while ~seeds_q.isempty()
        
        if with_iter
            if c > iterations
                break;
            end
            % Update waitbar
            progression = (c / iterations) * 100;
            if mod(progression,10) == 0
               waitbar(progression/100,f);
            end
        end
        
        % Extract current seed from the Queue and find its neighbors
        point = seeds_q.pop();
        neighbors = get_neighbors(point, mode);

        % Get the intensity of the seed point
        intensity_point =  intensity_matrix(point(1), point(2), point(3)); %intensity seed
        
        % Process all the neighbors
        for neigh = 1:length(neighbors)

            x=cell2mat(neighbors(neigh,1));
            y=cell2mat(neighbors(neigh,2));
            z=cell2mat(neighbors(neigh,3));

            % Check that the neighbors satisfies all the boundary conditions
            if (x <= m && y <= n && z <= k && x > 0 && y > 0 && z > 0 && mask(x,y,z) == 1 && outputMask(x,y,z) == 0)

                % Get the intensity of the neighbor
                intensity_neig = intensity_matrix(x,y,z);

                % If the difference of the intensities between the seed and
                % its neighbor is below the given threshold 
                if  intensity_neig <= (intensity_point + threshold) && intensity_neig >= (intensity_point - threshold)
                    outputMask(x,y,z) = outputMask(point(1), point(2), point(3));
                    seeds_q.push([x, y, z]);
                    if with_iter
                        c = c+1;
                    end
                end
            end
        end
    end
    close(f);
end

