function neighbors = get_neighbors(point, mode)
    if (mode == "26n")
    neighbors = {point(1)-1, point(2)-1, point(3)-1;        point(1)-1, point(2)-1, point(3);       point(1)-1, point(2)-1, point(3)+1;    % x-1 ; y-1 ; z+-1
                 point(1)-1, point(2),   point(3)-1;        point(1)-1, point(2),   point(3);       point(1)-1, point(2),   point(3)+1;    % x-1 ; y   ; z+-1
                 point(1)-1, point(2)+1, point(3)-1;        point(1)-1, point(2)+1, point(3);       point(1)-1, point(2)+1, point(3)+1;    % x-1 ; y+1 ; z+-1
                 
                 point(1),   point(2)-1, point(3)-1;        point(1),   point(2)-1, point(3);       point(1),   point(2)-1, point(3)+1;    % x   ; y-1 ; z+-1
                 point(1),   point(2),   point(3)-1;                                                point(1),   point(2),   point(3)+1;    % x   ; y   ; z+-1 
                 point(1),   point(2)+1, point(3)-1;        point(1),   point(2)+1, point(3);       point(1),   point(2)+1, point(3)+1;    % x   ; y+1 ; z+-1
                 
                 point(1)+1, point(2)-1, point(3)-1;        point(1)+1, point(2)-1, point(3);       point(1)+1, point(2)-1, point(3)+1;    % x+1 ; y+1 ; z+-1
                 point(1)+1, point(2),   point(3)-1;        point(1)+1, point(2),   point(3);       point(1)+1, point(2),   point(3)+1;    % x+1 ; y+1 ; z+-1
                 point(1)+1, point(2)+1, point(3)-1;        point(1)+1, point(2)+1, point(3);       point(1)+1, point(2)+1, point(3)+1};   % x+1 ; y+1 ; z+-1
    end
    if (mode == "6n")
            neighbors ={point(1), point(2), point(3)-1;
                        point(1), point(2), point(3)+1;
                        point(1), point(2)-1, point(3);
                        point(1), point(2)+1, point(3);
                        point(1)-1, point(2), point(3);
                        point(1)+1, point(2), point(3)};
    end
end
