function neighbors = get_neighbors(point, mode)
    if (mode == "26n")
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
    if (mode == "6n")
            neighbors ={point(1), point(2), point(3)-1;
                        point(1), point(2), point(3)+1;
                        point(1), point(2)-1, point(3);
                        point(1), point(2)+1, point(3);
                        point(1)-1, point(2), point(3);
                        point(1)+1, point(2), point(3)};
    end
end