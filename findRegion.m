function R = findRegion(x_coord, y_coord)
    %% Generate look-up table
    % Image dimensions
    x_size = 340;
    y_size = 180;
    x_blockSize = floor(x_size/3);
    y_blockSize = floor(y_size/3);

    % Construct regions
    grid = [0 1 2; 3 4 5; 6 7 8];
    lut = struct([]);
    for j=1:3
        for i=1:3
            lut = [lut struct(  'x_start', (i-1)*floor(x_blockSize), ...
                                'y_start', (j-1)*floor(y_blockSize), ...
                                'x_index', i, ...
                                'y_index', j ...
                                )];
        end
    end

    %% find region from starting position
    x = x_blockSize * floor(x_coord/x_blockSize);
    y = y_blockSize * floor(y_coord/y_blockSize);
    if x_coord == x_size || x_coord == x_size-1
        x = floor((2/3) * x_size);
    end
    if y_coord == y_size
        y = floor((2/3) * y_size);
    end
    
    x_pos = find([lut(:).x_start] == x);
    y_pos = find([lut(:).y_start] == y);
    index = intersect(x_pos, y_pos);
    
    R = grid(lut(index).y_index, lut(index).x_index);
end