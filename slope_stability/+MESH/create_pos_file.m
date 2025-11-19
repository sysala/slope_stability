function create_pos_file(coords, norm_E_node, filename)
    % CREATE_POS_FILE Generates a .pos file from coordinates and scalar values.
    % 
    % Inputs:
    %   coords - 3xN matrix containing x, y, z coordinates of the points.
    %   norm_E_node - 1xN matrix containing scalar values at each point.
    %   filename - Name of the output .pos file.
    %
    % Example usage:
    %   create_pos_file(coords, norm_E_node, 'solution_coarse.pos');
    
    % Transformation of norm_E_node
    % norm_E_node = (1 ./ norm_E_node * 3).^(1/3) * 1.56 / 2.175;
    
    % Extracting the coordinates
    x = coords(1, :);
    y = coords(2, :);
    z = coords(3, :);
    
    % Calculate transformed coordinates
    transformed_x = x;
    
    % Preallocate cell array for lines
    numPoints = size(coords, 2);
    pos_lines = cell(1, numPoints + 2);
    pos_lines{1} = sprintf('View "Field" {\n');
    
    % Generate lines for each point in a loop (avoiding vectorized string issue)
    pos_lines{2} = sprintf('  SP(%f, %f, %f) { %f };\n', [transformed_x; z; y; norm_E_node]);
    pos_lines{3} = sprintf('};\n');
    
    % Write all lines to the file at once
    fid = fopen(sprintf("%s",filename), 'w');
    fprintf(fid, '%s', pos_lines{:});
    fclose(fid);
end
