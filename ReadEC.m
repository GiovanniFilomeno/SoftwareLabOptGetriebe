function [EC,EC_length,EC_width,EC_height,EC_center, EC_SurfCenter] = ReadEC(file)

    if ~exist(file,'file')
        error(['File ''%s'' not found. If the file is not on MATLAB''s path' ...
               ', be sure to specify the full path to the file.'], file);
    end

    Cube=stlread(file);
    V = Cube.Points;
    F = Cube.ConnectivityList;
    N = faceNormal(Cube);
    EC = struct('faces',F,'vertices',V);

    EC_xmax = max(EC.vertices(:,1));
    EC_xmin = min(EC.vertices(:,1));
    EC_ymax = max(EC.vertices(:,2));
    EC_ymin = min(EC.vertices(:,2));
    EC_zmax = max(EC.vertices(:,3));
    EC_zmin = min(EC.vertices(:,3));
    EC_width = abs(EC_xmax - EC_xmin);
    EC_length = abs(EC_ymax - EC_ymin);
    EC_height = abs(EC_zmax - EC_zmin);
    
    EC_center = [(EC_xmax + EC_xmin)/2 (EC_ymax + EC_ymin)/2 (EC_zmax + EC_zmin)/2]; 
    
    EC_SurfCenter = zeros(4,3);


    EC_SurfCenter(1,:) = [EC_xmax (EC_ymax + EC_ymin)/2 (EC_zmax + EC_zmin)/2]; %Counterclockwise

    EC_SurfCenter(2,:) = [EC_xmin (EC_ymax + EC_ymin)/2 (EC_zmax + EC_zmin)/2]; %Clockwise

    EC_SurfCenter(3,:) = [(EC_xmax + EC_xmin)/2 (EC_ymax + EC_ymin)/2 EC_zmax];  
    
    EC_SurfCenter(4,:) = [(EC_xmax + EC_xmin)/2 (EC_ymax + EC_ymin)/2 EC_zmin];  
    
    
end