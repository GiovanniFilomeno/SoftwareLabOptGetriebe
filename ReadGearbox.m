function gearbox = ReadGearbox(file)

    if ~exist(file,'file')
        error(['File ''%s'' not found. If the file is not on MATLAB''s path' ...
               ', be sure to specify the full path to the file.'], file);
    end
    
    GEARBOX = stlread(file);
    V = GEARBOX.Points;
    V = 1000.*V;
    F = GEARBOX.ConnectivityList;
    N = faceNormal(GEARBOX);
    gearbox = struct('faces',F,'vertices',V);
end