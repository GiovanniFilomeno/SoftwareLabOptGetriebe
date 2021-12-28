function gearbox = ReadGearbox(file)
% Function that reads the Gearbox STL file and returns a struct with the
% vertices and faces of the geometry
%
% :param str file: The name of the STL file located in the path folder
%
% :return: 
%   *gearbox*   
%       - gearbox.vertices
%       - gearbox.faces
%
% :rtype: struct
%
% **Example in Code**
%
% .. code-block:: 
%
%   gearbox = ReadGearbox('Gearbox1.stl');
%



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