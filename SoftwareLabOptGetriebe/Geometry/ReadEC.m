function [EC,EC_length,EC_width,EC_height,EC_center, EC_SurfCenter] = ReadEC(file)
% Function that reads the Electric Component(EC) STL file and returns a struct with the
% vertices and faces of the geometry, as well as dimensions and centers of
% the geometry
%
% :param str file: The name of the STL file located in the path folder
%
% :return: 
%   *[EC, EC_length, EC_width, EC_height, EC_center, EC_SurfCenter]*
%       - EC: Struct of the geometry with EC.faces and EC.vertices.
%       - EC_length: Length of the component.
%       - EC_width: Width of the component.
%       - EC_height: Height of the component. 
%       - EC_center: Vector of coordinates of the body center of the
%         component, used for rotations of the component.
%       - EC_SurfCenter: Array of vectors with coordinates of surface centers of
%         the component, used to calculate distance between components.
%
% :rtype: [struct, double, double, double, double, double array, double array] 
%
% **Example in Code**
%
% .. code-block:: 
%
%   [EC, EC_length, EC_width, EC_height, EC_center, EC_SurfCenter] = ReadEC('EC.stl');
%
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