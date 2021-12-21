function vertex = Rotation(V, anchor,indice, angle)
% Function used to rotate the vertices of a geometry along a specified
% coordinate axis with regard of an anchor point for a defined value of
% radians.
%
% :param double V: Vertices of the geometry to be rotated.
% :param double anchor: Rotation Anchor point.
% :param double indice: Axis of rotation(1=x, 2=y, 3=z).
% :param double angle: Value in radians of the desired rotation.
%
% :return: **V**: Array of vertices with the applied rotation.
%
% :rtype: double array
%
% **Example in Code**
%
% .. code-block:: 
%
%   gearbox.vertices = Rotation(gearbox.vertices, Anchor, 2, pi/2);
%

    Rz = [ cos(angle), -sin(angle), 0 ;
          sin(angle), cos(angle), 0 ;
    0, 0, 1 ];

    Ry = [ cos(angle), 0, sin(angle) ;
    0, 1, 0 ;
          -sin(angle), 0, cos(angle) ];
      
    Rx = [ 1, 0, 0 ;
    0, cos(angle), -sin(angle);
    0, sin(angle), cos(angle) ];

    if(indice==1)
           vertex = (V-anchor)*Rx + anchor;
    end
    if(indice==2)
           vertex = (V-anchor)*Ry + anchor;
    end
    if(indice==3)
           vertex = (V-anchor)*Rz + anchor;
    end

end 