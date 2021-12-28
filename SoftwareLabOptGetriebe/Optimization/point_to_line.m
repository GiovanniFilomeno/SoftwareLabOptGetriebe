function d = point_to_line(pt,v1,v2)
% Function that calculates the perpendicular or shortest distance between an axis and
% a point of a geometry. Used in :func:`OptimizationFunction` to calculate
% distance between the center of the IGBT and the Electric Motor axis.
%
% :param array pt: Point of the geometry which its distance to an axis will
%                  be computed
% :param array v1: First or starting point of the axis.
% :param array v2: Second or ending point of the axis
%
% :return: 
%   **d** : Perpendicular or shortest distance between the point and the axis
%
% :rtype: double
%
% **Example in Code**
%
% .. code-block::
%
%   distance = point_to_line(EC1_center, EM_start,EM_end);
%


      a = v1 - v2;
      b = pt - v2;

      d = norm(cross(a,b)) / norm(a);
end