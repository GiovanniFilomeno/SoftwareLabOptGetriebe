clear
clc

hold on

Cube=stlread("cube.STL");
V = Cube.Points;
V2 = zeros(26,3);
V2(:,1) = 1.3477.*V(:,1);
V2(:,2) = 0.58428.*V(:,2);
V2(:,3) = 0.7778513.*V(:,3);
V2 = 0.25.*V2;
V3 = 0.25.*V +80;
V4 = 0.2.*V;
V5 = 0.2.*V;

F = Cube.ConnectivityList;
N = faceNormal(Cube);

EC1 = struct('faces',F,'vertices',V2);

EC2 = struct('faces',F,'vertices',V3);


% patch(EC1,'FaceColor',       [0.8 0 0], ...
%     'EdgeColor',       'green',        ...
%     'FaceLighting',    'gouraud',     ...
%     'AmbientStrength', 0.15);
% 
% patch(EC2,'FaceColor',       [0 0.8 0], ...
%     'EdgeColor',       'red',        ...
%     'FaceLighting',    'gouraud',     ...
%     'AmbientStrength', 0.15);
% hold on
% view([45 90 135]);
% xlabel('X') 
% ylabel('Y')
% zlabel('Z')
% %axis([0 10 0 10 0 10])
% camlight('headlight');
% material('default'); 
EC1_xmax = max(EC1.vertices(:,1));
EC1_xmin = min(EC1.vertices(:,1));
EC1_ymax = max(EC1.vertices(:,2));
EC1_ymin = min(EC1.vertices(:,2));
EC1_zmax = max(EC1.vertices(:,3));
EC1_zmin = min(EC1.vertices(:,3));
EC1_length = abs(EC1_xmax - EC1_xmin);
EC1_width = abs(EC1_ymax - EC1_ymin);
EC1_height = abs(EC1_zmax - EC1_zmin);
EC1_center = [(EC1_xmax + EC1_xmin)/2 (EC1_ymax + EC1_ymin)/2 (EC1_zmax + EC1_zmin)/2];

EC1_center1 = [(EC1_xmax + EC1_xmin)/2 EC1_ymax (EC1_zmax + EC1_zmin)/2];
p2 = [(EC1_xmax + EC1_xmin)/2 EC1_ymax + 25 (EC1_zmax + EC1_zmin)/2];

% scatter3(EC1_center1(1),EC1_center1(2),EC1_center1(3),'o','green')


EC2_xmax = max(EC2.vertices(:,1));
EC2_xmin = min(EC2.vertices(:,1));
EC2_ymax = max(EC2.vertices(:,2));
EC2_ymin = min(EC2.vertices(:,2));
EC2_zmax = max(EC2.vertices(:,3));
EC2_zmin = min(EC2.vertices(:,3));
EC2_length = abs(EC2_xmax - EC2_xmin);
EC2_width = abs(EC2_ymax - EC2_ymin);
EC2_height = abs(EC2_zmax - EC2_zmin);
EC2_center = [(EC2_xmax + EC2_xmin)/2 (EC2_ymax + EC2_ymin)/2  (EC2_zmax + EC2_zmin)/2];

p3=[(EC2_xmax + EC2_xmin)/2 EC2_ymin-25 (EC2_zmax + EC2_zmin)/2];
EC2_center2 = [(EC2_xmax + EC2_xmin)/2 EC2_ymin (EC2_zmax + EC2_zmin)/2];

% scatter3(EC2_center2(1),EC2_center2(2),EC2_center2(3),'o','red')

% xyz = [EC1_center1; p2; p3; EC2_center2];
% plot3(xyz(:,1),xyz(:,2),xyz(:,3),'ro','LineWidth',2);
% % text(xyz,[repmat('  ',4,1), num2str((1:4)')])
% ax = gca;
% ax.XTick = [];
% ax.YTick = [];
% ax.ZTick = [];
% % xyz = [x ;y ;z];
% hold on
% Spline = cscvn(xyz');
% xyz_t = xyz'
% Spline.coefs
% Spline.breaks
% Spline.pieces
% 
% D_Spline = fnder(Spline,1)
% D_Spline.coefs
% D_Spline.breaks
% D_Spline.pieces
% 
% D_Spline = fnder(Spline,2)
% D_Spline.coefs
% D_Spline.breaks
% D_Spline.pieces
% 
% 
% 
% fn = Spline;
% fnplt(fn,'r',2);
% fnprime = fnder(fn);
% Lfun = @(s) sqrt(sum(fnval(fnprime,s).^2,1));
% L = integral(Lfun,fn.breaks(1),fn.breaks(end))
% 
cable_radius = 2.3;
radius_of_curvature = 13.8;

[length,max_curvature, curvature, collisionwith1, collisionwith2] = spline_connection (EC1_center1, p2, p3, EC2_center2,EC1_center, EC2_center, EC1, EC2, cable_radius, radius_of_curvature);
length
max_curvature
curvature;
collisionwith1;
collisionwith2;



% EC1_center2 = [(EC1_xmax + EC1_xmin)/2 EC1_ymin (EC1_zmax + EC1_zmin)/2]; 
% EC1_center3 = [(EC1_xmax + EC1_xmin)/2 (EC1_ymax + EC1_ymin)/2 EC1_zmax]; 
% 
% EC2_center1 = [(EC2_xmax + EC2_xmin)/2 EC2_ymax  (EC2_zmax + EC2_zmin)/2];



% x= [-1 0 2 3];
% y= [2 2 4 4];
% z= [9 9 -3 -3];
% 
% plot3(x,y, z,'ro','LineWidth',2);
% text(x,y,z,[repmat('  ',4,1), num2str((1:4)')])
% ax = gca;
% ax.XTick = [];
% ax.YTick = [];
% ax.ZTick = [];
% xyz = [x ;y ;z];
% hold on
% Spline = cscvn(xyz);
% Spline.coefs
% Spline.breaks
% Spline.pieces
% 
% D_Spline = fnder(Spline,2)
% D_Spline.coefs
% D_Spline.breaks
% D_Spline.pieces
% 
% 
% fnplt(cscvn(xyz),'r',2)
% hold off
% 
