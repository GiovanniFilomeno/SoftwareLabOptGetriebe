clear
clc

hold on


%
%Read Gearbox
FV=stlread("Gearbox1.STL");
V = FV.Points;
V1 = 1000.*V;
F = FV.ConnectivityList;
N = faceNormal(FV);
fv = struct('faces',F,'vertices',V1);

%Read Assembly space
BR=stlread("Bauraum.STL");
V = BR.Points;
F = BR.ConnectivityList;
N = faceNormal(BR);
bauraum = struct('faces',F,'vertices',V);

%Read 4 Electric components and set Right dimensions
Cube=stlread("cube.STL");
V = Cube.Points;
V2 = zeros(26,3);
% V2(:,1) = 1.3477.*V(:,1);
% V2(:,2) = 0.58428.*V(:,2);
% V2(:,3) = 0.7778513.*V(:,3);
V2(:,2) = 0.566.*V(:,1);
V2(:,1) = 0.2045.*V(:,2);
V2(:,3) = 0.07629.*V(:,3);
V2 = V2 +[2300 -250 220];
% V2 = 0.25.*V2;
% V3 = 0.25.*V +80;
% V4 = 0.2.*V;
% V5 = 0.2.*V;

F = Cube.ConnectivityList;
N = faceNormal(Cube);

EC1 = struct('faces',F,'vertices',V2);

% EC2 = struct('faces',F,'vertices',V3);

%Plot Bauraum
BRaum = patch(bauraum,'FaceColor',       [0 0 1.0], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15);
%BRaum.FaceVertexAlphaData = 0.02;    % Set constant transparency 
BRaum.FaceAlpha = 0.2 ;          % Interpolate to find face transparency
hold on
view([45 90 135]);
xlabel('X') 
ylabel('Y')
zlabel('Z')
%axis([0 10 0 10 0 10])
camlight('headlight');
material('default');
% 
% 
% patch(EC1,'FaceColor',       [0.8 0 0], ...
%     'EdgeColor',       'green',        ...
%     'FaceLighting',    'gouraud',     ...
%     'AmbientStrength', 0.15);
% % 
% % patch(EC2,'FaceColor',       [0 0.8 0], ...
% %     'EdgeColor',       'red',        ...
% %     'FaceLighting',    'gouraud',     ...
% %     'AmbientStrength', 0.15);
% hold on
% view([45 90 135]);
% xlabel('X') 
% ylabel('Y')
% zlabel('Z')
% %axis([0 10 0 10 0 10])
% camlight('headlight');
% material('default'); 

%% Positioning(manually) of the gearbox inside the bauraum
[EM_start, EM_vector, Anchor, EM_radius] = FindingAxis("Gearbox1.stl");
EM_end = EM_start + EM_vector;
zero = [0 0 0];

%Definition(manually of the axis of the bauraum shaft 
point = [2539 14 -40];
vector = [0 -1 0];

%setting rotation angle values
if (EM_vector(1) ~= 0)
    if (EM_vector(1) > 0)
        anglex = 0;
        angley = 0;
        anglez = pi/2;
    else
        anglex = 0;
        angley = 0;
        anglez = -pi/2;
    end
end
    
if (EM_vector(2) ~= 0)
    if (EM_vector(2) > 0)
        anglex = 0;
        angley = 0;
        anglez = 0;
    else
        anglex = 0;
        angley = 0;
        anglez = 0;
    end
end
        
if (EM_vector(3) ~= 0)
    if (EM_vector(3) > 0)
        anglex = 0;
        angley = pi/2;
        anglez = -pi/2;
    else
        anglex = 0;
        angley = -pi/2;
        anglez = pi/2;
    end
end

%Positioning of the gearbox inside the Bauraum
fv.vertices = fv.vertices - Anchor + point;
EM_start = EM_start - Anchor + point;
Anchor = point;
shaftx = [2539 2539];
shafty = [14 -130];
shaftz = [-40 -40];
plot3(shaftx', shafty', shaftz','LineWidth',2,'Color','red')
hold on
pause(1)
fv.vertices = rotation(fv.vertices,Anchor, 1, anglex);
fv.vertices = rotation(fv.vertices,Anchor, 2, angley);
fv.vertices = rotation(fv.vertices,Anchor, 3, anglez);
EM_start = rotation(EM_start, Anchor, 1, anglex);
EM_start = rotation(EM_start, Anchor, 2, angley);
EM_start = rotation(EM_start, Anchor, 3, anglez);
EM_vector = rotation(EM_vector, zero, 1, anglex);
EM_vector = rotation(EM_vector, zero, 2, angley);
EM_vector = rotation(EM_vector, zero, 3, anglez);
EM_end = EM_start + EM_vector;

%Obtaining the angle range for the gearbox to rotate without escaping the
%bauraum
out_gearbox = [];
for i= 1:360
    
    fv.vertices = rotation(fv.vertices,Anchor, 2, pi/180);
    EM_start = rotation(EM_start, Anchor, 2, pi/180);
    EM_vector = rotation(EM_vector, zero, 2, pi/180);
    EM_end = EM_start + EM_vector;
    EM_shaftx = [EM_start(1,1) EM_start(1,1) + EM_vector(1,1)];
    EM_shafty = [EM_start(1,2) EM_start(1,2) + EM_vector(1,2)];
    EM_shaftz = [EM_start(1,3) EM_start(1,3) + EM_vector(1,3)];
    
    IN = inpolyhedron(bauraum,fv.vertices);
    out_gearbox(i) = sum(IN)/length(IN); 
end

valid_angles = [];
max_IN = max(out_gearbox)


for i= 1:360
    if (out_gearbox(i) == max_IN)
        valid_angles(end+1) = i;
        valid_radians = valid_angles.*(pi/180);
    end
end
out_gearbox = [];
angle_range = valid_angles(length(valid_angles)) - valid_angles(1);
rad_range = (length(valid_angles))*pi/180;

if angle_range == length(valid_angles) - 1
    first_angle = valid_radians(1);
else
    for i=1:length(valid_radians)
        if i==1
            valid_diff(i) = valid_radians(1) - valid_radians(length(valid_radians)); 
        else
            valid_diff(i) = valid_radians(i) - valid_radians(i-1);
        end
    end
    [srt,I] = sort(valid_diff, 'descend');
    first_angle = valid_radians(I(1,1));
end

fv.vertices = rotation(fv.vertices,Anchor, 2, first_angle);
EM_start = rotation(EM_start, Anchor, 2, first_angle);
EM_end = EM_start + EM_vector;
%% 
% 
% patch(fv,'FaceColor',       [0.8 0.8 1], ...
%                  'FaceLighting',    'gouraud',     ...
%                  'AmbientStrength', 0.15);
             
plot3(shaftx', shafty', shaftz','LineWidth',2,'Color','red')
            
%% 



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

%EC1_center1 = [(EC1_xmax + EC1_xmin)/2 EC1_ymax (EC1_zmax + EC1_zmin)/2];
EC1_center1 = [EC1_xmax  (EC1_ymax+EC1_ymin)/2 (EC1_zmax + EC1_zmin)/2]; 
%p2 = [(EC1_xmax + EC1_xmin)/2 EC1_ymax + 25 (EC1_zmax + EC1_zmin)/2];

EC1_center2 = [(EC1_xmax + EC1_xmin)/2 EC1_ymin (EC1_zmax + EC1_zmin)/2];
% scatter3(EC1_center1(1),EC1_center1(2),EC1_center1(3),'o','green')
EM_connection = EM_end;
clearance = 20;
rad_gb = EM_radius;
vector_pt1 = [EC1_center1(1)+(clearance) EM_connection(2) EC1_center1(3)];
vector_pt2 = EM_connection+[0 -clearance 0]; 
vector = vector_pt1 - vector_pt2;
vector = vector / norm(vector);
            %pt3 = [(EC1_xmax + EC1_xmin)/2 EM_connection(2)+(-clearance) EM_connection(3)+(50)];
pt3 = EM_connection+[0 -clearance 0]+ rad_gb*vector;
im = [EM_connection; EM_connection+[0 -clearance 0]; pt3; EC1_center2+[clearance 0 0]; EC1_center2];



% EC2_xmax = max(EC2.vertices(:,1));
% EC2_xmin = min(EC2.vertices(:,1));
% EC2_ymax = max(EC2.vertices(:,2));
% EC2_ymin = min(EC2.vertices(:,2));
% EC2_zmax = max(EC2.vertices(:,3));
% EC2_zmin = min(EC2.vertices(:,3));
% EC2_length = abs(EC2_xmax - EC2_xmin);
% EC2_width = abs(EC2_ymax - EC2_ymin);
% EC2_height = abs(EC2_zmax - EC2_zmin);
% EC2_center = [(EC2_xmax + EC2_xmin)/2 (EC2_ymax + EC2_ymin)/2  (EC2_zmax + EC2_zmin)/2];
% 
% p3=[(EC2_xmax + EC2_xmin)/2 EC2_ymin-25 (EC2_zmax + EC2_zmin)/2];
% EC2_center2 = [(EC2_xmax + EC2_xmin)/2 EC2_ymin (EC2_zmax + EC2_zmin)/2];

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

[length,max_curvature, curvature, collisionwith1, collisionwithgb] = spline_connection_gb (EM_connection, EM_connection+[0 -clearance 0], pt3, EC1_center1+[clearance 0 0], EC1_center1,EC1_center, EM_start, EM_end, fv, EC1);
length
max_curvature
curvature;
collisionwith1;
collisionwithgb;



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
function d = point_to_line(pt,v1,v2)
      a = v1 - v2;
      b = pt - v2;

      d = norm(cross(a,b)) / norm(a);
end

function vertex = rotation(V, anchor,indice, angle)

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