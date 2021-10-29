function [length,max_curvature, curvature, collisionwith1, collisionwith2] = spline_connection (pt1, pt2, pt3, pt4,EC1_center, EC2_center, EC1, EC2, cable_radius, radius_of_curvature)

collisionwith1 = 1;
collisionwith2 = 1;

vector = pt3-pt2;
pt2a = pt2 + 0.25 .* vector ;
pt2b = pt2 + 0.5 .* vector ;
pt2c = pt2 + 0.75 .* vector ;

xyz = [pt1; pt2; pt2a; pt2b; pt2c; pt3; pt4];
patch(EC1,'FaceColor',       [0.8 0 0], ...
    'EdgeColor',       'green',        ...
    'FaceLighting',    'gouraud',     ...
    'AmbientStrength', 0.15);

patch(EC2,'FaceColor',       [0 0.8 0], ...
    'EdgeColor',       'red',        ...
    'FaceLighting',    'gouraud',     ...
    'AmbientStrength', 0.15);
hold on
view([45 90 135]);
xlabel('X') 
ylabel('Y')
zlabel('Z')
%axis([0 10 0 10 0 10])
camlight('headlight');
material('default');
% scatter3(pt1(1),pt1(2),pt1(3),'o','y');
% scatter3(pt4(1),pt4(2),pt4(3),'o','y');

% plot3(xyz(:,1),xyz(:,2),xyz(:,3),'yo','LineWidth',2);
% text(xyz,[repmat('  ',4,1), num2str((1:4)')])
ax = gca;
ax.XTick = [];
ax.YTick = [];
ax.ZTick = [];
% xyz = [x ;y ;z];
hold on  

iteration = 0;
max_curvature = 0.5;
while (collisionwith1>0 || collisionwith2>0 || max_curvature<(1/radius_of_curvature) || iteration < 10)
%xyz = [pt1; pt2; pt3; pt4];
iteration = iteration+1

vector1 = pt2a-EC1_center;
 vector1 = vector1 / norm(vector1);
vector2 = pt2c-EC2_center;
 vector2 = vector2 / norm(vector2);
 midptvector = EC2_center - EC1_center;
% vector3 = vector3 / norm(vector3);
%vector3 = point_to_line(pt2b, EC1_center, EC2_center);
ptcentervector = pt2b - EC1_center;
C = dot(ptcentervector, midptvector);
C = (C / (norm(midptvector))^2)*midptvector;
ptonline = EC1_center + C;
vector3 = pt2b - ptonline;
vector3 = vector3 / norm(vector3);

% pt2a = pt2a +[1 -1 1] ;
% pt2b = pt2b +[0 0 1]  ;
% pt2c = pt2c +[-1 -1 1] ;
pt2a = pt2a + 0.5*vector1 ;
pt2b = pt2b + 0.5*vector3  ;
pt2c = pt2c + 0.5*vector2 ;

xyz = [pt1; pt2; pt2a; pt2b; pt2c; pt3; pt4];

% plot3(xyz(:,1),xyz(:,2),xyz(:,3),'ro','LineWidth',2);
% % text(xyz,[repmat('  ',4,1), num2str((1:4)')])
% ax = gca;
% ax.XTick = [];
% ax.YTick = [];
% ax.ZTick = [];
% % xyz = [x ;y ;z];
% hold on

%spline generation 
Spline = cscvn(xyz');
xyz_t = xyz';
Spline.coefs;
Spline.breaks;
Spline.pieces;
Spline.dim;

%first derivative
D_Spline = fnder(Spline,1);
D_Spline.coefs;
D_Spline.breaks;
D_Spline.pieces;

%second derivative
DD_Spline = fnder(Spline,2);
DD_Spline.coefs;
DD_Spline.breaks;
DD_Spline.pieces;

%length of spline
fn = Spline;
% if iteration==1
%     fnplt(fn,'y',2);
% end
fnprime = fnder(fn);
Lfun = @(s) sqrt(sum(fnval(fnprime,s).^2,1));
length = integral(Lfun,fn.breaks(1),fn.breaks(end));


curvature = [];
collision1 = [];
collision2 = [];



    for i=1:Spline.pieces
        
             x = Spline.coefs(i*Spline.dim-2,:);
             y = Spline.coefs(i*Spline.dim-1,:);
             z = Spline.coefs(i*Spline.dim,:);
            D_x = D_Spline.coefs(i*Spline.dim-2,:);
            D_y = D_Spline.coefs(i*Spline.dim-1,:);
            D_z = D_Spline.coefs(i*Spline.dim,:);
            DD_x = DD_Spline.coefs(i*Spline.dim-2,:);
            DD_y = DD_Spline.coefs(i*Spline.dim-1,:);
            DD_z = DD_Spline.coefs(i*Spline.dim,:);
            
            t1 = Spline.breaks(i);
            t2 = Spline.breaks(i+1);
            increment = (t2-t1)/19;
            for t=t1:increment:t2
                 x_t =  x(1)*(t-Spline.breaks(i))^3 + x(2)*(t-Spline.breaks(i))^2 + x(3)*(t-Spline.breaks(i)) + x(4);
                 y_t =  y(1)*(t-Spline.breaks(i))^3 + y(2)*(t-Spline.breaks(i))^2 + y(3)*(t-Spline.breaks(i)) + y(4);
                 z_t =  z(1)*(t-Spline.breaks(i))^3 + z(2)*(t-Spline.breaks(i))^2 + z(3)*(t-Spline.breaks(i)) + z(4);
                D_x_t =  D_x(1)*(t-Spline.breaks(i))^2 + D_x(2)*(t-Spline.breaks(i)) + D_x(3);
                D_y_t =  D_y(1)*(t-Spline.breaks(i))^2 + D_y(2)*(t-Spline.breaks(i)) + D_y(3);
                D_z_t =  D_z(1)*(t-Spline.breaks(i))^2 + D_z(2)*(t-Spline.breaks(i)) + D_z(3);
                DD_x_t =  DD_x(1)*(t-Spline.breaks(i)) + DD_x(2);
                DD_y_t =  DD_y(1)*(t-Spline.breaks(i)) + DD_y(2);
                DD_z_t =  DD_z(1)*(t-Spline.breaks(i)) + DD_z(2);
                               
                k = sqrt(((DD_z_t)*(D_y_t) - (DD_y_t)*(D_z_t))^2 +((DD_x_t)*(D_z_t) - (DD_z_t)*(D_x_t))^2 +((DD_y_t)*(D_x_t) - (DD_x_t)*(D_y_t))^2 )/((D_x_t)^2 + (D_y_t)^2 +(D_z_t)^2)^1.5;
                curvature(end+1) = k; 
                
                
                collision1(end+1) = inpolyhedron(EC1, [x_t y_t z_t]);
                collision2(end+1) = inpolyhedron(EC2, [x_t y_t z_t]);
            end

    end

max_curvature = max(curvature); 
max_curvature_index = find(curvature == max_curvature);
%max_curvature_positions

collisionwith1 = sum(collision1)/size(collision1,2)
collisionwith2 = sum(collision2)/size(collision2,2)

    if (collisionwith1<0.01 && collisionwith2<0.01 && max_curvature>(1/radius_of_curvature)) || iteration == 20
        pt2a = pt2a + cable_radius*vector1 ;
        pt2b = pt2b + cable_radius*vector3  ;
        pt2c = pt2c + cable_radius*vector2 ;

        xyz = [pt1; pt2; pt2a; pt2b; pt2c; pt3; pt4];
        fn = cscvn(xyz');
        break;
    end


end
xyz = [pt1; pt2; pt2a; pt2b; pt2c; pt3; pt4];
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
% scatter3(pt1(1),pt1(2),pt1(3),'o','red');
% scatter3(pt4(1),pt4(2),pt4(3),'o','red');
fnplt(fn,'r',2);
plot3(xyz(:,1),xyz(:,2),xyz(:,3),'ro','LineWidth',2);
% text(xyz,[repmat('  ',4,1), num2str((1:4)')])
ax = gca;
ax.XTick = [];
ax.YTick = [];
ax.ZTick = [];
% xyz = [x ;y ;z];
hold on  
end


