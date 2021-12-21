function [length,max_curvature_mid, curvature_mid, collisionwith1_mid, collisionwithgb_mid,fn_mid,fn_posy,fn_negy] = spline_connection_gb (pt1, pt2, pt3, pt4, pt5, pt6, pt7, pt8, pt9, pt10, pt11, pt12, pt13, pt14, pt15, EC1_center, EM_start, EM_end, fv, EC1,  cable_radius, radius_of_curvature)
% Function that creates a spline between the gearbox and the IGBT. Spline generation can be 
% successful or fail depending if it collides with one of the bodies. This result
% as well as the resultant spline length will have impact in the fitness
% value in the current individual in :func:`EvaluateFitness`.
%
% :param array pt1: 
% :param array pt2: 
% :param array pt3: 
% :param array pt4: 
% :param array pt5: 
% :param array pt6: 
% :param array pt7: 
% :param array pt8: 
% :param array pt9: 
% :param array pt10: 
% :param array pt11: 
% :param array pt12: 
% :param array pt13: 
% :param array pt14: 
% :param array pt15: 
% :param array EC1_center: Vector of coordinates of the body center of the
%                          IGBT.
% :param array EM_start: Coordinates of the center of the electric motor.
% :param array EM_end: Coordinates of the end of the electric motor opposite to EM_start.
% :param struct fv: struct that represents the Gearbox geometry.
% :param struct EC1: struct that represents the IGBT geometry.
% :param double cale_radius:
% :param double radius_of_curvature:
%
% :return:
%   *[length,max_curvature, curvature, collisionwith1, collisionwith2, Penalty_spline, fn]*
%      - **length** :
%      - **max_curvature** :
%      - **curvature** :
%      - **collisionwith1** :
%      - **collisionwith2** :
%      - **Penalty_spline,** :
%      - **fn** :
%
%
% **Example in Code**
%
% .. code-block::
%
%    [splinelength,max_curvature, curvature, collisionwith1, collisionwithgb,fn_mid,fn_posy,fn_negy] = spline_connection_gb (EM_connection, EM_connection+[0 -clearance_mid 0], pt3, EC1_SurfCenter(2,:)+(clearance_mid.*vector1)/norm(vector1), EC1_SurfCenter(2,:), pt6, pt7, pt8, pt9, pt10, pt11, pt12, pt13, pt14, pt15,EC1_center, EM_start, EM_end, gearbox, EC1,  cable_radius, radius_of_curvature);
% 
            
            
            
            

collisionwith1_mid = 1;
collisionwithgb_mid = 1;

vector_mid = pt4-pt3;
pt3a = pt3 + 0.25 .* vector_mid ;
pt3b = pt3 + 0.5 .* vector_mid ;
pt3c = pt3 + 0.75 .* vector_mid ;

vector_posy = pt9-pt8;
pt8a = pt8 + 0.25 .* vector_posy ;
pt8b = pt8 + 0.5 .* vector_posy ;
pt8c = pt8 + 0.75 .* vector_posy ;

vector_negy = pt14-pt13;
pt13a = pt13 + 0.25 .* vector_negy ;
pt13b = pt13 + 0.5 .* vector_negy;
pt13c = pt13 + 0.75 .* vector_negy ;

xyz_mid = [pt1; pt2; pt3; pt3a; pt3b; pt3c; pt4; pt5];
xyz_posy = [pt6; pt7; pt8; pt8a; pt8b; pt8c; pt9; pt10];
xyz_negy = [pt11; pt12; pt13; pt13a; pt13b; pt13c; pt14; pt15];

% patch(EC1,'FaceColor',       [0.8 0 0], ...
%     'EdgeColor',       'green',        ...
%     'FaceLighting',    'gouraud',     ...
%     'AmbientStrength', 0.15);
% patch(fv,'FaceColor',       [0.8 0.8 1], ...
%                  'FaceLighting',    'gouraud',     ...
%                  'AmbientStrength', 0.15);
             
%plot3(shaftx', shafty', shaftz','LineWidth',2,'Color','red')
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
% scatter3(pt1(1),pt1(2),pt1(3),'o','y');
% scatter3(pt4(1),pt4(2),pt4(3),'o','y');

% plot3(xyz_mid(:,1),xyz_mid(:,2),xyz_mid(:,3),'yo','LineWidth',2);
% plot3(xyz_posy(:,1),xyz_posy(:,2),xyz_posy(:,3),'yo','LineWidth',2);
% plot3(xyz_negy(:,1),xyz_negy(:,2),xyz_negy(:,3),'yo','LineWidth',2);
% text(xyz,[repmat('  ',4,1), num2str((1:4)')])
ax = gca;
ax.XTick = [];
ax.YTick = [];
ax.ZTick = [];
% xyz = [x ;y ;z];
hold on  
EM_connection = EM_end;
iteration = 0;
max_curvature_mid = 0.5;
while (collisionwith1_mid>0.0072 || collisionwithgb_mid>0.0072 || max_curvature_mid<(1/radius_of_curvature) || iteration < 10)
%xyz = [pt1; pt2; pt3; pt4];
iteration = iteration+1;
EM_center = EM_start + 0.5*(EM_end-EM_start);
vector1 = pt3a-EM_center;
vector1 = vector1 / norm(vector1);
vector2 = pt3c-EC1_center;
vector2 = vector2 / norm(vector2);
midptvector = EC1_center - EM_connection;
% vector3 = vector3 / norm(vector3);
%vector3 = point_to_line(pt2b, EC1_center, EC2_center);
ptcentervector = pt3b - EM_connection;
C = dot(ptcentervector, midptvector);
C = (C / (norm(midptvector))^2)*midptvector;
ptonline = EM_connection + C;
vector3 = pt3b - ptonline;
vector3 = vector3 / norm(vector3);

% pt2a = pt2a +[1 -1 1] ;
% pt2b = pt2b +[0 0 1]  ;
% pt2c = pt2c +[-1 -1 1] ;
pt3a = pt3a + 0.5*vector1 ;
pt3b = pt3b + 0.5*vector3  ;
pt3c = pt3c + 0.5*vector2 ;

pt8a = pt8a + 0.5*vector1 ;
pt8b = pt8b + 0.5*vector3  ;
pt8c = pt8c + 0.5*vector2 ;

pt13a = pt13a + 0.5*vector1 ;
pt13b = pt13b + 0.5*vector3  ;
pt13c = pt13c + 0.5*vector2 ;

xyz_mid = [pt1; pt2; pt3; pt3a; pt3b; pt3c; pt4; pt5];
xyz_posy = [pt6; pt7; pt8; pt8a; pt8b; pt8c; pt9; pt10];
xyz_negy = [pt11; pt12; pt13; pt13a; pt13b; pt13c; pt14; pt15];
% plot3(xyz(:,1),xyz(:,2),xyz(:,3),'ro','LineWidth',2);
% % text(xyz,[repmat('  ',4,1), num2str((1:4)')])
% ax = gca;
% ax.XTick = [];
% ax.YTick = [];
% ax.ZTick = [];
% % xyz = [x ;y ;z];
% hold on

%spline generation 
Spline_mid = cscvn(xyz_mid');
Spline_mid.coefs;
Spline_mid.breaks;
Spline_mid.pieces;
Spline_mid.dim;

%first derivative
D_Spline_mid = fnder(Spline_mid,1);
D_Spline_mid.coefs;
D_Spline_mid.breaks;
D_Spline_mid.pieces;

%second derivative
DD_Spline_mid = fnder(Spline_mid,2);
DD_Spline_mid.coefs;
DD_Spline_mid.breaks;
DD_Spline_mid.pieces;

%length of spline
fn_mid = Spline_mid;
fn_posy = cscvn(xyz_posy');
fn_negy = cscvn(xyz_negy');

% if iteration==1
%     fnplt(fn_mid,'y',2);
%     fnplt(fn_posy,'y',2);
%     fnplt(fn_negy,'y',2);
% end
fnprime = fnder(fn_mid);
Lfun = @(s) sqrt(sum(fnval(fnprime,s).^2,1));
length = integral(Lfun,fn_mid.breaks(1),fn_mid.breaks(end));


curvature_mid = [];
collision1_mid = [];
collisiongb_mid = [];



    for i=1:Spline_mid.pieces
        
             x = Spline_mid.coefs(i*Spline_mid.dim-2,:);
             y = Spline_mid.coefs(i*Spline_mid.dim-1,:);
             z = Spline_mid.coefs(i*Spline_mid.dim,:);
            D_x = D_Spline_mid.coefs(i*Spline_mid.dim-2,:);
            D_y = D_Spline_mid.coefs(i*Spline_mid.dim-1,:);
            D_z = D_Spline_mid.coefs(i*Spline_mid.dim,:);
            DD_x = DD_Spline_mid.coefs(i*Spline_mid.dim-2,:);
            DD_y = DD_Spline_mid.coefs(i*Spline_mid.dim-1,:);
            DD_z = DD_Spline_mid.coefs(i*Spline_mid.dim,:);
            
            t1 = Spline_mid.breaks(i);
            t2 = Spline_mid.breaks(i+1);
            increment = (t2-t1)/19;
            for t=t1:increment:t2
                 x_t =  x(1)*(t-Spline_mid.breaks(i))^3 + x(2)*(t-Spline_mid.breaks(i))^2 + x(3)*(t-Spline_mid.breaks(i)) + x(4);
                 y_t =  y(1)*(t-Spline_mid.breaks(i))^3 + y(2)*(t-Spline_mid.breaks(i))^2 + y(3)*(t-Spline_mid.breaks(i)) + y(4);
                 z_t =  z(1)*(t-Spline_mid.breaks(i))^3 + z(2)*(t-Spline_mid.breaks(i))^2 + z(3)*(t-Spline_mid.breaks(i)) + z(4);
                D_x_t =  D_x(1)*(t-Spline_mid.breaks(i))^2 + D_x(2)*(t-Spline_mid.breaks(i)) + D_x(3);
                D_y_t =  D_y(1)*(t-Spline_mid.breaks(i))^2 + D_y(2)*(t-Spline_mid.breaks(i)) + D_y(3);
                D_z_t =  D_z(1)*(t-Spline_mid.breaks(i))^2 + D_z(2)*(t-Spline_mid.breaks(i)) + D_z(3);
                DD_x_t =  DD_x(1)*(t-Spline_mid.breaks(i)) + DD_x(2);
                DD_y_t =  DD_y(1)*(t-Spline_mid.breaks(i)) + DD_y(2);
                DD_z_t =  DD_z(1)*(t-Spline_mid.breaks(i)) + DD_z(2);
                               
                k = sqrt(((DD_z_t)*(D_y_t) - (DD_y_t)*(D_z_t))^2 +((DD_x_t)*(D_z_t) - (DD_z_t)*(D_x_t))^2 +((DD_y_t)*(D_x_t) - (DD_x_t)*(D_y_t))^2 )/((D_x_t)^2 + (D_y_t)^2 +(D_z_t)^2)^1.5;
                curvature_mid(end+1) = k; 
                
                
                collision1_mid(end+1) = inpolyhedron(EC1, [x_t y_t z_t]);
                collisiongb_mid(end+1) = inpolyhedron(fv, [x_t y_t z_t]);
            end

    end

max_curvature_mid = max(curvature_mid); 
%max_curvature_index = find(curvature_mid == max_curvature_mid);
%max_curvature_positions

collisionwith1_mid = sum(collision1_mid)/size(collision1_mid,2);
collisionwithgb_mid = sum(collisiongb_mid)/size(collisiongb_mid,2);

    if (collisionwith1_mid<=0.0072 && collisionwithgb_mid<=0.0072 && max_curvature_mid>(1/radius_of_curvature)) || iteration == 10
        pt3a = pt3a + cable_radius*vector1 ;
        pt3b = pt3b + cable_radius*vector3  ;
        pt3c = pt3c + cable_radius*vector2 ;

        xyz_mid = [pt1; pt2; pt3; pt3a; pt3b; pt3c; pt4; pt5];
        fn_mid = cscvn(xyz_mid');
        
        pt8a = pt8a + cable_radius*vector1 ;
        pt8b = pt8b + cable_radius*vector3  ;
        pt8c = pt8c + cable_radius*vector2 ;
        xyz_posy = [pt6; pt7; pt8; pt8a; pt8b; pt8c; pt9; pt10];
        fn_posy = cscvn(xyz_posy');
        
        pt13a = pt13a + cable_radius*vector1 ;
        pt13b = pt13b + cable_radius*vector3  ;
        pt13c = pt13c + cable_radius*vector2 ;
        xyz_negy = [pt11; pt12; pt13; pt13a; pt13b; pt13c; pt14; pt15];
        fn_negy = cscvn(xyz_negy');
        break;
    end


end

xyz_mid = [pt1; pt2; pt3; pt3a; pt3b; pt3c; pt4; pt5];
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
% fnplt(fn_mid,'r',2);
% plot3(xyz_mid(:,1),xyz_mid(:,2),xyz_mid(:,3),'ro','LineWidth',2);
% fnplt(fn_posy,'r',2);
% plot3(xyz_posy(:,1),xyz_posy(:,2),xyz_posy(:,3),'go','LineWidth',2);
% fnplt(fn_negy,'r',2);
% plot3(xyz_negy(:,1),xyz_negy(:,2),xyz_negy(:,3),'bo','LineWidth',2);
% text(xyz,[repmat('  ',4,1), num2str((1:4)')])
ax = gca;
ax.XTick = [];
ax.YTick = [];
ax.ZTick = [];
% xyz = [x ;y ;z];
% hold on  
end


