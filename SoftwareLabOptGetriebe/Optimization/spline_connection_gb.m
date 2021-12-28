function [length,max_curvature_mid, curvature_mid, collisionwith1_mid, collisionwithgb_mid,fn_mid,fn_posy,fn_negy] = spline_connection_gb (pt1, pt2, pt3, pt4, pt5, pt6, pt7, pt8, pt9, pt10, pt11, pt12, pt13, pt14, pt15, EC1_center, EM_start, EM_end, fv, EC1,  cable_radius, radius_of_curvature)
% Function that creates a spline between the gearbox and the IGBT. Spline generation can be 
% successful or fail depending if it collides with one of the bodies. This result
% as well as the resultant spline length will have impact in the fitness
% value in the current individual in :func:`EvaluateFitness`.
%
% :param array pt1: Connection point of middle cable to Electric Motor
% :param array pt2: Second control point of middle cable ensuring perpendicularity
%                   of cable to the Electric Motor
% :param array pt3: Third control point of middle cable ensuring turnaround
%                   the motor without collision
% :param array pt4: Fourth control point of middle cable ensuring perpendicularity
%                   of cable to the component
% :param array pt5: Connection point of middle cable to component
% :param array pt6: Connection point of right cable to Electric Motor
% :param array pt7: Second control point of right cable ensuring perpendicularity
%                   of cable to the Electric Motor
% :param array pt8: Third control point of right cable ensuring turnaround
%                   the motor without collision
% :param array pt9: Fourth control point of right cable ensuring perpendicularity
%                   of cable to the component
% :param array pt10: Connection point of right cable to component 
% :param array pt11: Connection point of left cable to Electric Motor 
% :param array pt12: Second control point of left cable ensuring perpendicularity
%                   of cable to the Electric Motor
% :param array pt13: Third control point of left cable ensuring turnaround
%                   the motor without collision
% :param array pt14: Fourth control point of left cable ensuring perpendicularity
%                   of cable to the component
% :param array pt15: Connection point of left cable to component 
% :param array EC1_center: Vector of coordinates of the body center of the
%                          IGBT.
% :param array EM_start: Coordinates of the center of the electric motor.
% :param array EM_end: Coordinates of the end of the electric motor opposite to EM_start.
% :param struct fv: struct that represents the Gearbox geometry.
% :param struct EC1: struct that represents the IGBT geometry.
% :param double cale_radius: Radius of the prescribed cables
% :param double radius_of_curvature: Minimum allowable radius of curvature 
%                                    of the prescribed cables
%
% :return:
%   *[length,max_curvature, curvature, collisionwith1, collisionwith2, Penalty_spline, fn]*
%      - **length** : Length of the spline.
%      - **max_curvature_mid** : Max curvature of the spline.
%      - **curvature_mid** : Array with curvature values of segments of the
%        spline.
%      - **collisionwith1_mid** : Percentage of collision between spline and
%        first EC.
%      - **collisionwithgb_mid** : Percentage of collision between spline and
%        second EC.
%      - **fn_mid** : Function that describes the middle spline so it can be plotted
%        with the :func:`Visualize`.
%      - **fn_posy** : Function that describes the left spline so it can be plotted
%        with the :func:`Visualize`.
%      - **fn_negy** : Function that describes the right spline so it can be plotted
%        with the :func:`Visualize`.
%
%
% **Example in Code**
%
% .. code-block::
%
%    [splinelength,max_curvature, curvature, collisionwith1, collisionwithgb,fn_mid,fn_posy,fn_negy] = spline_connection_gb (EM_connection, EM_connection+[0 -clearance_mid 0], pt3, EC1_SurfCenter(2,:)+(clearance_mid.*vector1)/norm(vector1), EC1_SurfCenter(2,:), pt6, pt7, pt8, pt9, pt10, pt11, pt12, pt13, pt14, pt15,EC1_center, EM_start, EM_end, gearbox, EC1,  cable_radius, radius_of_curvature);
% 
            
            
            
            
% initial set values for collision
collisionwith1_mid = 1;
collisionwithgb_mid = 1;
% determining intermediate control points for middle spline
vector_mid = pt4-pt3;
pt3a = pt3 + 0.25 .* vector_mid ;
pt3b = pt3 + 0.5 .* vector_mid ;
pt3c = pt3 + 0.75 .* vector_mid ;
% determining intermediate control points for right spline
vector_posy = pt9-pt8;
pt8a = pt8 + 0.25 .* vector_posy ;
pt8b = pt8 + 0.5 .* vector_posy ;
pt8c = pt8 + 0.75 .* vector_posy ;
% determining intermediate control points for left spline
vector_negy = pt14-pt13;
pt13a = pt13 + 0.25 .* vector_negy ;
pt13b = pt13 + 0.5 .* vector_negy;
pt13c = pt13 + 0.75 .* vector_negy ;
% control points vectors
xyz_mid = [pt1; pt2; pt3; pt3a; pt3b; pt3c; pt4; pt5];
xyz_posy = [pt6; pt7; pt8; pt8a; pt8b; pt8c; pt9; pt10];
xyz_negy = [pt11; pt12; pt13; pt13a; pt13b; pt13c; pt14; pt15];

ax = gca;
ax.XTick = [];
ax.YTick = [];
ax.ZTick = [];

hold on  

EM_connection = EM_end;
% initial iteration value
iteration = 0;
% set initial maximum curvature
max_curvature_mid = 0.5;
% Enter iteration loop in case of collisions or curvature violation
while (collisionwith1_mid>0.0072 || collisionwithgb_mid>0.0072 || max_curvature_mid<(1/radius_of_curvature) || iteration < 10)

% iteration update
iteration = iteration+1;
% centre of electric motor
EM_center = EM_start + 0.5*(EM_end-EM_start);
% vector direction to shift control point 3a in direction parallel to line
% joining point and the electric motor centre
vector1 = pt3a-EM_center;
vector1 = vector1 / norm(vector1);
% vector direction to shift control point 3c in direction parallel to line
% joining point and the component centre
vector2 = pt3c-EC1_center;
vector2 = vector2 / norm(vector2);
% vector direction to shift control point 3b in direction perpendicular to
% line joining the centres of the electric mototr and component
midptvector = EC1_center - EM_connection;
ptcentervector = pt3b - EM_connection;
C = dot(ptcentervector, midptvector);
C = (C / (norm(midptvector))^2)*midptvector;
ptonline = EM_connection + C;
vector3 = pt3b - ptonline;
vector3 = vector3 / norm(vector3);

% Update the control points locations of middle spline
pt3a = pt3a + 0.5*vector1 ;
pt3b = pt3b + 0.5*vector3  ;
pt3c = pt3c + 0.5*vector2 ;
% Update the control points locations of right spline
pt8a = pt8a + 0.5*vector1 ;
pt8b = pt8b + 0.5*vector3  ;
pt8c = pt8c + 0.5*vector2 ;
% Update the control points locations of left spline
pt13a = pt13a + 0.5*vector1 ;
pt13b = pt13b + 0.5*vector3  ;
pt13c = pt13c + 0.5*vector2 ;
% updated control points vectors
xyz_mid = [pt1; pt2; pt3; pt3a; pt3b; pt3c; pt4; pt5];
xyz_posy = [pt6; pt7; pt8; pt8a; pt8b; pt8c; pt9; pt10];
xyz_negy = [pt11; pt12; pt13; pt13a; pt13b; pt13c; pt14; pt15];

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

fnprime = fnder(fn_mid);
Lfun = @(s) sqrt(sum(fnval(fnprime,s).^2,1));
length = integral(Lfun,fn_mid.breaks(1),fn_mid.breaks(end));

% Vectors to store collision and curvature values
curvature_mid = [];
collision1_mid = [];
collisiongb_mid = [];


% Discretize the spline to check for collision and curvature
    for i=1:Spline_mid.pieces
             %Obtaining the cartesian coordinates from parametric form
             x = Spline_mid.coefs(i*Spline_mid.dim-2,:);
             y = Spline_mid.coefs(i*Spline_mid.dim-1,:);
             z = Spline_mid.coefs(i*Spline_mid.dim,:);
             %Obtaining the first derivative from parametric form
            D_x = D_Spline_mid.coefs(i*Spline_mid.dim-2,:);
            D_y = D_Spline_mid.coefs(i*Spline_mid.dim-1,:);
            D_z = D_Spline_mid.coefs(i*Spline_mid.dim,:);
            %Obtaining the second derivative from parametric form
            DD_x = DD_Spline_mid.coefs(i*Spline_mid.dim-2,:);
            DD_y = DD_Spline_mid.coefs(i*Spline_mid.dim-1,:);
            DD_z = DD_Spline_mid.coefs(i*Spline_mid.dim,:);
            
            % discretization
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
                % curvature calculation               
                k = sqrt(((DD_z_t)*(D_y_t) - (DD_y_t)*(D_z_t))^2 +((DD_x_t)*(D_z_t) - (DD_z_t)*(D_x_t))^2 +((DD_y_t)*(D_x_t) - (DD_x_t)*(D_y_t))^2 )/((D_x_t)^2 + (D_y_t)^2 +(D_z_t)^2)^1.5;
                curvature_mid(end+1) = k; 
                
                % collision check at each point
                collision1_mid(end+1) = inpolyhedron(EC1, [x_t y_t z_t]);
                collisiongb_mid(end+1) = inpolyhedron(fv, [x_t y_t z_t]);
            end

    end
% saving curvature values
max_curvature_mid = max(curvature_mid); 
%max_curvature_index = find(curvature_mid == max_curvature_mid);

% total collision
collisionwith1_mid = sum(collision1_mid)/size(collision1_mid,2);
collisionwithgb_mid = sum(collisiongb_mid)/size(collisiongb_mid,2);

% Loop exit in case of no collisions and allowable curvature or maximum iterations
    if (collisionwith1_mid<=0.0072 && collisionwithgb_mid<=0.0072 && max_curvature_mid>(1/radius_of_curvature)) || iteration == 10
        % Accomodating cable radius
        pt3a = pt3a + cable_radius*vector1 ;
        pt3b = pt3b + cable_radius*vector3  ;
        pt3c = pt3c + cable_radius*vector2 ;
        % final middle spline control points 
        xyz_mid = [pt1; pt2; pt3; pt3a; pt3b; pt3c; pt4; pt5];
        fn_mid = cscvn(xyz_mid');
        % Accomodating cable radius
        pt8a = pt8a + cable_radius*vector1 ;
        pt8b = pt8b + cable_radius*vector3  ;
        pt8c = pt8c + cable_radius*vector2 ;
        % final right spline control points 
        xyz_posy = [pt6; pt7; pt8; pt8a; pt8b; pt8c; pt9; pt10];
        fn_posy = cscvn(xyz_posy');
        % Accomodating cable radius
        pt13a = pt13a + cable_radius*vector1 ;
        pt13b = pt13b + cable_radius*vector3  ;
        pt13c = pt13c + cable_radius*vector2 ;
        % final left spline control points 
        xyz_negy = [pt11; pt12; pt13; pt13a; pt13b; pt13c; pt14; pt15];
        fn_negy = cscvn(xyz_negy');
        break;
    end


end
% final middle spline control points
xyz_mid = [pt1; pt2; pt3; pt3a; pt3b; pt3c; pt4; pt5];

ax = gca;
ax.XTick = [];
ax.YTick = [];
ax.ZTick = [];

end


