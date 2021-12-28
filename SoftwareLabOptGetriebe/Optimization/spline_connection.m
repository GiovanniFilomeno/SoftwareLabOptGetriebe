function [length,max_curvature, curvature, collisionwith1, collisionwith2, Penalty_spline, fn] = spline_connection (pt1, pt2, pt3, pt4,EC1_center, EC2_center, EC1, EC2, bauraum, cable_radius, radius_of_curvature)
% Function that creates a spline (cable) between two Electrical components. Spline generation can be 
% successful or unsuccesful depending on whether it collides with one of the bodies. This result
% as well as the resultant spline length will have an impact on the fitness
% value of the current individual in :func:`EvaluateFitness`.
%
% :param array pt1: Connection point to component 1      
% :param array pt2: Second control point of cable ensuring perpendicularity
%                   of cable to the component 1
% :param array pt3: Third control point of cable ensuring perpendicularity
%                   of cable to the component 2
% :param array pt4: Connection point to component 2
% :param array EC1_center: Vector of coordinates of the body center of the
%                          IGBT.
% :param array EC2_center: Vector of coordinates of the body center of the
%                          Capacitor.
% :param struct EC1: struct that represents the IGBT geometry.
% :param struct EC2: struct that represents the Capacitor geometry.
% :param struct bauraum: struct that represents the bauraum geometry.
% :param double cable_radius: Radius of the prescribed cable
% :param double radius_of_curvature: Minimum allowable radius of curvature 
%                                    of the prescribed cable 
%
% :return:
%   *[length,max_curvature, curvature, collisionwith1, collisionwith2, Penalty_spline, fn]*
%      - **length** : Length of the spline.
%      - **max_curvature** : Max curvature of the spline.
%      - **curvature** : Array with curvature values of segments of the
%        spline.
%      - **collisionwith1** : Percentage of collision between spline and
%        first EC.
%      - **collisionwith2** : Percentage of collision between spline and
%        second EC.
%      - **Penalty_spline,** : Penalty to be applied to the fitness of that
%        individual. Can be either 0 or 1000 depending if the spline was
%        successfully created
%      - **fn** : Function that describes the spline so it can be plotted
%        with the :func:`Visualize`.
%
%
% **Example in Code**
%
% .. code-block::
%
%    [splinelength,max_curvature, curvature, collisionwith1, collisionwith2, Penalty_spline, fn] = spline_connection (EC1_SurfCenter(1,:), pt2, pt3, EC2_SurfCenter(2,:),EC1_center, EC2_center, EC1, EC2, bauraum, cable_radius, radius_of_curvature);
%



% initial set values for collision
collisionwith1 = 1;
collisionwith2 = 1;
collisionwithbauraum = 0;

% determining intermediate control points
vector = pt3-pt2;
pt2a = pt2 + 0.25 .* vector ;
pt2b = pt2 + 0.5 .* vector ;
pt2c = pt2 + 0.75 .* vector ;

% control points vector
xyz = [pt1; pt2; pt2a; pt2b; pt2c; pt3; pt4];

ax = gca;
ax.XTick = [];
ax.YTick = [];
ax.ZTick = [];

% initial iteration value
iteration = 0;
% set initial maximum curvature
max_curvature = 0.5;
% Enter iteration loop in case of collisions or curvature violation
while (collisionwith1>0 || collisionwith2>0 || collisionwithbauraum<1 || max_curvature<(1/radius_of_curvature) || iteration < 30)

% iteration update
iteration = iteration+1;
% vector direction to shift control point 2a in direction parallel to line
% joining point and the component 1 centre
vector1 = pt2a-EC1_center;
 vector1 = vector1 / norm(vector1);
% vector direction to shift control point 2c in direction parallel to line
% joining point and the component 2 centre
vector2 = pt2c-EC2_center;
 vector2 = vector2 / norm(vector2);
% vector direction to shift control point 2b in direction perpendicular to
% line joining the centres of the components
 midptvector = EC2_center - EC1_center;
ptcentervector = pt2b - EC1_center;
C = dot(ptcentervector, midptvector);
C = (C / (norm(midptvector))^2)*midptvector;
ptonline = EC1_center + C;
vector3 = pt2b - ptonline;
vector3 = vector3 / norm(vector3);

% Update the control points locations
pt2a = pt2a + 0.5*vector1 ;
pt2b = pt2b + 0.5*vector3  ;
pt2c = pt2c + 0.5*vector2 ;
% updated control points vector
xyz = [pt1; pt2; pt2a; pt2b; pt2c; pt3; pt4];

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

% Vectors to store collision and curvature values
curvature = [];
collision1 = [];
collision2 = [];
collisionbauraum = [];


% Discretize the spline to check for collision and curvature
    for i=1:Spline.pieces
            %Obtaining the cartesian coordinates from parametric form
             x = Spline.coefs(i*Spline.dim-2,:);
             y = Spline.coefs(i*Spline.dim-1,:);
             z = Spline.coefs(i*Spline.dim,:);
            %Obtaining the first derivative from parametric form
            D_x = D_Spline.coefs(i*Spline.dim-2,:);
            D_y = D_Spline.coefs(i*Spline.dim-1,:);
            D_z = D_Spline.coefs(i*Spline.dim,:);
            %Obtaining the second derivative from parametric form
            DD_x = DD_Spline.coefs(i*Spline.dim-2,:);
            DD_y = DD_Spline.coefs(i*Spline.dim-1,:);
            DD_z = DD_Spline.coefs(i*Spline.dim,:);
            
            % discretization
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
                % curvature calculation               
                k = sqrt(((DD_z_t)*(D_y_t) - (DD_y_t)*(D_z_t))^2 +((DD_x_t)*(D_z_t) - (DD_z_t)*(D_x_t))^2 +((DD_y_t)*(D_x_t) - (DD_x_t)*(D_y_t))^2 )/((D_x_t)^2 + (D_y_t)^2 +(D_z_t)^2)^1.5;
                curvature(end+1) = k; 
                
                % collision check at each point
                collision1(end+1) = inpolyhedron(EC1, [x_t y_t z_t]);
                collision2(end+1) = inpolyhedron(EC2, [x_t y_t z_t]);
                collisionbauraum(end+1) = inpolyhedron(bauraum, [x_t y_t z_t]);
            end

    end
% saving curvature values
max_curvature = max(curvature); 
max_curvature_index = find(curvature == max_curvature);

% total collision
collisionwith1 = sum(collision1)/size(collision1,2);
collisionwith2 = sum(collision2)/size(collision2,2);
collisionwithbauraum = sum(collisionbauraum)/size(collisionbauraum,2);

% Loop exit in case of no collisions and allowable curvature or maximum iterations
    if (collisionwith1<0.01 && collisionwith2<0.01 && collisionwithbauraum>0.99 && max_curvature>(1/radius_of_curvature) || iteration == 30)
        Penalty_spline = 0;
        if iteration == 10
            Penalty_spline = 1000;
        end
        % Accomodating cable radius
        pt2a = pt2a + cable_radius*vector1 ;
        pt2b = pt2b + cable_radius*vector3  ;
        pt2c = pt2c + cable_radius*vector2 ;
        % final spline control points
        xyz = [pt1; pt2; pt2a; pt2b; pt2c; pt3; pt4];
        fn = cscvn(xyz');
        
        break;
    end


end
% final spline control points
xyz = [pt1; pt2; pt2a; pt2b; pt2c; pt3; pt4];

ax = gca;
ax.XTick = [];
ax.YTick = [];
ax.ZTick = [];

hold on  
end


