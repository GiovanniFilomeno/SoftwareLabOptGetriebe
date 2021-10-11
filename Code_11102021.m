clear all
clc
tStart = cputime;

%% Loading all the solids
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


% Defining the limits of the bauraum to generate the population inside.
xmax = max(bauraum.vertices(:,1));
xmin = min(bauraum.vertices(:,1));
ymax = max(bauraum.vertices(:,2));
ymin = min(bauraum.vertices(:,2));
zmax = max(bauraum.vertices(:,3));
zmin = min(bauraum.vertices(:,3));

%Moving the bauraum to the origin
bauraum_mov = [xmin ymin zmin];
bauraum.vertices = bauraum.vertices - bauraum_mov;

xmax = xmax-xmin;
ymax = ymax-ymin;
zmax = zmax-zmin;
xmin=0;
ymin=0;
zmin=0;

%Read 4 Electric components and set Right dimensions
Cube=stlread("cube.STL");
V = Cube.Points;

V2(:,2) = 0.566.*V(:,1);
V2(:,1) = 0.2045.*V(:,2);
V2(:,3) = 0.07629.*V(:,3);

V3(:,2) = 1.12736.*V(:,1);
V3(:,1) = 0.24335.*V(:,2);
V3(:,3) = 0.12534.*V(:,3);

V4(:,2) = 0.3302.*V(:,1);
V4(:,1) = 0.06135.*V(:,2);
V4(:,3) = 0.06812.*V(:,3);

V5(:,2) = 1.4151.*V(:,1);
V5(:,1) = 0.3681.*V(:,2);
V5(:,3) = 0.0545.*V(:,3);

F = Cube.ConnectivityList;
N = faceNormal(Cube);

%IGBT
EC1 = struct('faces',F,'vertices',V2);

%DC-Link Capacitor
EC2 = struct('faces',F,'vertices',V3);

%EMI-Filter
EC3 = struct('faces',F,'vertices',V4);

%Control board
EC4 = struct('faces',F,'vertices',V5);



%Plot Bauraum
BRaum = patch(bauraum,'FaceColor',       [0 0 1.0], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15);
BRaum.FaceVertexAlphaData = 0.2;    % Set constant transparency 
BRaum.FaceAlpha = 'flat' ;          % Interpolate to find face transparency
hold on
view([45 90 135]);
xlabel('X') 
ylabel('Y')
zlabel('Z')
%axis([0 10 0 10 0 10])
camlight('headlight');
material('default');

%% Positioning(manually) of the gearbox inside the bauraum
[EM_start, EM_vector, Anchor, EM_radius] = FindingAxis("Gearbox1.stl");
EM_end = EM_start + EM_vector;
zero = [0 0 0];

%Definition(manually of the axis of the bauraum shaft 
point = [2539 14 -40] - bauraum_mov;
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
shaftx = [point(1) point(1)];
shafty = [point(2) point(2) - 119];
shaftz = [point(3) point(3)];
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


%% Defining the limits and centers of the Electric components

%IGBT
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

%CENTER POINTS OF SURFACE OF COMPONENT 1
EC1_center1 = [EC1_xmax (EC1_ymax + EC1_ymin)/2 (EC1_zmax + EC1_zmin)/2]; %Counterclockwise

EC1_center2 = [EC1_xmin (EC1_ymax + EC1_ymin)/2 (EC1_zmax + EC1_zmin)/2]; %Clockwise

%To use with Control board
EC1_center3 = [(EC1_xmax + EC1_xmin)/2 (EC1_ymax + EC1_ymin)/2 EC1_zmax];  


%Capacitor
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

%CENTER POINTS OF SURFACE OF COMPONENT 2
EC2_center1 = [EC2_xmax (EC2_ymax + EC2_ymin)/2  (EC2_zmax + EC2_zmin)/2];

EC2_center2 = [EC2_xmin (EC2_ymax + EC2_ymin)/2 (EC2_zmax + EC2_zmin)/2];


%EMI-filter
EC3_xmax = max(EC3.vertices(:,1));
EC3_xmin = min(EC3.vertices(:,1));
EC3_ymax = max(EC3.vertices(:,2));
EC3_ymin = min(EC3.vertices(:,2));
EC3_zmax = max(EC3.vertices(:,3));
EC3_zmin = min(EC3.vertices(:,3));
EC3_length = abs(EC3_xmax - EC3_xmin);
EC3_width = abs(EC3_ymax - EC3_ymin);
EC3_height = abs(EC3_zmax - EC3_zmin);
EC3_center = [(EC3_xmax + EC3_xmin)/2 (EC3_ymax + EC3_ymin)/2 (EC3_zmax + EC3_zmin)/2];

%CENTER POINTS OF SURFACE OF COMPONENT 3
EC3_center1 = [EC3_xmax (EC3_ymax + EC3_xmin)/2 (EC3_zmax + EC3_zmin)/2]; 

EC3_center2 = [EC3_xmin (EC3_ymax + EC3_xmin)/2 (EC3_zmax + EC3_zmin)/2];

%Control board
EC4_xmax = max(EC4.vertices(:,1));
EC4_xmin = min(EC4.vertices(:,1));
EC4_ymax = max(EC4.vertices(:,2));
EC4_ymin = min(EC4.vertices(:,2));
EC4_zmax = max(EC4.vertices(:,3));
EC4_zmin = min(EC4.vertices(:,3));
EC4_length = abs(EC4_xmax - EC4_xmin);
EC4_width = abs(EC4_ymax - EC4_ymin);
EC4_height = abs(EC4_zmax - EC4_zmin);
EC4_center = [(EC4_xmax + EC4_xmin)/2 (EC4_ymax + EC4_ymin)/2 (EC4_zmax + EC4_zmin)/2]; 

%To use with IGBT
EC4_center1 = [(EC4_xmax + EC4_xmin)/2 (EC4_ymax + EC4_ymin)/2 EC4_zmin]; 
     
     
%% Genetic Algorithm implementation
pop = 500;
iter = 20;


%Create first sample population and evaluate fitness

%Orientation of the gearbox
gb_rot = rad_range.*rand(pop,1);
gb_mov = zeros(pop,1);
for i=1:pop
    if i==1
        gb_mov(i) = gb_rot(i);
    else
        gb_mov(i) = gb_rot(i) - gb_rot(i-1);
    end
end

%IGBT
r1 = EM_radius + (EC1_height/2) + (EC1_height).*rand(pop,1);
y1= EM_start(2) + (EM_vector(2)).*rand(pop,1);
theta1 = -2*pi + 2*pi.*rand(pop,1);
EC1_rot = -pi/12 + (pi/6).*rand(pop,1);
theta1_mov = zeros(pop,1);
EC1_mov = zeros(pop,1);
for i=1:pop
    if i==1
        EC1_mov(i) = theta1(i) + EC1_rot(i);
    else
        EC1_mov(i) = theta1(i) + EC1_rot(i) - theta1(i-1) - EC1_rot(i-1);
    end
end


%Condensator
r2 = EM_radius + (EC2_height/2) + (EC2_height).*rand(pop,1);
y2= EM_start(2) + (EM_vector(2)).*rand(pop,1);
angle_range2 = (EC1_length + EC2_length)/(2*(EM_radius)); 
theta2_direction = rand(pop,1);
EC2_placement = round(theta2_direction);
theta2 = zeros(pop,1);
for i=1:pop
    if theta2_direction(i) > 0.5
        theta2(i) = theta1(i) + angle_range2 + (pi/15)*rand(); %Component counterclockwise from 1
    else
        theta2(i) = theta1(i) - angle_range2 - (pi/15)*rand(); %Component clockwise from 1
    end
end   

EC2_rot = -pi/12 + (pi/6).*rand(pop,1);
EC2_mov = zeros(pop,1);
for i=1:pop
    if i==1
        EC2_mov(i) = theta2(i) + EC2_rot(i);
    else
        EC2_mov(i) = theta2(i) + EC2_rot(i) - theta2(i-1) - EC2_rot(i-1);
    end
end


%EMI-Filter
r3 = EM_radius + (EC3_height/2) + (EC3_height).*rand(pop,1);
y3= EM_start(2) + (EM_vector(2)).*rand(pop,1);
angle_range3 = (EC2_length + EC3_length)/(2*(EM_radius)); 
theta3 = zeros(pop,1);
for i=1:pop
    if theta2(i) > theta1(i)
        theta3(i) = theta2(i) + angle_range3 + (pi/15)*rand();
    else
        theta3(i) = theta2(i) - angle_range3 - (pi/15)*rand();
    end
end   
EC3_rot = -pi/12 + (pi/6).*rand(pop,1);
EC3_mov = zeros(pop,1);
for i=1:pop
    if i==1
        EC3_mov(i) = theta3(i) + EC3_rot(i);
    else
        EC3_mov(i) = theta3(i) + EC3_rot(i) - theta3(i-1) - EC3_rot(i-1);
    end
end



%Control board
r4 = r1 + (EC1_height) + (EC1_height).*rand(pop,1);
y4= EM_start(2) + (EM_vector(2)).*rand(pop,1);
theta4 = theta1;
% EC4_rot = -pi/12 + (pi/6).*rand(pop,1);
EC4_rot = EC1_rot;
theta4_mov = zeros(pop,1);
EC4_mov = zeros(pop,1);
for i=1:pop
    if i==1
        EC4_mov(i) = theta4(i) + EC4_rot(i);
    else
        EC4_mov(i) = theta4(i) + EC4_rot(i) - theta4(i-1) - EC4_rot(i-1);
    end
end



%Array that contain the fitness of the chromosome
d = zeros(pop,1);


for i=1:pop
    
    %Move the components and associated variables to the position of that chromosome
    %Gearbox
    fv.vertices = rotation(fv.vertices, Anchor, 2, gb_mov(i,1));
    EM_start = rotation(EM_start, Anchor, 2, gb_mov(i,1));
    EM_end = EM_start + EM_vector;
    
    coords1(i,:) = [EM_start(1)+ r1(i).*sin(theta1(i)) y1(i) EM_start(3) + r1(i).*cos(theta1(i))];
    coords2(i,:) = [EM_start(1)+ r2(i).*sin(theta2(i)) y2(i) EM_start(3) + r2(i).*cos(theta2(i))];
    coords3(i,:) = [EM_start(1)+ r3(i).*sin(theta3(i)) y3(i) EM_start(3) + r3(i).*cos(theta3(i))];
    coords4(i,:) = [EM_start(1)+ r4(i).*sin(theta4(i)) y4(i) EM_start(3) + r4(i).*cos(theta4(i))];
    
    %IGBT
    EC1.vertices = EC1.vertices - EC1_center + coords1(i,:);
    EC1_center1 = EC1_center1 - EC1_center + coords1(i,:);
    EC1_center2 = EC1_center2 - EC1_center + coords1(i,:);
    EC1_center3 = EC1_center3 - EC1_center + coords1(i,:);
    EC1_center = coords1(i,:);
    EC1.vertices = rotation(EC1.vertices, EC1_center, 2, -EC1_mov(i,1));
    EC1_center1 = rotation(EC1_center1, EC1_center, 2, -EC1_mov(i,1));
    EC1_center2 = rotation(EC1_center2, EC1_center, 2, -EC1_mov(i,1));
    EC1_center3 = rotation(EC1_center3, EC1_center, 2, -EC1_mov(i,1));
    
    %Capacitor
    EC2.vertices = EC2.vertices - EC2_center + coords2(i,:);
    EC2_center1 = EC2_center1 - EC2_center + coords2(i,:);
    EC2_center2 = EC2_center2 - EC2_center + coords2(i,:);
    EC2_center = coords2(i,:);
    EC2.vertices = rotation(EC2.vertices, EC2_center, 2, -EC2_mov(i,1));
    EC2_center1 = rotation(EC2_center1, EC2_center, 2, -EC2_mov(i,1));
    EC2_center2 = rotation(EC2_center2, EC2_center, 2, -EC2_mov(i,1));
    
    %EMI-filter
    EC3.vertices = EC3.vertices - EC3_center + coords3(i,:);
    EC3_center1 = EC3_center1 - EC3_center + coords3(i,:);
    EC3_center = coords3(i,:);
    EC3.vertices = rotation(EC3.vertices, EC3_center, 2, -EC3_mov(i,1));
    EC3_center1 = rotation(EC3_center1, EC3_center, 2, -EC3_mov(i,1));
    
    %Control board
    EC4.vertices = EC4.vertices - EC4_center + coords4(i,:);
    EC4_center1 = EC4_center1 - EC4_center + coords4(i,:);
    EC4_center = coords4(i,:);
    EC4.vertices = rotation(EC4.vertices, EC4_center, 2, -EC4_mov(i,1));
    EC4_center1 = rotation(EC4_center1, EC4_center, 2, -EC4_mov(i,1));
    

    %Checking for collisions
    
    %Components inside the bauraum
    IN = inpolyhedron(bauraum,EC1.vertices);
    in1_bauraum = sum(IN)/length(IN);
    IN = inpolyhedron(bauraum,EC2.vertices);
    in2_bauraum = sum(IN)/length(IN);
    IN = inpolyhedron(bauraum,EC3.vertices);
    in3_bauraum = sum(IN)/length(IN);
    IN = inpolyhedron(bauraum,EC4.vertices);
    in4_bauraum = sum(IN)/length(IN);
    
    %Components outside the gearbox/EM
    IN = inpolyhedron(fv,EC1.vertices);
    out1_gearbox = sum(IN)/length(IN);
    IN = inpolyhedron(fv,EC2.vertices);
    out2_gearbox = sum(IN)/length(IN);
    IN = inpolyhedron(fv,EC3.vertices);
    out3_gearbox = sum(IN)/length(IN);
    IN = inpolyhedron(fv,EC4.vertices);
    out4_gearbox = sum(IN)/length(IN);
    
    
    %Components with themselves
    IN = inpolyhedron(EC1,EC2.vertices);
    collision12 = sum(IN)/length(IN);
    IN = inpolyhedron(EC1,EC3.vertices);
    collision13 = sum(IN)/length(IN);
    IN = inpolyhedron(EC1,EC4.vertices);
    collision14 = sum(IN)/length(IN);
    IN = inpolyhedron(EC2,EC3.vertices);
    collision23 = sum(IN)/length(IN);
    IN = inpolyhedron(EC2,EC4.vertices);
    collision24 = sum(IN)/length(IN);
    IN = inpolyhedron(EC3,EC4.vertices);
    collision34 = sum(IN)/length(IN);
    

    d(i,1) = optimization_function(EC1_center, EC1_center1, EC1_center2, EC1_center3, EC2_center1,EC2_center2, EC3_center1, EC3_center2, EC4_center1, EM_start, EM_end, in1_bauraum, in2_bauraum, in3_bauraum, in4_bauraum, out1_gearbox,  out2_gearbox, out3_gearbox, out4_gearbox, collision12, collision13, collision14, collision23, collision24, collision34, EC2_placement(i));
    
%     hold on     
%      
%     figure(i)
%     
%     BRaum = patch(bauraum,'FaceColor',       [0 0 1.0], ...
%          'EdgeColor',       'none',        ...
%          'FaceLighting',    'gouraud',     ...
%          'AmbientStrength', 0.15);
%     BRaum.FaceVertexAlphaData = 0.2;    % Set constant transparency 
%     BRaum.FaceAlpha = 'flat' ;          % Interpolate to find face transparency
%     hold on
%     view([45 90 135]);
%     xlabel('X') 
%     ylabel('Y')
%     zlabel('Z')
%     %axis([0 10 0 10 0 10])
%     camlight('headlight');
%     material('default');
%     
%     patch(fv,'FaceColor',       [0.8 0.8 1], ...
%          'FaceLighting',    'gouraud',     ...
%          'AmbientStrength', 0.15);
%     plot3(shaftx', shafty', shaftz','LineWidth',2,'Color','red')
% 
% 
%     patch(EC1,'FaceColor',       [0.8 0 0], ...
%              'EdgeColor',       'green',        ...
%              'FaceLighting',    'gouraud',     ...
%              'AmbientStrength', 0.15);
%     patch(EC2,'FaceColor',       [0 0.8 0], ...
%          'EdgeColor',       'red',        ...
%          'FaceLighting',    'gouraud',     ...
%          'AmbientStrength', 0.15);
%      patch(EC3,'FaceColor',       [1  1 0], ...
%              'EdgeColor',       'black',        ...
%              'FaceLighting',    'gouraud',     ...
%              'AmbientStrength', 0.15);
%     patch(EC4,'FaceColor',       [1 0 1], ...
%          'EdgeColor',       'yellow',        ...
%          'FaceLighting',    'gouraud',     ...
%          'AmbientStrength', 0.15);
%     
%     
end

%Return the gearbox to initial angle
fv.vertices = rotation(fv.vertices, Anchor, 2, -gb_rot(pop,1));
EM_start = rotation(EM_start, Anchor, 2, -gb_rot(pop,1));
EM_end = EM_start + EM_vector;

%Return components to initial angle
EC1.vertices = rotation(EC1.vertices, EC1_center, 2, theta1(pop,1)+EC1_rot(pop,1));
EC1_center1 = rotation(EC1_center1, EC1_center, 2, EC1_rot(pop,1));
EC1_center2 = rotation(EC1_center2, EC1_center, 2, EC1_rot(pop,1));
EC1_center3 = rotation(EC1_center3, EC1_center, 2, EC1_rot(pop,1));

EC2.vertices = rotation(EC2.vertices, EC2_center, 2, theta2(pop,1)+EC2_rot(pop,1));
EC2_center1 = rotation(EC2_center1, EC2_center, 2, EC2_rot(pop,1));
EC2_center2 = rotation(EC2_center2, EC2_center, 2, EC2_rot(pop,1));

EC3.vertices = rotation(EC3.vertices, EC3_center, 2, theta3(pop,1)+EC3_rot(pop,1));
EC3_center1 = rotation(EC3_center1, EC3_center, 2, EC3_rot(pop,1));

EC4.vertices = rotation(EC4.vertices, EC4_center, 2, theta4(pop,1)+EC4_rot(pop,1));
EC4_center1 = rotation(EC4_center1, EC4_center, 2, EC4_rot(pop,1));

%Sort chromosomes based on their fitness values
[srt,I]=sort(d);


fprintf('Fitness of Generation 1: %d \n',srt(1,1));
pause(1)

%Select the best half of chromosomes from this generation
bestgb = zeros(pop/2,1);
best1 = zeros(pop/2,3);
best2 = zeros(pop/2,3);
best3 = zeros(pop/2,3);
best4 = zeros(pop/2,3);
best_rot1 = zeros(pop/2,1);
best_rot2 = zeros(pop/2,1);
best_rot3 = zeros(pop/2,1);
best_rot4 = zeros(pop/2,1);
bestcoords1 = zeros(pop/2,3);
bestcoords2 = zeros(pop/2,3);
bestcoords3 = zeros(pop/2,3);
bestcoords4 = zeros(pop/2,3);
bestd = zeros(pop/2,1);
best_placement = zeros(pop,2/1);
for i=1:pop/2
    bestgb(i,1) = gb_rot(I(i,1),1);
    best1(i,1) = r1(I(i,1),:);
    best2(i,1) = r2(I(i,1),:);
    best3(i,1) = r3(I(i,1),:);
    best4(i,1) = r4(I(i,1),:);
    best1(i,2) = y1(I(i,1),:);
    best2(i,2) = y2(I(i,1),:);
    best3(i,2) = y3(I(i,1),:);
    best4(i,2) = y4(I(i,1),:);
    best1(i,3) = theta1(I(i,1),:);
    best2(i,3) = theta2(I(i,1),:);
    best3(i,3) = theta3(I(i,1),:);
    best4(i,3) = theta4(I(i,1),:);
    bestcoords1(i,:) = coords1(I(i,1),:);
    bestcoords2(i,:) = coords2(I(i,1),:);
    bestcoords3(i,:) = coords3(I(i,1),:);
    bestcoords4(i,:) = coords4(I(i,1),:);
    best_rot1(i,1) = EC1_rot(I(i,1),1);
    best_rot2(i,1) = EC2_rot(I(i,1),1);
    best_rot3(i,1) = EC3_rot(I(i,1),1);
    best_rot4(i,1) = EC4_rot(I(i,1),1);
    bestd(i,1) = d(I(i,1),1);
    best_placement(i,1) = EC2_placement(I(i,1),1);
end



 for g = 2:iter
      for k = 1: pop/2
          
          
        %Selection
        sel_coeff = d(I(pop/2,1),1);
        parents1 = zeros(2,3);
        parents2 = zeros(2,4);
        parents3 = zeros(2,3);
        parents4 = zeros(2,3);
        parents_rot = zeros(2,5);
        for j=1:2
            rand_sel = 1 + round(pop.*rand());
            for i = rand_sel:pop
                if d(i,1)<sel_coeff
                    parents1(j,1) = r1(i,1);
                    parents1(j,2) = y1(i,1);
                    parents1(j,3) = theta1(i,1);
                    parents2(j,1) = r2(i,1);
                    parents2(j,2) = y2(i,1);
                    parents2(j,3) = theta2(i,1);
                    parents2(j,4) = EC2_placement(i,1);
                    parents3(j,1) = r3(i,1);
                    parents3(j,2) = y3(i,1);
                    parents3(j,3) = theta3(i,1);
                    parents4(j,1) = r4(i,1);
                    parents4(j,2) = y4(i,1);
                    parents4(j,3) = theta4(i,1);
                    parents_rot(j,1) = gb_rot(i,1);
                    parents_rot(j,2) = EC1_rot(i,1);
                    parents_rot(j,3) = EC2_rot(i,1);
                    parents_rot(j,4) = EC3_rot(i,1);
                    parents_rot(j,5) = EC4_rot(i,1);
                    
                    break
                end
            end
        end
 
        %Crossover
        child1 = zeros(2,3);
        child2 = zeros(2,4);
        child3 = zeros(2,3);
        child4 = zeros(2,3);
        child_rot = zeros(2,5);
        alpha = 0.15;
        beta = 0.85;
        child1(1,1) = alpha.*parents1(1,1) + (1-alpha).*parents1(2,1);
        child1(2,1) = beta.*parents1(1,1) + (1-beta).*parents1(2,1);
        child2(1,1) = alpha.*parents2(1,1) + (1-alpha).*parents2(2,1);
        child2(2,1) = beta.*parents2(1,1) + (1-beta).*parents2(2,1);
        child3(1,1) = alpha.*parents3(1,1) + (1-alpha).*parents3(2,1);
        child3(2,1) = beta.*parents3(1,1) + (1-beta).*parents3(2,1);
        child4(1,1) = alpha.*parents4(1,1) + (1-alpha).*parents4(2,1);
        child4(2,1) = beta.*parents4(1,1) + (1-beta).*parents4(2,1);
        
        child1(1,2) = alpha.*parents1(1,2) + (1-alpha).*parents1(2,2);
        child1(2,2) = beta.*parents1(1,2) + (1-beta).*parents1(2,2);
        child2(1,2) = alpha.*parents2(1,2) + (1-alpha).*parents2(2,2);
        child2(2,2) = beta.*parents2(1,2) + (1-beta).*parents2(2,2);
        child3(1,2) = alpha.*parents3(1,2) + (1-alpha).*parents3(2,2);
        child3(2,2) = beta.*parents3(1,2) + (1-beta).*parents3(2,2);
        child4(1,2) = alpha.*parents4(1,2) + (1-alpha).*parents4(2,2);
        child4(2,2) = beta.*parents4(1,2) + (1-beta).*parents4(2,2);
        
        child1(1,3) = parents1(1,3);
        child1(2,3) = parents1(2,3);
        child2(1,3) = parents2(1,3);
        child2(2,3) = parents2(2,3);
        child3(1,3) = parents3(1,3);
        child3(2,3) = parents3(2,3);
        child4(1,3) = parents4(1,3);
        child4(2,3) = parents4(2,3);
        
        child2(1,4) = parents2(1,4);
        child2(2,4) = parents2(1,4);
        
        child_rot(1,:) = alpha*parents_rot(1,:) + (1-alpha)*parents_rot(2,:);
        child_rot(2,:) = beta*parents_rot(1,:) + (1-beta)*parents_rot(2,:);
        


%       Mutation
        for i = 1:2
            mut_prob = rand();
            if (mut_prob>=0.75)
                child1(i,2) = child1(i,2) + 5;
                child2(i,2) = child2(i,2) + 5;
                child3(i,2) = child3(i,2) + 5;
                child4(i,2) = child4(i,2) + 5;
            elseif (mut_prob>=0.5) && (mut_prob<0.75)
                child1(i,2) = child1(i,2) - 5;
                child2(i,2) = child2(i,2) - 5;
                child3(i,2) = child3(i,2) - 5;
                child4(i,2) = child4(i,2) - 5;
            end
            mut_prob = rand();
            if (mut_prob>=0.75)
                child1(i,3) = child1(i,3) + (pi/9)*rand();
                child2(i,3) = child2(i,3) + (pi/9)*rand();
                child3(i,3) = child3(i,3) + (pi/9)*rand();
                child4(i,3) = child1(i,3);
            elseif (mut_prob>=0.5) && (mut_prob<0.75)
                child1(i,3) = child1(i,3) - (pi/9)*rand();
                child2(i,3) = child2(i,3) - (pi/9)*rand();
                child3(i,3) = child3(i,3) - (pi/9)*rand();
                child4(i,3) = child1(i,3);
            end
            for j=1:5
                mut_prob = rand();
                if mut_prob > 0.5 
                    child_rot(i,j) = child_rot(i,j) + (pi/36);
                else
                    child_rot(i,j) = child_rot(i,j) - (pi/36);
                end
                 
            end
        end
        
        %Write first half of new generation 
        r1(k,1) = child1(1,1);
        y1(k,1) = child1(1,2);
        theta1(k,1) = child1(1,3);
        r1(k+1,1) = child1(2,1);
        y1(k+1,1) = child1(2,2);
        theta1(k+1,1) = child1(2,3);
        r2(k,1) = child2(1,1);
        y2(k,1) = child2(1,2);
        theta2(k,1) = child2(1,3);
        EC2_placement(k,1) = child2(1,4);
        r2(k+1,1) = child2(2,1);
        y2(k+1,1) = child2(2,2);
        theta2(k+1,1) = child2(2,3);
        EC2_placement(k+1,1) = child2(2,4);
        r3(k,1) = child3(1,1);
        y3(k,1) = child3(1,2);
        theta3(k,1) = child3(1,3);
        r3(k+1,1) = child3(2,1);
        y3(k+1,1) = child3(2,2);
        theta3(k+1,1) = child3(2,3);
        r4(k,1) = child4(1,1);
        y4(k,1) = child4(1,2);
        theta4(k,1) = child4(1,3);
        r4(k+1,1) = child4(2,1);
        y4(k+1,1) = child4(2,2);
        theta4(k+1,1) = child4(2,3);
        gb_rot(k,1) = child_rot(1,1);
        gb_rot(k+1,1) = child_rot(2,1);
        EC1_rot(k,1) = child_rot(1,2);
        EC1_rot(k+1,1) = child_rot(2,2);
        EC2_rot(k,1) = child_rot(1,3);
        EC2_rot(k+1,1) = child_rot(2,3);
        EC3_rot(k,1) = child_rot(1,4);
        EC3_rot(k+1,1) = child_rot(2,4);
        EC4_rot(k,1) = child_rot(1,5);
        EC4_rot(k+1,1) = child_rot(2,5);
      end
      
%Write other half of new population from the best of previous generation     
      for i=(pop/2+1):pop
          r1(i,1) = best1(i-(pop/2),1);
          r2(i,1) = best2(i-(pop/2),1);
          r3(i,1) = best3(i-(pop/2),1);
          r4(i,1) = best4(i-(pop/2),1);
          y1(i,1) = best1(i-(pop/2),2);
          y2(i,1) = best2(i-(pop/2),2);
          y3(i,1) = best3(i-(pop/2),2);
          y4(i,1) = best4(i-(pop/2),2);
          theta1(i,1) = best1(i-(pop/2),3);
          theta2(i,1) = best2(i-(pop/2),3);
          theta3(i,1) = best3(i-(pop/2),3);
          theta4(i,1) = best4(i-(pop/2),3);
          EC2_placement(i,1) = best_placement(i-(pop/2),1); 
          gb_rot(i,1) = bestgb(i-(pop/2),1);
          EC1_rot(i,1) = best_rot1(i-(pop/2),1);
          EC2_rot(i,1) = best_rot2(i-(pop/2),1);
          EC3_rot(i,1) = best_rot3(i-(pop/2),1);
          EC4_rot(i,1) = best_rot4(i-(pop/2),1);
          coords1(i,:) = bestcoords1(i-(pop/2),1);
          coords2(i,:) = bestcoords2(i-(pop/2),1);
          coords3(i,:) = bestcoords3(i-(pop/2),1);
          coords4(i,:) = bestcoords4(i-(pop/2),1);
          d(i,1) = bestd(i-(pop/2),1);
      end
     
      for i=1:pop
        if i==1
            gb_mov(i) = gb_rot(i);
        else
            gb_mov(i) = gb_rot(i) - gb_rot(i-1);
        end
      end
      
      for i=1:pop
        if i==1
            EC1_mov(i) = theta1(i) + EC1_rot(i);
        else
            EC1_mov(i) = theta1(i) + EC1_rot(i) - theta1(i-1) - EC1_rot(i-1);
        end
      end
      
      for i=1:pop
        if i==1
            EC2_mov(i) = theta2(i) + EC2_rot(i);
        else
            EC2_mov(i) = theta2(i) + EC2_rot(i) - theta2(i-1) - EC2_rot(i-1);
        end
      end
      
      for i=1:pop
        if i==1
            EC3_mov(i) = theta3(i) + EC3_rot(i);
        else
            EC3_mov(i) = theta3(i) + EC3_rot(i) - theta3(i-1) - EC3_rot(i-1);
        end
      end
      
      for i=1:pop
        if i==1
            EC4_mov(i) = theta4(i) + EC4_rot(i);
        else
            EC4_mov(i) = theta4(i) + EC4_rot(i) - theta4(i-1) - EC4_rot(i-1);
        end
      end
      
%Test new generation
        for i=1:pop

            %Move the components and associated variables to the position of that chromosome
            %Gearbox
            fv.vertices = rotation(fv.vertices, Anchor, 2, gb_mov(i,1));
            EM_start = rotation(EM_start, Anchor, 2, gb_mov(i,1));
            EM_end = EM_start + EM_vector;

            coords1(i,:) = [EM_start(1)+ r1(i).*sin(theta1(i)) y1(i) EM_start(3) + r1(i).*cos(theta1(i))];
            coords2(i,:) = [EM_start(1)+ r2(i).*sin(theta2(i)) y2(i) EM_start(3) + r2(i).*cos(theta2(i))];
            coords3(i,:) = [EM_start(1)+ r3(i).*sin(theta3(i)) y3(i) EM_start(3) + r3(i).*cos(theta3(i))];
            coords4(i,:) = [EM_start(1)+ r4(i).*sin(theta4(i)) y4(i) EM_start(3) + r4(i).*cos(theta4(i))];

            %IGBT
            EC1.vertices = EC1.vertices - EC1_center + coords1(i,:);
            EC1_center1 = EC1_center1 - EC1_center + coords1(i,:);
            EC1_center2 = EC1_center2 - EC1_center + coords1(i,:);
            EC1_center3 = EC1_center3 - EC1_center + coords1(i,:);
            EC1_center = coords1(i,:);
            EC1.vertices = rotation(EC1.vertices, EC1_center, 2, -EC1_mov(i,1));
            EC1_center1 = rotation(EC1_center1, EC1_center, 2, -EC1_mov(i,1));
            EC1_center2 = rotation(EC1_center2, EC1_center, 2, -EC1_mov(i,1));
            EC1_center3 = rotation(EC1_center3, EC1_center, 2, -EC1_mov(i,1));

            %Capacitor
            EC2.vertices = EC2.vertices - EC2_center + coords2(i,:);
            EC2_center1 = EC2_center1 - EC2_center + coords2(i,:);
            EC2_center2 = EC2_center2 - EC2_center + coords2(i,:);
            EC2_center = coords2(i,:);
            EC2.vertices = rotation(EC2.vertices, EC2_center, 2, -EC2_mov(i,1));
            EC2_center1 = rotation(EC2_center1, EC2_center, 2, -EC2_mov(i,1));
            EC2_center2 = rotation(EC2_center2, EC2_center, 2, -EC2_mov(i,1));

            %EMI-filter
            EC3.vertices = EC3.vertices - EC3_center + coords3(i,:);
            EC3_center1 = EC3_center1 - EC3_center + coords3(i,:);
            EC3_center = coords3(i,:);
            EC3.vertices = rotation(EC3.vertices, EC3_center, 2, -EC3_mov(i,1));
            EC3_center1 = rotation(EC3_center1, EC3_center, 2, -EC3_mov(i,1));

            %Control board
            EC4.vertices = EC4.vertices - EC4_center + coords4(i,:);
            EC4_center1 = EC4_center1 - EC4_center + coords4(i,:);
            EC4_center = coords4(i,:);
            EC4.vertices = rotation(EC4.vertices, EC4_center, 2, -EC4_mov(i,1));
            EC4_center1 = rotation(EC4_center1, EC4_center, 2, -EC4_mov(i,1));


            %Checking for collisions

            %Components inside the bauraum
            IN = inpolyhedron(bauraum,EC1.vertices);
            in1_bauraum = sum(IN)/length(IN);
            IN = inpolyhedron(bauraum,EC2.vertices);
            in2_bauraum = sum(IN)/length(IN);
            IN = inpolyhedron(bauraum,EC3.vertices);
            in3_bauraum = sum(IN)/length(IN);
            IN = inpolyhedron(bauraum,EC4.vertices);
            in4_bauraum = sum(IN)/length(IN);

            %Components outside the gearbox/EM
            IN = inpolyhedron(fv,EC1.vertices);
            out1_gearbox = sum(IN)/length(IN);
            IN = inpolyhedron(fv,EC2.vertices);
            out2_gearbox = sum(IN)/length(IN);
            IN = inpolyhedron(fv,EC3.vertices);
            out3_gearbox = sum(IN)/length(IN);
            IN = inpolyhedron(fv,EC4.vertices);
            out4_gearbox = sum(IN)/length(IN);


            %Components with themselves
            IN = inpolyhedron(EC1,EC2.vertices);
            collision12 = sum(IN)/length(IN);
            IN = inpolyhedron(EC1,EC3.vertices);
            collision13 = sum(IN)/length(IN);
            IN = inpolyhedron(EC1,EC4.vertices);
            collision14 = sum(IN)/length(IN);
            IN = inpolyhedron(EC2,EC3.vertices);
            collision23 = sum(IN)/length(IN);
            IN = inpolyhedron(EC2,EC4.vertices);
            collision24 = sum(IN)/length(IN);
            IN = inpolyhedron(EC3,EC4.vertices);
            collision34 = sum(IN)/length(IN);


            d(i,1) = optimization_function(EC1_center, EC1_center1, EC1_center2, EC1_center3, EC2_center1,EC2_center2, EC3_center1, EC3_center2, EC4_center1, EM_start, EM_end, in1_bauraum, in2_bauraum, in3_bauraum, in4_bauraum, out1_gearbox,  out2_gearbox, out3_gearbox, out4_gearbox, collision12, collision13, collision14, collision23, collision24, collision34, EC2_placement(i));

            hold on     
%      
%             figure(i)
% 
%             BRaum = patch(bauraum,'FaceColor',       [0 0 1.0], ...
%                  'EdgeColor',       'none',        ...
%                  'FaceLighting',    'gouraud',     ...
%                  'AmbientStrength', 0.15);
%             BRaum.FaceVertexAlphaData = 0.2;    % Set constant transparency 
%             BRaum.FaceAlpha = 'flat' ;          % Interpolate to find face transparency
%             hold on
%             view([45 90 135]);
%             xlabel('X') 
%             ylabel('Y')
%             zlabel('Z')
%             %axis([0 10 0 10 0 10])
%             camlight('headlight');
%             material('default');
% 
%             patch(fv,'FaceColor',       [0.8 0.8 1], ...
%                  'FaceLighting',    'gouraud',     ...
%                  'AmbientStrength', 0.15);
%             plot3(shaftx', shafty', shaftz','LineWidth',2,'Color','red')
% 
% 
%             patch(EC1,'FaceColor',       [0.8 0 0], ...
%                      'EdgeColor',       'green',        ...
%                      'FaceLighting',    'gouraud',     ...
%                      'AmbientStrength', 0.15);
%             patch(EC2,'FaceColor',       [0 0.8 0], ...
%                  'EdgeColor',       'red',        ...
%                  'FaceLighting',    'gouraud',     ...
%                  'AmbientStrength', 0.15);
%              patch(EC3,'FaceColor',       [1  1 0], ...
%                      'EdgeColor',       'black',        ...
%                      'FaceLighting',    'gouraud',     ...
%                      'AmbientStrength', 0.15);
%             patch(EC4,'FaceColor',       [1 0 1], ...
%                  'EdgeColor',       'yellow',        ...
%                  'FaceLighting',    'gouraud',     ...
%                  'AmbientStrength', 0.15);
% 
% 
%             
        end
        
        
        %Return the gearbox to initial angle
        fv.vertices = rotation(fv.vertices, Anchor, 2, -gb_rot(pop,1));
        EM_start = rotation(EM_start, Anchor, 2, -gb_rot(pop,1));
        EM_end = EM_start + EM_vector;

         %Return components to initial angle
        EC1.vertices = rotation(EC1.vertices, EC1_center, 2, theta1(pop,1)+EC1_rot(pop,1));
        EC1_center1 = rotation(EC1_center1, EC1_center, 2, EC1_rot(pop,1));
        EC1_center2 = rotation(EC1_center2, EC1_center, 2, EC1_rot(pop,1));
        EC1_center3 = rotation(EC1_center3, EC1_center, 2, EC1_rot(pop,1));

        EC2.vertices = rotation(EC2.vertices, EC2_center, 2, theta2(pop,1)+EC2_rot(pop,1));
        EC2_center1 = rotation(EC2_center1, EC2_center, 2, EC2_rot(pop,1));
        EC2_center2 = rotation(EC2_center2, EC2_center, 2, EC2_rot(pop,1));

        EC3.vertices = rotation(EC3.vertices, EC3_center, 2, theta3(pop,1)+EC3_rot(pop,1));
        EC3_center1 = rotation(EC3_center1, EC3_center, 2, EC3_rot(pop,1));

        EC4.vertices = rotation(EC4.vertices, EC4_center, 2, theta4(pop,1)+EC4_rot(pop,1));
        EC4_center1 = rotation(EC4_center1, EC4_center, 2, EC4_rot(pop,1));

        %Sort chromosomes based on their fitness values
        [srt,I]=sort(d);

%         
        for i=1:pop/2
            bestgb(i,1) = gb_rot(I(i,1),1);
            best1(i,1) = r1(I(i,1),:);
            best2(i,1) = r2(I(i,1),:);
            best3(i,1) = r3(I(i,1),:);
            best4(i,1) = r4(I(i,1),:);
            best1(i,2) = y1(I(i,1),:);
            best2(i,2) = y2(I(i,1),:);
            best3(i,2) = y3(I(i,1),:);
            best4(i,2) = y4(I(i,1),:);
            best1(i,3) = theta1(I(i,1),:);
            best2(i,3) = theta2(I(i,1),:);
            best3(i,3) = theta3(I(i,1),:);
            best4(i,3) = theta4(I(i,1),:);
            best_rot1(i,1) = EC1_rot(I(i,1),1);
            best_rot2(i,1) = EC2_rot(I(i,1),1);
            best_rot3(i,1) = EC3_rot(I(i,1),1);
            best_rot4(i,1) = EC4_rot(I(i,1),1);
            bestd(i,1) = d(I(i,1),1);
            best_placement(i,1) = EC2_placement(I(i,1),1);
        end
        
            if (g==iter)
        %Prints best solution(of last generation)
            %Gearbox
            fv.vertices = rotation(fv.vertices, Anchor, 2, gb_rot(I(1,1),1));
            EM_start = rotation(EM_start, Anchor, 2, gb_rot(I(1,1),1));
            EM_end = EM_start + EM_vector;

            %IGBT
            EC1.vertices = EC1.vertices - EC1_center + coords1(I(1,1),:);
            EC1_center1 = EC1_center1 - EC1_center + coords1(I(1,1),:);
            EC1_center2 = EC1_center2 - EC1_center + coords1(I(1,1),:);
            EC1_center3 = EC1_center3 - EC1_center + coords1(I(1,1),:);
            EC1_center = coords1(I(1,1),:);
            EC1.vertices = rotation(EC1.vertices, EC1_center, 2, -theta1(I(1,1),1) - EC1_rot(I(1,1),1));
            EC1_center1 = rotation(EC1_center1, EC1_center, 2, -theta1(I(1,1),1) - EC1_rot(I(1,1),1));
            EC1_center2 = rotation(EC1_center2, EC1_center, 2, -theta1(I(1,1),1) - EC1_rot(I(1,1),1));
            EC1_center3 = rotation(EC1_center3, EC1_center, 2, -theta1(I(1,1),1) - EC1_rot(I(1,1),1));

            %Capacitor
            EC2.vertices = EC2.vertices - EC2_center + coords2(I(1,1),:);
            EC2_center1 = EC2_center1 - EC2_center + coords2(I(1,1),:);
            EC2_center2 = EC2_center2 - EC2_center + coords2(I(1,1),:);
            EC2_center = coords2(I(1,1),:);
            EC2.vertices = rotation(EC2.vertices, EC2_center, 2, -theta2(I(1,1),1) - EC2_rot(I(1,1),1));
            EC2_center1 = rotation(EC2_center1, EC2_center, 2, -theta2(I(1,1),1) - EC2_rot(I(1,1),1));
            EC2_center2 = rotation(EC2_center2, EC2_center, 2, -theta2(I(1,1),1) - EC2_rot(I(1,1),1));

            %EMI-filter
            EC3.vertices = EC3.vertices - EC3_center + coords3(I(1,1),:);
            EC3_center1 = EC3_center1 - EC3_center + coords3(I(1,1),:);
            EC3_center = coords3(I(1,1),:);
            EC3.vertices = rotation(EC3.vertices, EC3_center, 2, -theta3(I(1,1),1) - EC3_rot(I(1,1),1));
            EC3_center1 = rotation(EC3_center1, EC3_center, 2, -theta3(I(1,1),1) - EC3_rot(I(1,1),1));

            %Control board
            EC4.vertices = EC4.vertices - EC4_center + coords4(I(1,1),:);
            EC4_center1 = EC4_center1 - EC4_center + coords4(I(1,1),:);
            EC4_center = coords4(I(1,1),:);
            EC4.vertices = rotation(EC4.vertices, EC4_center, 2, -theta4(I(1,1),1) - EC4_rot(I(1,1),1));
            EC4_center1 = rotation(EC4_center1, EC4_center, 2, -theta4(I(1,1),1) - EC4_rot(I(1,1),1));
            

            hold on    
            figure(pop+1)
            BRaum = patch(bauraum,'FaceColor',       [0 0 1.0], ...
                     'EdgeColor',       'none',        ...
                     'FaceLighting',    'gouraud',     ...
                     'AmbientStrength', 0.15);
            BRaum.FaceVertexAlphaData = 0.2;    % Set constant transparency 
            BRaum.FaceAlpha = 'flat' ;          % Interpolate to find face transparency
            hold on
            view([45 90 135]);
            xlabel('X') 
            ylabel('Y')
            zlabel('Z')
            
            
            hold on     
            patch(fv,'FaceColor',       [0.8 0.8 1], ...
                 'FaceLighting',    'gouraud',     ...
                 'AmbientStrength', 0.15);
            plot3(shaftx', shafty', shaftz','LineWidth',2,'Color','red')
            
            
            patch(EC1,'FaceColor',       [0.8 0 0], ...
                     'EdgeColor',       'green',        ...
                     'FaceLighting',    'gouraud',     ...
                     'AmbientStrength', 0.15);
            patch(EC2,'FaceColor',       [0 0.8 0], ...
                 'EdgeColor',       'red',        ...
                 'FaceLighting',    'gouraud',     ...
                 'AmbientStrength', 0.15);
             patch(EC3,'FaceColor',       [1  1 0], ...
                     'EdgeColor',       'black',        ...
                     'FaceLighting',    'gouraud',     ...
                     'AmbientStrength', 0.15);
            patch(EC4,'FaceColor',       [1 0 1], ...
                 'EdgeColor',       'yellow',        ...
                 'FaceLighting',    'gouraud',     ...
                 'AmbientStrength', 0.15); 
            end
       fprintf('Selection coefficient of Generation %d: %d \n',g-1,sel_coeff);      
       fprintf('Fitness of Generation %d: %d \n' ,g,srt(1,1));

 end

tEnd = cputime - tStart

function fitness = optimization_function(EC1_center, EC1_center1, EC1_center2, EC1_center3, EC2_center1,EC2_center2, EC3_center1, EC3_center2, EC4_center1, E_start, E_end, in1_bauraum, in2_bauraum, in3_bauraum, in4_bauraum, out1_gearbox,  out2_gearbox, out3_gearbox, out4_gearbox, collision12, collision13, collision14, collision23, collision24, collision34, EC2_placement) 
    Pen_gearbox = 10000*5;
    Pen_bauraum = 5000*5;
    Pen_collision = 5000*5;

    %IGBT to Electric motor
    if EC1_center(1,2) > E_start(1,2) || EC1_center(1,2) < E_end(1,2)
         distance = point_to_line(EC1_center, E_start,E_end);
         fitness1 = distance + (1-in1_bauraum)*Pen_bauraum +(out1_gearbox)*Pen_gearbox;
     else
         distance = point_to_line(EC1_center, E_start,E_end);
         fitness1 = distance - 135.85 + (1-in1_bauraum)*Pen_bauraum +(out1_gearbox)*Pen_gearbox;
    end
 
    %B to Electric Motor
%     if EC2_center(1,2) > E_start(1,2) || EC2_center(1,2) < E_end(1,2)
%          distance = point_to_line(EC2_center, E_start,E_end);
%          fitness2 = distance + (1-in2_bauraum+out2_gearbox)*Pen;
%      else
%          distance = point_to_line(EC2_center, E_start,E_end);
%          fitness2 = distance - 135.85 + (1-in2_bauraum+out2_gearbox)*Pen;
%     end  
    
    %Counterclockwise placement of Condensator
    if EC2_placement == 1       
        %IGBT to condensator
            distance = norm(EC1_center1 - EC2_center2);
            fitness2 = distance + (1-in2_bauraum)*Pen_bauraum +(out2_gearbox)*Pen_gearbox;
        %Condensator to EMI-filter    
            distance = norm(EC2_center1 - EC3_center2);
            fitness3 = distance + (1-in3_bauraum)*Pen_bauraum +(out3_gearbox)*Pen_gearbox;
            
    %Clockwise placement of Condensator
    else
        %IGBT to condensator
            distance = norm(EC1_center2 - EC2_center1);
            fitness2 = distance + (1-in2_bauraum)*Pen_bauraum +(out2_gearbox)*Pen_gearbox;
        %Condensator to EMI-filter    
            distance = norm(EC2_center2 - EC3_center1);
            fitness3 = distance + (1-in3_bauraum)*Pen_bauraum +(out3_gearbox)*Pen_gearbox;
    end
        
%IGBT to Control board    
    distance = norm(EC1_center3 - EC4_center1);
    fitness4 = distance + (1-in4_bauraum)*Pen_bauraum +(out4_gearbox)*Pen_gearbox;
    
    
    fitness = fitness1 + fitness2 + fitness3 + fitness4 + Pen_collision*collision12 + Pen_collision*collision13 + Pen_collision*collision14 + Pen_collision*collision23 + Pen_collision*collision24 + Pen_collision*collision34;

 end

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