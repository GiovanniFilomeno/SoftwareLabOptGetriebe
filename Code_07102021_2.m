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

%% Defining the limits of the bauraum to generate the population inside.
xmax = max(bauraum.vertices(:,1));
xmin = min(bauraum.vertices(:,1));
ymax = max(bauraum.vertices(:,2));
ymin = min(bauraum.vertices(:,2));
zmax = max(bauraum.vertices(:,3));
zmin = min(bauraum.vertices(:,3));


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

%To use with Electric motor
EC1_center1 = [EC1_xmax  (EC1_ymax+EC1_ymin)/2 (EC1_zmax + EC1_zmin)/2]; 

%To use with Condensator
EC1_center2 = [EC1_xmin  (EC1_ymax+EC1_ymin)/2 (EC1_zmax + EC1_zmin)/2]; 

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

%To use with IGBT
EC2_center1 = [EC2_xmax (EC2_ymax+EC2_ymin)/2 (EC2_zmax + EC2_zmin)/2];

%To use with EMI-filter
EC2_center2 = [EC2_xmin (EC2_ymax+EC2_ymin)/2 (EC2_zmax + EC2_zmin)/2];


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

%To use with Condensator
EC3_center1 = [EC3_xmax (EC3_ymax+EC3_ymin)/2 (EC3_zmax + EC3_zmin)/2]; 

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
pop = 200;
iter = 10;


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
x1 = xmin + (xmax-xmin).*rand(pop,1);
y1= ymin + (ymax-ymin).*rand(pop,1);
z1 = zmin + (zmax-zmin).*rand(pop,1);
coords1 = [x1,y1,z1];
EC1_rot = -pi/4 + (pi/2).*rand(pop,1);
EC1_mov = zeros(pop,1);
for i=1:pop
    if i==1
        EC1_mov(i) = EC1_rot(i);
    else
        EC1_mov(i) = EC1_rot(i) - EC1_rot(i-1);
    end
end


%Condensator
x2 = xmin + (xmax-xmin).*rand(pop,1);
y2= ymin + (ymax-ymin).*rand(pop,1);
z2 = zmin + (zmax-zmin).*rand(pop,1);
coords2 = [x2,y2,z2];
EC2_rot = -pi/4 + (pi/2).*rand(pop,1);
EC2_mov = zeros(pop,1);
for i=1:pop
    if i==1
        EC2_mov(i) = EC2_rot(i);
    else
        EC2_mov(i) = EC2_rot(i) - EC2_rot(i-1);
    end
end


%EMI-Filter
x3 = xmin + (xmax-xmin).*rand(pop,1);
y3= ymin + (ymax-ymin).*rand(pop,1);
z3 = zmin + (zmax-zmin).*rand(pop,1);
coords3 = [x3,y3,z3];
EC3_rot = -pi/4 + (pi/2).*rand(pop,1);
EC3_mov = zeros(pop,1);
for i=1:pop
    if i==1
        EC3_mov(i) = EC3_rot(i);
    else
        EC3_mov(i) = EC3_rot(i) - EC3_rot(i-1);
    end
end


%Control board
x4 = xmin + (xmax-xmin).*rand(pop,1);
y4= ymin + (ymax-ymin).*rand(pop,1);
z4 = zmin + (zmax-zmin).*rand(pop,1);
coords4 = [x4,y4,z4];
EC4_rot = -pi/4 + (pi/2).*rand(pop,1);
EC4_mov = zeros(pop,1);
for i=1:pop
    if i==1
        EC4_mov(i) = EC4_rot(i);
    else
        EC4_mov(i) = EC4_rot(i) - EC4_rot(i-1);
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
    
    %IGBT
    EC1.vertices = EC1.vertices - EC1_center + coords1(i,:);
    EC1_center1 = EC1_center1 - EC1_center + coords1(i,:);
    EC1_center2 = EC1_center2 - EC1_center + coords1(i,:);
    EC1_center3 = EC1_center3 - EC1_center + coords1(i,:);
    EC1_center = coords1(i,:);
    EC1.vertices = rotation(EC1.vertices, EC1_center, 2, EC1_mov(i,1));
    EC1_center1 = rotation(EC1_center1, EC1_center, 2, EC1_mov(i,1));
    EC1_center2 = rotation(EC1_center2, EC1_center, 2, EC1_mov(i,1));
    EC1_center3 = rotation(EC1_center3, EC1_center, 2, EC1_mov(i,1));
    
    %Capacitor
    EC2.vertices = EC2.vertices - EC2_center + coords2(i,:);
    EC2_center1 = EC2_center1 - EC2_center + coords2(i,:);
    EC2_center2 = EC2_center2 - EC2_center + coords2(i,:);
    EC2_center = coords2(i,:);
    EC2.vertices = rotation(EC2.vertices, EC2_center, 2, EC2_mov(i,1));
    EC2_center1 = rotation(EC2_center1, EC2_center, 2, EC2_mov(i,1));
    EC2_center2 = rotation(EC2_center2, EC2_center, 2, EC2_mov(i,1));
    
    %EMI-filter
    EC3.vertices = EC3.vertices - EC3_center + coords3(i,:);
    EC3_center1 = EC3_center1 - EC3_center + coords3(i,:);
    EC3_center = coords3(i,:);
    EC3.vertices = rotation(EC3.vertices, EC3_center, 2, EC3_mov(i,1));
    EC3_center1 = rotation(EC3_center1, EC3_center, 2, EC3_mov(i,1));
    
    %Control board
    EC4.vertices = EC4.vertices - EC4_center + coords4(i,:);
    EC4_center1 = EC4_center1 - EC4_center + coords4(i,:);
    EC4_center = coords4(i,:);
    EC4.vertices = rotation(EC4.vertices, EC4_center, 2, EC4_mov(i,1));
    EC4_center1 = rotation(EC4_center1, EC4_center, 2, EC4_mov(i,1));
    

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
    

    d(i,1) = optimization_function(EC1_center1, EC1_center2, EC1_center3, EC2_center1,EC2_center2, EC3_center1, EC4_center1, EM_start, EM_end, in1_bauraum, in2_bauraum, in3_bauraum, in4_bauraum, out1_gearbox,  out2_gearbox, out3_gearbox, out4_gearbox, collision12, collision13, collision14, collision23, collision24, collision34);

end

%Return the gearbox to initial angle
fv.vertices = rotation(fv.vertices, Anchor, 2, -gb_rot(pop,1));
EM_start = rotation(EM_start, Anchor, 2, -gb_rot(pop,1));
EM_end = EM_start + EM_vector;

%Return components to initial angle
EC1.vertices = rotation(EC1.vertices, EC1_center, 2, -EC1_rot(pop,1));
EC1_center1 = rotation(EC1_center1, EC1_center, 2, -EC1_rot(pop,1));
EC1_center2 = rotation(EC1_center2, EC1_center, 2, -EC1_rot(pop,1));
EC1_center3 = rotation(EC1_center3, EC1_center, 2, -EC1_rot(pop,1));

EC2.vertices = rotation(EC2.vertices, EC2_center, 2, -EC2_rot(pop,1));
EC2_center1 = rotation(EC2_center1, EC2_center, 2, -EC2_rot(pop,1));
EC2_center2 = rotation(EC2_center2, EC2_center, 2, -EC2_rot(pop,1));

EC3.vertices = rotation(EC3.vertices, EC3_center, 2, -EC3_rot(pop,1));
EC3_center1 = rotation(EC3_center1, EC3_center, 2, -EC3_rot(pop,1));

EC4.vertices = rotation(EC4.vertices, EC4_center, 2, -EC4_rot(pop,1));
EC4_center1 = rotation(EC4_center1, EC4_center, 2, -EC4_rot(pop,1));

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
for i=1:pop/2
    bestgb(i,1) = gb_rot(I(i,1),1);
    best1(i,:) = coords1(I(i,1),:);
    best2(i,:) = coords2(I(i,1),:);
    best3(i,:) = coords3(I(i,1),:);
    best4(i,:) = coords4(I(i,1),:);
    best_rot1(i,1) = EC1_rot(I(i,1),1);
    best_rot2(i,1) = EC2_rot(I(i,1),1);
    best_rot3(i,1) = EC3_rot(I(i,1),1);
    best_rot4(i,1) = EC4_rot(I(i,1),1);
end



 for g = 2:iter
      for k = 1: pop/2
          
          
        %Selection
        sel_coeff = mean(d);
        parents1 = zeros(2,3);
        parents2 = zeros(2,3);
        parents3 = zeros(2,3);
        parents4 = zeros(2,3);
        parents_rot = zeros(2,5);
        for j=1:2
            rand_sel = 1 + round(pop.*rand());
            for i = rand_sel:pop
                if d(i,1)<sel_coeff
                    parents1(j,1) = coords1(i,1);
                    parents1(j,2) = coords1(i,2);
                    parents1(j,3) = coords1(i,3);
                    parents2(j,1) = coords2(i,1);
                    parents2(j,2) = coords2(i,2);
                    parents2(j,3) = coords2(i,3);
                    parents3(j,1) = coords3(i,1);
                    parents3(j,2) = coords3(i,2);
                    parents3(j,3) = coords3(i,3);
                    parents4(j,1) = coords4(i,1);
                    parents4(j,2) = coords4(i,2);
                    parents4(j,3) = coords4(i,3);
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
        child2 = zeros(2,3);
        child3 = zeros(2,3);
        child4 = zeros(2,3);
        child_rot = zeros(2,5);
        alpha = 0.15;
        beta = 0.85;
        child1(1,:) = alpha.*parents1(1,:) + (1-alpha).*parents1(2,:);
        child1(2,:) = beta.*parents1(1,:) + (1-beta).*parents1(2,:);
        child2(1,:) = alpha.*parents2(1,:) + (1-alpha).*parents2(2,:);
        child2(2,:) = beta.*parents2(1,:) + (1-beta).*parents2(2,:);
        child3(1,:) = alpha.*parents3(1,:) + (1-alpha).*parents3(2,:);
        child3(2,:) = beta.*parents3(1,:) + (1-beta).*parents3(2,:);
        child4(1,:) = alpha.*parents4(1,:) + (1-alpha).*parents4(2,:);
        child4(2,:) = beta.*parents4(1,:) + (1-beta).*parents4(2,:);
        child_rot(1,:) = alpha*parents_rot(1,:) + (1-alpha)*parents_rot(2,:);
        child_rot(2,:) = beta*parents_rot(1,:) + (1-beta)*parents_rot(2,:);
        


%       Mutation
        for i = 1:2
             child1_start = EM_start - child1(i,:);
             child1_end = EM_end - child1(i,:);
             for j = 1:3
                mut_prob = rand();
                if mut_prob>=0.75
                    if child1(i,2) > EM_start(1,2)
                        child1(i,j) = child1(i,j) + 0.05*child1_start(1,j);
                    elseif child1(i,2) < EM_end(1,2)
                        child1(i,j) = child1(i,j) + 0.05*child1_end(1,j);
                    else
                        child1(i,j) = child1(i,j) + 0.05*child1_end(1,j);
                        child1(i,2) =child1(i,2) - 0.05*child1_end(1,2); 
                    end
                    
                    child2(i,j) = child2(i,j) + 0.05*(child1(i,j)-child2(i,j));
                    child3(i,j) = child3(i,j) + 0.05*(child2(i,j)-child3(i,j));
                    child4(i,j) = child4(i,j) + 0.05*(child1(i,j)-child4(i,j));
                end
             end
             for j=1:5
                 mut_prob = rand();
                 if mut_prob > 0.5 
                     child_rot(i,j) = child_rot(i,j) + (pi/180);
                 else
                     child_rot(i,j) = child_rot(i,j) - (pi/180);
                 end
                 
             end
        end
        
        %Write first half of new generation 
        coords1(k,1) = child1(1,1);
        coords1(k,2) = child1(1,2);
        coords1(k,3) = child1(1,3);
        coords1(k+1,1) = child1(2,1);
        coords1(k+1,2) = child1(2,2);
        coords1(k+1,3) = child1(2,3);
        coords2(k,1) = child2(1,1);
        coords2(k,2) = child2(1,2);
        coords2(k,3) = child2(1,3);
        coords2(k+1,1) = child2(2,1);
        coords2(k+1,2) = child2(2,2);
        coords2(k+1,3) = child2(2,3);
        coords3(k,1) = child3(1,1);
        coords3(k,2) = child3(1,2);
        coords3(k,3) = child3(1,3);
        coords3(k+1,1) = child3(2,1);
        coords3(k+1,2) = child3(2,2);
        coords3(k+1,3) = child3(2,3);
        coords4(k,1) = child4(1,1);
        coords4(k,2) = child4(1,2);
        coords4(k,3) = child4(1,3);
        coords4(k+1,1) = child4(2,1);
        coords4(k+1,2) = child4(2,2);
        coords4(k+1,3) = child4(2,3);
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
      
%Write othter half of new population from the best of previous generation     
      for i=(pop/2+1):pop
          coords1(i,:) = best1(i-(pop/2),:);
          coords2(i,:) = best2(i-(pop/2),:);
          coords3(i,:) = best3(i-(pop/2),:);
          coords4(i,:) = best4(i-(pop/2),:);
          gb_rot(i,1) = bestgb(i-(pop/2),1);
          EC1_rot(i,1) = best_rot1(i-(pop/2),1);
          EC2_rot(i,1) = best_rot2(i-(pop/2),1);
          EC3_rot(i,1) = best_rot3(i-(pop/2),1);
          EC4_rot(i,1) = best_rot4(i-(pop/2),1);
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
            EC1_mov(i) = EC1_rot(i);
        else
            EC1_mov(i) = EC1_rot(i) - EC1_rot(i-1);
        end
      end
      
      for i=1:pop
        if i==1
            EC2_mov(i) = EC2_rot(i);
        else
            EC2_mov(i) = EC2_rot(i) - EC2_rot(i-1);
        end
      end
      
      for i=1:pop
        if i==1
            EC3_mov(i) = EC3_rot(i);
        else
            EC3_mov(i) = EC3_rot(i) - EC3_rot(i-1);
        end
      end
      
      for i=1:pop
        if i==1
            EC4_mov(i) = EC4_rot(i);
        else
            EC4_mov(i) = EC4_rot(i) - EC4_rot(i-1);
        end
      end
      
%Test new generation
      for i=1:pop
            %Gearbox
            fv.vertices = rotation(fv.vertices, Anchor, 2, gb_mov(i,1));
            EM_start = rotation(EM_start, Anchor, 2, gb_mov(i,1));
            EM_end = EM_start + EM_vector;

            %IGBT
            EC1.vertices = EC1.vertices - EC1_center + coords1(i,:);
            EC1_center1 = EC1_center1 - EC1_center + coords1(i,:);
            EC1_center2 = EC1_center2 - EC1_center + coords1(i,:);
            EC1_center3 = EC1_center3 - EC1_center + coords1(i,:);
            EC1_center = coords1(i,:);
            EC1.vertices = rotation(EC1.vertices, EC1_center, 2, EC1_mov(i,1));
            EC1_center1 = rotation(EC1_center1, EC1_center, 2, EC1_mov(i,1));
            EC1_center2 = rotation(EC1_center2, EC1_center, 2, EC1_mov(i,1));
            EC1_center3 = rotation(EC1_center3, EC1_center, 2, EC1_mov(i,1));

            %Capacitor
            EC2.vertices = EC2.vertices - EC2_center + coords2(i,:);
            EC2_center1 = EC2_center1 - EC2_center + coords2(i,:);
            EC2_center2 = EC2_center2 - EC2_center + coords2(i,:);
            EC2_center = coords2(i,:);
            EC2.vertices = rotation(EC2.vertices, EC2_center, 2, EC2_mov(i,1));
            EC2_center1 = rotation(EC2_center1, EC2_center, 2, EC2_mov(i,1));
            EC2_center2 = rotation(EC2_center2, EC2_center, 2, EC2_mov(i,1));

            %EMI-filter
            EC3.vertices = EC3.vertices - EC3_center + coords3(i,:);
            EC3_center1 = EC3_center1 - EC3_center + coords3(i,:);
            EC3_center = coords3(i,:);
            EC3.vertices = rotation(EC3.vertices, EC3_center, 2, EC3_mov(i,1));
            EC3_center1 = rotation(EC3_center1, EC3_center, 2, EC3_mov(i,1));

            %Control board
            EC4.vertices = EC4.vertices - EC4_center + coords4(i,:);
            EC4_center1 = EC4_center1 - EC4_center + coords4(i,:);
            EC4_center = coords4(i,:);
            EC4.vertices = rotation(EC4.vertices, EC4_center, 2, EC4_mov(i,1));
            EC4_center1 = rotation(EC4_center1, EC4_center, 2, EC4_mov(i,1));


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
            
            d(i,1) = optimization_function(EC1_center1, EC1_center2, EC1_center3, EC2_center1,EC2_center2, EC3_center1, EC4_center1, EM_start, EM_end, in1_bauraum, in2_bauraum, in3_bauraum, in4_bauraum, out1_gearbox,  out2_gearbox, out3_gearbox, out4_gearbox, collision12, collision13, collision14, collision23, collision24, collision34);

      end
      
        %Return the gearbox to initial angle
        fv.vertices = rotation(fv.vertices, Anchor, 2, -gb_rot(pop,1));
        EM_start = rotation(EM_start, Anchor, 2, -gb_rot(pop,1));
        EM_end = EM_start + EM_vector;

        %Return components to initial angle
        EC1.vertices = rotation(EC1.vertices, EC1_center, 2, -EC1_rot(pop,1));
        EC1_center1 = rotation(EC1_center1, EC1_center, 2, -EC1_rot(pop,1));
        EC1_center2 = rotation(EC1_center2, EC1_center, 2, -EC1_rot(pop,1));
        EC1_center3 = rotation(EC1_center3, EC1_center, 2, -EC1_rot(pop,1));

        EC2.vertices = rotation(EC2.vertices, EC2_center, 2, -EC2_rot(pop,1));
        EC2_center1 = rotation(EC2_center1, EC2_center, 2, -EC2_rot(pop,1));
        EC2_center2 = rotation(EC2_center2, EC2_center, 2, -EC2_rot(pop,1));

        EC3.vertices = rotation(EC3.vertices, EC3_center, 2, -EC3_rot(pop,1));
        EC3_center1 = rotation(EC3_center1, EC3_center, 2, -EC3_rot(pop,1));

        EC4.vertices = rotation(EC4.vertices, EC4_center, 2, -EC4_rot(pop,1));
        EC4_center1 = rotation(EC4_center1, EC4_center, 2, -EC4_rot(pop,1));
      
        [srt,I]=sort(d);
        
        for i=1:pop/2
            bestgb(i,1) = gb_rot(I(i,1),1);
            best1(i,:) = coords1(I(i,1),:);
            best2(i,:) = coords2(I(i,1),:);
            best3(i,:) = coords3(I(i,1),:);
            best4(i,:) = coords4(I(i,1),:);
            best_rot1(i,1) = EC1_rot(I(i,1),1);
            best_rot2(i,1) = EC2_rot(I(i,1),1);
            best_rot3(i,1) = EC3_rot(I(i,1),1);
            best_rot4(i,1) = EC4_rot(I(i,1),1);
        end
        
        
        %Prints best solution(of last generation)
        if g == iter
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
            EC1.vertices = rotation(EC1.vertices, EC1_center, 2, EC1_rot(I(1,1),1));
            EC1_center1 = rotation(EC1_center1, EC1_center, 2, EC1_rot(I(1,1),1));
            EC1_center2 = rotation(EC1_center2, EC1_center, 2, EC1_rot(I(1,1),1));
            EC1_center3 = rotation(EC1_center3, EC1_center, 2, EC1_rot(I(1,1),1));

            %Capacitor
            EC2.vertices = EC2.vertices - EC2_center + coords2(I(1,1),:);
            EC2_center1 = EC2_center1 - EC2_center + coords2(I(1,1),:);
            EC2_center2 = EC2_center2 - EC2_center + coords2(I(1,1),:);
            EC2_center = coords2(I(1,1),:);
            EC2.vertices = rotation(EC2.vertices, EC2_center, 2, EC2_rot(I(1,1),1));
            EC2_center1 = rotation(EC2_center1, EC2_center, 2, EC2_rot(I(1,1),1));
            EC2_center2 = rotation(EC2_center2, EC2_center, 2, EC2_rot(I(1,1),1));

            %EMI-filter
            EC3.vertices = EC3.vertices - EC3_center + coords3(I(1,1),:);
            EC3_center1 = EC3_center1 - EC3_center + coords3(I(1,1),:);
            EC3_center = coords3(I(1,1),:);
            EC3.vertices = rotation(EC3.vertices, EC3_center, 2, EC3_rot(I(1,1),1));
            EC3_center1 = rotation(EC3_center1, EC3_center, 2, EC3_rot(I(1,1),1));

            %Control board
            EC4.vertices = EC4.vertices - EC4_center + coords4(I(1,1),:);
            EC4_center1 = EC4_center1 - EC4_center + coords4(I(1,1),:);
            EC4_center = coords4(I(1,1),:);
            EC4.vertices = rotation(EC4.vertices, EC4_center, 2, EC4_rot(I(1,1),1));
            EC4_center1 = rotation(EC4_center1, EC4_center, 2, EC4_rot(I(1,1),1));
            
            
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
scatter3(EM_end(1),EM_end(2),EM_end(3));
rad_gb = EM_radius;
%             %igbt
%             EC1_xmax = max(EC1.vertices(:,1));
%             EC1_xmin = min(EC1.vertices(:,1));
%             EC1_ymax = max(EC1.vertices(:,2));
%             EC1_ymin = min(EC1.vertices(:,2));
%             EC1_zmax = max(EC1.vertices(:,3));
%             EC1_zmin = min(EC1.vertices(:,3));
%             
%             %condenser
%             EC2_xmax = max(EC2.vertices(:,1));
%             EC2_xmin = min(EC2.vertices(:,1));
%             EC2_ymax = max(EC2.vertices(:,2));
%             EC2_ymin = min(EC2.vertices(:,2));
%             EC2_zmax = max(EC2.vertices(:,3));
%             EC2_zmin = min(EC2.vertices(:,3));
%             
%             %EMI-filter
%             EC3_xmax = max(EC3.vertices(:,1));
%             EC3_xmin = min(EC3.vertices(:,1));
%             EC3_ymax = max(EC3.vertices(:,2));
%             EC3_ymin = min(EC3.vertices(:,2));
%             EC3_zmax = max(EC3.vertices(:,3));
%             EC3_zmin = min(EC3.vertices(:,3));
%             
%             %Control board
%             EC4_xmax = max(EC4.vertices(:,1));
%             EC4_xmin = min(EC4.vertices(:,1));
%             EC4_ymax = max(EC4.vertices(:,2));
%             EC4_ymin = min(EC4.vertices(:,2));
%             EC4_zmax = max(EC4.vertices(:,3));
%             EC4_zmin = min(EC4.vertices(:,3));
            
            clearance = 5;
            % igbt to motor
            EM_connection = EM_end;
            %EC1_center1 = [(EC1_xmax + EC1_xmin)/2 EC1_ymax (EC1_zmax + EC1_zmin)/2];
            %EC1_center2 = [(EC1_xmax + EC1_xmin)/2 EC1_ymin (EC1_zmax + EC1_zmin)/2];
            %pt3
            vector_pt1 = [EC1_center2(1)+(-clearance) EM_connection(2) EC1_center2(3)];
            vector_pt2 = EM_connection+[0 -clearance 0]; 
            vector = vector_pt1 - vector_pt2;
            vector = vector / norm(vector);
            %pt3 = [(EC1_xmax + EC1_xmin)/2 EM_connection(2)+(-clearance) EM_connection(3)+(50)];
            pt3 = EM_connection+[0 -clearance 0]+rad_gb*vector;
            im = [EM_connection; EM_connection+[0 -clearance 0]; pt3; EC1_center2+[-clearance 0 0]; EC1_center2];
            plot3(im(:,1),im(:,2), im(:,3),'yo','LineWidth',2);
%             text(x,y,z,[repmat('  ',4,1), num2str((1:4)')])
            ax = gca;
            ax.XTick = [];
            ax.YTick = [];
            ax.ZTick = [];
            hold on
            im_spline = cscvn(im');
            fnplt(im_spline,'y',2)
            im_fnprime = fnder(im_spline);
            Lfun = @(s) sqrt(sum(fnval(im_fnprime,s).^2,1));
            L_im = integral(Lfun,im_spline.breaks(1),im_spline.breaks(end));
            hold on
            
            %connection between igbt and condensator
            %EC1_center1 = [(EC1_xmax + EC1_xmin)/2 EC1_ymax (EC1_zmax + EC1_zmin)/2];
            %EC1_center2 = [(EC1_xmax + EC1_xmin)/2 EC1_ymin (EC1_zmax + EC1_zmin)/2];
            %EC2_center1 = [(EC2_xmax + EC2_xmin)/2 EC2_ymax  (EC2_zmax + EC2_zmin)/2];
            ic = [EC1_center1; EC1_center1+[clearance 0 0]; EC2_center1+[clearance 0 0]; EC2_center1];
            plot3(ic(:,1),ic(:,2), ic(:,3),'ro','LineWidth',2);
%             text(x,y,z,[repmat('  ',4,1), num2str((1:4)')])
            ax = gca;
            ax.XTick = [];
            ax.YTick = [];
            ax.ZTick = [];
            hold on
            ic_spline = cscvn(ic');
            fnplt(ic_spline,'r',2)
            ic_fnprime = fnder(ic_spline);
            Lfun = @(s) sqrt(sum(fnval(ic_fnprime,s).^2,1));
            L_ic = integral(Lfun,ic_spline.breaks(1),ic_spline.breaks(end));
            hold on
            
            %connection between condensator and emi
            %EC2_center2 = [(EC2_xmax + EC2_xmin)/2 EC2_ymin (EC2_zmax + EC2_zmin)/2];
            %EC3_center1 = [(EC3_xmax + EC3_xmin)/2 EC3_ymax (EC3_zmax + EC3_zmin)/2];
            ce = [EC2_center2; EC2_center2+[-clearance 0 0]; EC3_center1+[clearance 0 0]; EC3_center1];
            plot3(ce(:,1),ce(:,2), ce(:,3),'go','LineWidth',2);
%             text(x,y,z,[repmat('  ',4,1), num2str((1:4)')])
            ax = gca;
            ax.XTick = [];
            ax.YTick = [];
            ax.ZTick = [];
            hold on
            ce_spline = cscvn(ce');
            fnplt(ce_spline,'g',2)
            ce_fnprime = fnder(ce_spline);
            Lfun = @(s) sqrt(sum(fnval(ce_fnprime,s).^2,1));
            L_ce = integral(Lfun,ce_spline.breaks(1),ce_spline.breaks(end));
            hold on
            
            %connection between igbt and control board
            %EC1_center3 = [(EC1_xmax + EC1_xmin)/2 (EC1_ymax + EC1_ymin)/2 EC1_zmax]; 
            %EC4_center1 = [(EC4_xmax + EC4_xmin)/2 (EC4_ymax + EC4_ymin)/2 EC4_zmin]; 
            icb = [EC1_center3; EC1_center3+[0 0 +clearance]; EC4_center1+[0 0 -clearance]; EC4_center1];
            plot3(icb(:,1),icb(:,2), icb(:,3),'bo','LineWidth',2);
%             text(x,y,z,[repmat('  ',4,1), num2str((1:4)')])
            ax = gca;
            ax.XTick = [];
            ax.YTick = [];
            ax.ZTick = [];
            hold on
            icb_spline = cscvn(icb');
            fnplt(icb_spline,'b',2)
            icb_fnprime = fnder(icb_spline);
            Lfun = @(s) sqrt(sum(fnval(icb_fnprime,s).^2,1));
            L_icb = integral(Lfun,icb_spline.breaks(1),icb_spline.breaks(end));
            hold on
            
            L_Total = L_im+L_ic+L_ce+L_icb

tEnd = cputime - tStart

function fitness = optimization_function(EC1_center1, EC1_center2, EC1_center3, EC2_center1,EC2_center2, EC3_center1, EC4_center1, E_start, E_end, in1_bauraum, in2_bauraum, in3_bauraum, in4_bauraum, out1_gearbox,  out2_gearbox, out3_gearbox, out4_gearbox, collision12, collision13, collision14, collision23, collision24, collision34) 
    Pen = 5000;

    %IGBT to Electric motor
    if EC1_center1(1,2) > E_start(1,2) || EC1_center1(1,2) < E_end(1,2)
         distance = point_to_line(EC1_center1, E_start,E_end);
         fitness1 = distance + (1-in1_bauraum+out1_gearbox)*Pen;
     else
         distance = point_to_line(EC1_center1, E_start,E_end);
         fitness1 = distance - 135.85 + (1-in1_bauraum+out1_gearbox)*Pen;
    end
 
    %B to Electric Motor
%     if EC2_center(1,2) > E_start(1,2) || EC2_center(1,2) < E_end(1,2)
%          distance = point_to_line(EC2_center, E_start,E_end);
%          fitness2 = distance + (1-in2_bauraum+out2_gearbox)*Pen;
%      else
%          distance = point_to_line(EC2_center, E_start,E_end);
%          fitness2 = distance - 135.85 + (1-in2_bauraum+out2_gearbox)*Pen;
%     end  
    
%IGBT to condensator
    distance = norm(EC1_center2 - EC2_center1);
    fitness2 = distance + (1-in2_bauraum+out2_gearbox)*Pen;

%Condensator to EMI-filter    
    distance = norm(EC2_center2 - EC3_center1);
    fitness3 = distance + (1-in3_bauraum+out3_gearbox)*Pen;
    
%IGBT to Control board    
    distance = norm(EC1_center3 - EC4_center1);
    fitness4 = distance + (1-in4_bauraum+out4_gearbox)*Pen;
    
    
    fitness = fitness1 + fitness2 + fitness3 + fitness4 + Pen*collision12 + Pen*collision13 + Pen*collision14 + Pen*collision23 + Pen*collision24 + Pen*collision34;

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