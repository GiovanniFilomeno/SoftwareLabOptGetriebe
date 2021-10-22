clear all
clc
tStart = cputime;

%% Loading all the solids
%Read Gearbox
gearbox = ReadGearbox('Gearbox1.stl');

%Read Assembly space
bauraum = ReadBauraum('Bauraum.stl');

%Read Electric components and define its limits
[EC1,EC1_length,EC1_width,EC1_height,EC1_center, EC1_SurfCenter] = ReadEC('EC1.stl');
[EC2,EC2_length,EC2_width,EC2_height,EC2_center, EC2_SurfCenter] = ReadEC('EC2.stl');
[EC3,EC3_length,EC3_width,EC3_height,EC3_center, EC3_SurfCenter] = ReadEC('EC3.stl');
[EC4,EC4_length,EC4_width,EC4_height,EC4_center, EC4_SurfCenter] = ReadEC('EC4.stl');

%% Positioning of the gearbox inside the bauraum
[EM_start, EM_vector, Anchor, EM_radius] = FindingAxis("Gearbox1.stl");
EM_end = EM_start + EM_vector;
zero = [0 0 0];

%Definition(manually) of the axis of the bauraum shaft 
point = [390 409 136.1488];

%Positioning of the gearbox inside the bauraum and obtaining angle range
[gearbox, Anchor, EM_start, EM_vector, EM_end, rad_range] = PositionGearbox(bauraum,gearbox,EM_start, EM_vector, Anchor, point, zero);

%% Genetic Algorithm implementation
pop = 500;
iter = 20;

%Create initial random population
[r, y, theta, rotations, EC2_placement ]=InitialPopulation(rad_range,pop,EM_start,EM_vector,EM_radius,EC1_height,EC1_width,EC2_height,EC2_width,EC3_height,EC3_width);

%Array that contain the fitness of the chromosome
Fitness = zeros(pop,1);

%Evaluate fitness of the first generation
Fitness = EvaluateFitness(pop,bauraum,gearbox,EC1,EC1_center,EC1_SurfCenter,EC2,EC2_center,EC2_SurfCenter,EC3,EC3_center,EC3_SurfCenter,EC4,EC4_center,EC4_SurfCenter,Anchor,EM_start,EM_vector,r,y,theta,rotations,EC2_placement, Fitness);

VisualizeBest(1, bauraum, gearbox, Anchor, EC1, EC1_center, EC2, EC2_center, EC3, EC3_center, EC4, EC4_center,EM_start,Fitness,r,y,theta,rotations)

% All_Saved_r = [];
% All_Saved_y = [];
% All_Saved_theta = [];
% All_Saved_rotations = [];

for g=2:iter

    %Make Selection, Crossover and Mutation for new generation
    [r,y,theta,rotations,EC2_placement, Fitness] = NewGeneration(pop,g,r,y,theta,rotations,EC2_placement, Fitness);

%     All_Saved_r(g-1,:)= r(:);
%     All_Saved_y(g-1,:) = y(:);
%     All_Saved_theta(g-1,:) = theta(:);
%     All_Saved_rotations(g-1,:) = rotations(:);

    
    %Evaluate new generation
    Fitness = EvaluateFitness(pop/2,bauraum,gearbox, EC1,EC1_center,EC1_SurfCenter,EC2,EC2_center,EC2_SurfCenter,EC3,EC3_center,EC3_SurfCenter,EC4,EC4_center,EC4_SurfCenter,Anchor,EM_start,EM_vector,r,y,theta,rotations,EC2_placement, Fitness);
    
    VisualizeBest(g, bauraum, gearbox, Anchor, EC1, EC1_center, EC2, EC2_center, EC3, EC3_center, EC4, EC4_center,EM_start,Fitness,r,y,theta,rotations)
end

[srt,I]=sort(Fitness);

fprintf('Fitness of Generation %d: %d \n',iter,srt(1,1));

tEnd = cputime - tStart