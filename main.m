clear all
clc
tStart = cputime;

%% Loading all the solids
%Read Gearbox
gearbox = ReadGearbox('Gearbox1.stl');

%Read Assembly space
bauraum = ReadBauraum('Bauraum4.stl');

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
% disp("Enter population size = ")
prompt = 'Enter population size = ';
pop = input(prompt);
prompt2 = 'Enter number of generations = ';
iter = input(prompt2);
% pop = 200;
% iter = 20;

Best = zeros(iter,1);
Mean = zeros(iter,1);
Worst = zeros(iter,1);
Mean_fitness = 0;

%Create initial random population
[r, y, theta, rotations, EC2_placement ]=InitialPopulation(rad_range,pop,EM_start,EM_vector,EM_radius,EC1_height,EC1_width,EC2_height,EC2_width,EC3_height,EC3_width);

%Array that contain the fitness of the chromosome
Fitness = zeros(pop,1);
g = 0;
%Evaluate fitness of the first generation
Fitness = EvaluateFitness(Mean_fitness, iter, g,pop,bauraum,gearbox,EC1, EC1_width, EC1_center,EC1_SurfCenter,EC2,EC2_center,EC2_SurfCenter,EC3,EC3_center,EC3_SurfCenter,EC4,EC4_center,EC4_SurfCenter,Anchor,EM_start,EM_vector,r,y,theta,rotations,EC2_placement, Fitness);

VisualizeBest(1, bauraum, gearbox, Anchor, EC1,EC1_width, EC1_center,EC1_SurfCenter, EC2, EC2_center, EC2_SurfCenter, EC3, EC3_center,EC3_SurfCenter, EC4, EC4_center, EC4_SurfCenter, EM_start,EM_vector, Fitness,r,y,theta,rotations,EC2_placement)

% All_Saved_r = [];
% All_Saved_y = [];
% All_Saved_theta = [];
% All_Saved_rotations = [];


for g=2:iter

    %Make Selection, Crossover and Mutation for new generation
    
    [r,y,theta,rotations,EC2_placement, Fitness, Best_fitness, Mean_fitness, Worst_fitness] = NewGeneration(pop,g,r,y,theta,rotations,EC2_placement, Fitness);

    Best(g-1,1) = Best_fitness;
    Mean(g-1,1) = Mean_fitness;
    Worst(g-1,1) = Worst_fitness;
%     All_Saved_r(g-1,:)= r(:);
%     All_Saved_y(g-1,:) = y(:);
%     All_Saved_theta(g-1,:) = theta(:);
%     All_Saved_rotations(g-1,:) = rotations(:);

    
    %Evaluate new generation
    Fitness = EvaluateFitness(Mean_fitness, iter, g,pop/2,bauraum,gearbox, EC1, EC1_width, EC1_center,EC1_SurfCenter,EC2,EC2_center,EC2_SurfCenter,EC3,EC3_center,EC3_SurfCenter,EC4,EC4_center,EC4_SurfCenter,Anchor,EM_start,EM_vector,r,y,theta,rotations,EC2_placement, Fitness);
    
    
   
    VisualizeBest(g, bauraum, gearbox, Anchor, EC1, EC1_width, EC1_center,EC1_SurfCenter, EC2, EC2_center, EC2_SurfCenter, EC3, EC3_center,EC3_SurfCenter, EC4, EC4_center, EC4_SurfCenter, EM_start,EM_vector,Fitness,r,y,theta,rotations,EC2_placement)

end

[srt,I]=sort(Fitness);
Best_fitness = srt(1,1);
Mean_fitness = mean(Fitness);
Worst_fitness = srt(pop,1);
Best(iter,1) = Best_fitness;
Mean(iter,1) = Mean_fitness;
Worst(iter,1) = Worst_fitness;
Iterations = 1:20;
figure(iter+1)
hold on
plot(Iterations,Best,Iterations,Mean);
fprintf('Fitness of Generation %d: %d \n',iter,srt(1,1));

tEnd = cputime - tStart