% Main script of the project, run to start optimization process. The script
% will ask the user to enter the population size as well as the number of
% generations that the algorithm will run. Once it runs, the script will
% display every time a generation is finished along with the fitness of the
% best individual of that generation.
% 
% Once all generations are completed, the script will generate 5 figures
% showing the best five individuals of the final generation, which are
% considered the best solutions for the optimization problem.
%
% **Example in Code**
%
% .. code-block::
%
%   prompt = 'Enter population size = ';
%   pop = input(prompt);
%   prompt2 = 'Enter number of generations = ';
%   iter = input(prompt2);

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
% disp("Enter population size = ")
prompt = 'Enter population size = ';
pop = input(prompt);
prompt2 = 'Enter number of generations = ';
iter = input(prompt2);

Best = zeros(iter,1);
Mean = zeros(iter,1);
Mean_fitness = 20000;

%Create initial random population
[r, y, theta, rotations, EC2_placement ]=InitialPopulation(rad_range,pop,EM_start,EM_vector,EM_radius,EC1_height,EC1_width,EC2_height,EC2_width,EC3_height,EC3_width);

%Array that contain the fitness of the chromosome
Fitness = zeros(pop,1);
g = 0;
%Evaluate fitness of the first generation
Fitness = EvaluateFitness(Mean_fitness, pop,bauraum,gearbox,EC1, EC1_width, EC1_center,EC1_SurfCenter,EC2,EC2_center,EC2_SurfCenter,EC3,EC3_center,EC3_SurfCenter,EC4,EC4_center,EC4_SurfCenter,Anchor,EM_start,EM_vector,r,y,theta,rotations,EC2_placement, Fitness);

for g=2:iter

    %Make Selection, Crossover and Mutation for new generation
    
    [r,y,theta,rotations,EC2_placement, Fitness, Best_fitness, Mean_fitness] = NewGeneration(pop,g,r,y,theta,rotations,EC2_placement, Fitness);

    Best(g-1,1) = Best_fitness;
    Mean(g-1,1) = Mean_fitness;

    
    %Evaluate new generation
    Fitness = EvaluateFitness(Mean_fitness, pop/2,bauraum,gearbox, EC1, EC1_width, EC1_center,EC1_SurfCenter,EC2,EC2_center,EC2_SurfCenter,EC3,EC3_center,EC3_SurfCenter,EC4,EC4_center,EC4_SurfCenter,Anchor,EM_start,EM_vector,r,y,theta,rotations,EC2_placement, Fitness);
    
    
    if g==iter
        VisualizeMany(5,g, bauraum, gearbox, Anchor, EC1, EC1_width, EC1_center,EC1_SurfCenter, EC2, EC2_center, EC2_SurfCenter, EC3, EC3_center,EC3_SurfCenter, EC4, EC4_center, EC4_SurfCenter, EM_start,EM_vector,Fitness,r,y,theta,rotations,EC2_placement)
    end
end

[srt,I]=sort(Fitness);
Best_fitness = srt(1,1);
Mean_fitness = mean(Fitness);
Worst_fitness = srt(pop,1);
Best(iter,1) = Best_fitness;
Mean(iter,1) = Mean_fitness;
Iterations = 1:iter;
fprintf('Fitness of Generation %d: %d \n',iter,srt(1,1));
%figure(iter+1)
%hold on
plot(Iterations,Best,Iterations,Mean);


tEnd = cputime - tStart