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

%Evaluate fitness of the first generation
Fitness = EvaluateFitness(pop,bauraum,gearbox,EC1,EC1_center,EC1_SurfCenter,EC2,EC2_center,EC2_SurfCenter,EC3,EC3_center,EC3_SurfCenter,EC4,EC4_center,EC4_SurfCenter,Anchor,EM_start,EM_vector,r,y,theta,rotations,EC2_placement);

for g=2:iter

    %Make Selection, Crossover and Mutation for new generation
    [r,y,theta,rotations,EC2_placement, Fitness] = NewGeneration(pop,g,r,y,theta,rotations,EC2_placement, Fitness);

    %Evaluate new generation
    Fitness = EvaluateFitness(pop,bauraum,gearbox, EC1,EC1_center,EC1_SurfCenter,EC2,EC2_center,EC2_SurfCenter,EC3,EC3_center,EC3_SurfCenter,EC4,EC4_center,EC4_SurfCenter,Anchor,EM_start,EM_vector,r,y,theta,rotations,EC2_placement);
    
    VisualizeBest(g, bauraum, gearbox, Anchor, EC1, EC1_center, EC2, EC2_center, EC3, EC3_center, EC4, EC4_center,EM_start,Fitness,r,y,theta,rotations)
end




% bestgb = zeros(pop/2,1);
% best1 = zeros(pop/2,3);
% best2 = zeros(pop/2,3);
% best3 = zeros(pop/2,3);
% best4 = zeros(pop/2,3);
% best_rot1 = zeros(pop/2,1);
% best_rot2 = zeros(pop/2,1);
% best_rot3 = zeros(pop/2,1);
% best_rot4 = zeros(pop/2,1);
% bestcoords1 = zeros(pop/2,3);
% bestcoords2 = zeros(pop/2,3);
% bestcoords3 = zeros(pop/2,3);
% bestcoords4 = zeros(pop/2,3);
% bestd = zeros(pop/2,1);
% best_placement = zeros(pop,2/1);
% for i=1:pop/2
%     bestgb(i,1) = gb_rot(I(i,1),1);
%     best1(i,1) = r1(I(i,1),:);
%     best2(i,1) = r2(I(i,1),:);
%     best3(i,1) = r3(I(i,1),:);
%     best4(i,1) = r4(I(i,1),:);
%     best1(i,2) = y1(I(i,1),:);
%     best2(i,2) = y2(I(i,1),:);
%     best3(i,2) = y3(I(i,1),:);
%     best4(i,2) = y4(I(i,1),:);
%     best1(i,3) = theta1(I(i,1),:);
%     best2(i,3) = theta2(I(i,1),:);
%     best3(i,3) = theta3(I(i,1),:);
%     best4(i,3) = theta4(I(i,1),:);
%     bestcoords1(i,:) = coords1(I(i,1),:);
%     bestcoords2(i,:) = coords2(I(i,1),:);
%     bestcoords3(i,:) = coords3(I(i,1),:);
%     bestcoords4(i,:) = coords4(I(i,1),:);
%     best_rot1(i,1) = EC1_rot(I(i,1),1);
%     best_rot2(i,1) = EC2_rot(I(i,1),1);
%     best_rot3(i,1) = EC3_rot(I(i,1),1);
%     best_rot4(i,1) = EC4_rot(I(i,1),1);
%     bestd(i,1) = d(I(i,1),1);
%     best_placement(i,1) = EC2_placement(I(i,1),1);
% end



%  for g = 2:iter
%       for k = 1:2: pop/2
%           
%           
%         %Selection
%         sel_coeff = d(I(pop/2,1),1);
%         parents1 = zeros(2,3);
%         parents2 = zeros(2,4);
%         parents3 = zeros(2,3);
%         parents4 = zeros(2,3);
%         parents_rot = zeros(2,5);
%         for j=1:2
%             rand_sel = 1 + round(pop.*rand());
%             for i = rand_sel:pop
%                 if d(i,1)<sel_coeff
%                     parents1(j,1) = r1(i,1);
%                     parents1(j,2) = y1(i,1);
%                     parents1(j,3) = theta1(i,1);
%                     parents2(j,1) = r2(i,1);
%                     parents2(j,2) = y2(i,1);
%                     parents2(j,3) = theta2(i,1);
%                     parents2(j,4) = EC2_placement(i,1);
%                     parents3(j,1) = r3(i,1);
%                     parents3(j,2) = y3(i,1);
%                     parents3(j,3) = theta3(i,1);
%                     parents4(j,1) = r4(i,1);
%                     parents4(j,2) = y4(i,1);
%                     parents4(j,3) = theta4(i,1);
%                     parents_rot(j,1) = gb_rot(i,1);
%                     parents_rot(j,2) = EC1_rot(i,1);
%                     parents_rot(j,3) = EC2_rot(i,1);
%                     parents_rot(j,4) = EC3_rot(i,1);
%                     parents_rot(j,5) = EC4_rot(i,1);
%                     
%                     break
%                 end
%             end
%         end
%  
%         %Crossover
%         child1 = zeros(2,3);
%         child2 = zeros(2,4);
%         child3 = zeros(2,3);
%         child4 = zeros(2,3);
%         child_rot = zeros(2,5);
%         alpha = 0.15;
%         beta = 0.85;
%         child1(1,1) = alpha.*parents1(1,1) + (1-alpha).*parents1(2,1);
%         child1(2,1) = beta.*parents1(1,1) + (1-beta).*parents1(2,1);
%         child2(1,1) = alpha.*parents2(1,1) + (1-alpha).*parents2(2,1);
%         child2(2,1) = beta.*parents2(1,1) + (1-beta).*parents2(2,1);
%         child3(1,1) = alpha.*parents3(1,1) + (1-alpha).*parents3(2,1);
%         child3(2,1) = beta.*parents3(1,1) + (1-beta).*parents3(2,1);
%         child4(1,1) = alpha.*parents4(1,1) + (1-alpha).*parents4(2,1);
%         child4(2,1) = beta.*parents4(1,1) + (1-beta).*parents4(2,1);
%         
%         child1(1,2) = alpha.*parents1(1,2) + (1-alpha).*parents1(2,2);
%         child1(2,2) = beta.*parents1(1,2) + (1-beta).*parents1(2,2);
%         child2(1,2) = alpha.*parents2(1,2) + (1-alpha).*parents2(2,2);
%         child2(2,2) = beta.*parents2(1,2) + (1-beta).*parents2(2,2);
%         child3(1,2) = alpha.*parents3(1,2) + (1-alpha).*parents3(2,2);
%         child3(2,2) = beta.*parents3(1,2) + (1-beta).*parents3(2,2);
%         child4(1,2) = alpha.*parents4(1,2) + (1-alpha).*parents4(2,2);
%         child4(2,2) = beta.*parents4(1,2) + (1-beta).*parents4(2,2);
%         
%         child1(1,3) = parents1(1,3);
%         child1(2,3) = parents1(2,3);
%         child2(1,3) = parents2(1,3);
%         child2(2,3) = parents2(2,3);
%         child3(1,3) = parents3(1,3);
%         child3(2,3) = parents3(2,3);
%         child4(1,3) = parents4(1,3);
%         child4(2,3) = parents4(2,3);
%         
%         child2(1,4) = parents2(1,4);
%         child2(2,4) = parents2(1,4);
%         
%         child_rot(1,:) = alpha*parents_rot(1,:) + (1-alpha)*parents_rot(2,:);
%         child_rot(2,:) = beta*parents_rot(1,:) + (1-beta)*parents_rot(2,:);
%         
% 
% 
% %       Mutation
%         for i = 1:2
%             mut_prob = rand();
%             if (mut_prob>=0.75)
%                 child1(i,2) = child1(i,2) + 5;
%                 child2(i,2) = child2(i,2) + 5;
%                 child3(i,2) = child3(i,2) + 5;
%                 child4(i,2) = child4(i,2) + 5;
%             elseif (mut_prob>=0.5) && (mut_prob<0.75)
%                 child1(i,2) = child1(i,2) - 5;
%                 child2(i,2) = child2(i,2) - 5;
%                 child3(i,2) = child3(i,2) - 5;
%                 child4(i,2) = child4(i,2) - 5;
%             end
%             mut_prob = rand();
%             if (mut_prob>=0.75)
%                 child1(i,3) = child1(i,3) + (pi/9)*rand();
%                 child2(i,3) = child2(i,3) + (pi/9)*rand();
%                 child3(i,3) = child3(i,3) + (pi/9)*rand();
%                 child4(i,3) = child1(i,3);
%             elseif (mut_prob>=0.5) && (mut_prob<0.75)
%                 child1(i,3) = child1(i,3) - (pi/9)*rand();
%                 child2(i,3) = child2(i,3) - (pi/9)*rand();
%                 child3(i,3) = child3(i,3) - (pi/9)*rand();
%                 child4(i,3) = child1(i,3);
%             end
%             for j=1:5
%                 mut_prob = rand();
%                 if mut_prob > 0.5 
%                     child_rot(i,j) = child_rot(i,j) + (pi/36);
%                 else
%                     child_rot(i,j) = child_rot(i,j) - (pi/36);
%                 end
%                  
%             end
%         end
%         
%         %Write first half of new generation 
%         r1(k,1) = child1(1,1);
%         y1(k,1) = child1(1,2);
%         theta1(k,1) = child1(1,3);
%         r1(k+1,1) = child1(2,1);
%         y1(k+1,1) = child1(2,2);
%         theta1(k+1,1) = child1(2,3);
%         r2(k,1) = child2(1,1);
%         y2(k,1) = child2(1,2);
%         theta2(k,1) = child2(1,3);
%         EC2_placement(k,1) = child2(1,4);
%         r2(k+1,1) = child2(2,1);
%         y2(k+1,1) = child2(2,2);
%         theta2(k+1,1) = child2(2,3);
%         EC2_placement(k+1,1) = child2(2,4);
%         r3(k,1) = child3(1,1);
%         y3(k,1) = child3(1,2);
%         theta3(k,1) = child3(1,3);
%         r3(k+1,1) = child3(2,1);
%         y3(k+1,1) = child3(2,2);
%         theta3(k+1,1) = child3(2,3);
%         r4(k,1) = child4(1,1);
%         y4(k,1) = child4(1,2);
%         theta4(k,1) = child4(1,3);
%         r4(k+1,1) = child4(2,1);
%         y4(k+1,1) = child4(2,2);
%         theta4(k+1,1) = child4(2,3);
%         gb_rot(k,1) = child_rot(1,1);
%         gb_rot(k+1,1) = child_rot(2,1);
%         EC1_rot(k,1) = child_rot(1,2);
%         EC1_rot(k+1,1) = child_rot(2,2);
%         EC2_rot(k,1) = child_rot(1,3);
%         EC2_rot(k+1,1) = child_rot(2,3);
%         EC3_rot(k,1) = child_rot(1,4);
%         EC3_rot(k+1,1) = child_rot(2,4);
%         EC4_rot(k,1) = child_rot(1,5);
%         EC4_rot(k+1,1) = child_rot(2,5);
%       end
%       
% %Write other half of new population from the best of previous generation     
%       for i=(pop/2+1):pop
%           r1(i,1) = best1(i-(pop/2),1);
%           r2(i,1) = best2(i-(pop/2),1);
%           r3(i,1) = best3(i-(pop/2),1);
%           r4(i,1) = best4(i-(pop/2),1);
%           y1(i,1) = best1(i-(pop/2),2);
%           y2(i,1) = best2(i-(pop/2),2);
%           y3(i,1) = best3(i-(pop/2),2);
%           y4(i,1) = best4(i-(pop/2),2);
%           theta1(i,1) = best1(i-(pop/2),3);
%           theta2(i,1) = best2(i-(pop/2),3);
%           theta3(i,1) = best3(i-(pop/2),3);
%           theta4(i,1) = best4(i-(pop/2),3);
%           EC2_placement(i,1) = best_placement(i-(pop/2),1); 
%           gb_rot(i,1) = bestgb(i-(pop/2),1);
%           EC1_rot(i,1) = best_rot1(i-(pop/2),1);
%           EC2_rot(i,1) = best_rot2(i-(pop/2),1);
%           EC3_rot(i,1) = best_rot3(i-(pop/2),1);
%           EC4_rot(i,1) = best_rot4(i-(pop/2),1);
%           coords1(i,:) = bestcoords1(i-(pop/2),1);
%           coords2(i,:) = bestcoords2(i-(pop/2),1);
%           coords3(i,:) = bestcoords3(i-(pop/2),1);
%           coords4(i,:) = bestcoords4(i-(pop/2),1);
%           d(i,1) = bestd(i-(pop/2),1);
%       end
%      
%       for i=1:pop
%         if i==1
%             gb_mov(i) = gb_rot(i);
%         else
%             gb_mov(i) = gb_rot(i) - gb_rot(i-1);
%         end
%       end
%       
%       for i=1:pop
%         if i==1
%             EC1_mov(i) = theta1(i) + EC1_rot(i);
%         else
%             EC1_mov(i) = theta1(i) + EC1_rot(i) - theta1(i-1) - EC1_rot(i-1);
%         end
%       end
%       
%       for i=1:pop
%         if i==1
%             EC2_mov(i) = theta2(i) + EC2_rot(i);
%         else
%             EC2_mov(i) = theta2(i) + EC2_rot(i) - theta2(i-1) - EC2_rot(i-1);
%         end
%       end
%       
%       for i=1:pop
%         if i==1
%             EC3_mov(i) = theta3(i) + EC3_rot(i);
%         else
%             EC3_mov(i) = theta3(i) + EC3_rot(i) - theta3(i-1) - EC3_rot(i-1);
%         end
%       end
%       
%       for i=1:pop
%         if i==1
%             EC4_mov(i) = theta4(i) + EC4_rot(i);
%         else
%             EC4_mov(i) = theta4(i) + EC4_rot(i) - theta4(i-1) - EC4_rot(i-1);
%         end
%       end




%Test new generation
%         for i=1:pop
% 
%             %Move the components and associated variables to the position of that chromosome
%             %Gearbox
%             gearbox.vertices = Rotation(gearbox.vertices, Anchor, 2, gb_mov(i,1));
%             EM_start = Rotation(EM_start, Anchor, 2, gb_mov(i,1));
%             EM_end = EM_start + EM_vector;
% 
%             coords1(i,:) = [EM_start(1)+ r1(i).*sin(theta1(i)) y1(i) EM_start(3) + r1(i).*cos(theta1(i))];
%             coords2(i,:) = [EM_start(1)+ r2(i).*sin(theta2(i)) y2(i) EM_start(3) + r2(i).*cos(theta2(i))];
%             coords3(i,:) = [EM_start(1)+ r3(i).*sin(theta3(i)) y3(i) EM_start(3) + r3(i).*cos(theta3(i))];
%             coords4(i,:) = [EM_start(1)+ r4(i).*sin(theta4(i)) y4(i) EM_start(3) + r4(i).*cos(theta4(i))];
% 
%             %IGBT
%             EC1.vertices = EC1.vertices - EC1_center + coords1(i,:);
%             EC1_center1 = EC1_center1 - EC1_center + coords1(i,:);
%             EC1_center2 = EC1_center2 - EC1_center + coords1(i,:);
%             EC1_center3 = EC1_center3 - EC1_center + coords1(i,:);
%             EC1_center = coords1(i,:);
%             EC1.vertices = Rotation(EC1.vertices, EC1_center, 2, -EC1_mov(i,1));
%             EC1_center1 = Rotation(EC1_center1, EC1_center, 2, -EC1_mov(i,1));
%             EC1_center2 = Rotation(EC1_center2, EC1_center, 2, -EC1_mov(i,1));
%             EC1_center3 = Rotation(EC1_center3, EC1_center, 2, -EC1_mov(i,1));
% 
%             %Capacitor
%             EC2.vertices = EC2.vertices - EC2_center + coords2(i,:);
%             EC2_center1 = EC2_center1 - EC2_center + coords2(i,:);
%             EC2_center2 = EC2_center2 - EC2_center + coords2(i,:);
%             EC2_center = coords2(i,:);
%             EC2.vertices = Rotation(EC2.vertices, EC2_center, 2, -EC2_mov(i,1));
%             EC2_center1 = Rotation(EC2_center1, EC2_center, 2, -EC2_mov(i,1));
%             EC2_center2 = Rotation(EC2_center2, EC2_center, 2, -EC2_mov(i,1));
% 
%             %EMI-filter
%             EC3.vertices = EC3.vertices - EC3_center + coords3(i,:);
%             EC3_center1 = EC3_center1 - EC3_center + coords3(i,:);
%             EC3_center = coords3(i,:);
%             EC3.vertices = Rotation(EC3.vertices, EC3_center, 2, -EC3_mov(i,1));
%             EC3_center1 = Rotation(EC3_center1, EC3_center, 2, -EC3_mov(i,1));
% 
%             %Control board
%             EC4.vertices = EC4.vertices - EC4_center + coords4(i,:);
%             EC4_center1 = EC4_center1 - EC4_center + coords4(i,:);
%             EC4_center = coords4(i,:);
%             EC4.vertices = Rotation(EC4.vertices, EC4_center, 2, -EC4_mov(i,1));
%             EC4_center1 = Rotation(EC4_center1, EC4_center, 2, -EC4_mov(i,1));
% 
% 
%             %Checking for collisions
% 
%             %Components inside the bauraum
%             IN = inpolyhedron(bauraum,EC1.vertices);
%             in1_bauraum = sum(IN)/length(IN);
%             IN = inpolyhedron(bauraum,EC2.vertices);
%             in2_bauraum = sum(IN)/length(IN);
%             IN = inpolyhedron(bauraum,EC3.vertices);
%             in3_bauraum = sum(IN)/length(IN);
%             IN = inpolyhedron(bauraum,EC4.vertices);
%             in4_bauraum = sum(IN)/length(IN);
% 
%             %Components outside the gearbox/EM
%             IN = inpolyhedron(gearbox,EC1.vertices);
%             out1_gearbox = sum(IN)/length(IN);
%             IN = inpolyhedron(gearbox,EC2.vertices);
%             out2_gearbox = sum(IN)/length(IN);
%             IN = inpolyhedron(gearbox,EC3.vertices);
%             out3_gearbox = sum(IN)/length(IN);
%             IN = inpolyhedron(gearbox,EC4.vertices);
%             out4_gearbox = sum(IN)/length(IN);
% 
% 
%             %Components with themselves
%             IN = inpolyhedron(EC1,EC2.vertices);
%             collision12 = sum(IN)/length(IN);
%             IN = inpolyhedron(EC1,EC3.vertices);
%             collision13 = sum(IN)/length(IN);
%             IN = inpolyhedron(EC1,EC4.vertices);
%             collision14 = sum(IN)/length(IN);
%             IN = inpolyhedron(EC2,EC3.vertices);
%             collision23 = sum(IN)/length(IN);
%             IN = inpolyhedron(EC2,EC4.vertices);
%             collision24 = sum(IN)/length(IN);
%             IN = inpolyhedron(EC3,EC4.vertices);
%             collision34 = sum(IN)/length(IN);
% 
% 
%             d(i,1) = optimization_function(EC1_center, EC1_center1, EC1_center2, EC1_center3, EC2_center1,EC2_center2, EC3_center1, EC3_center2, EC4_center1, EM_start, EM_end, in1_bauraum, in2_bauraum, in3_bauraum, in4_bauraum, out1_gearbox,  out2_gearbox, out3_gearbox, out4_gearbox, collision12, collision13, collision14, collision23, collision24, collision34, EC2_placement(i));
% 
%             hold on     
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
%         end
%         
%         
%         %Return the gearbox to initial angle
%         gearbox.vertices = Rotation(gearbox.vertices, Anchor, 2, -gb_rot(pop,1));
%         EM_start = Rotation(EM_start, Anchor, 2, -gb_rot(pop,1));
%         EM_end = EM_start + EM_vector;
% 
%          %Return components to initial angle
%         EC1.vertices = Rotation(EC1.vertices, EC1_center, 2, theta1(pop,1)+EC1_rot(pop,1));
%         EC1_center1 = Rotation(EC1_center1, EC1_center, 2, EC1_rot(pop,1));
%         EC1_center2 = Rotation(EC1_center2, EC1_center, 2, EC1_rot(pop,1));
%         EC1_center3 = Rotation(EC1_center3, EC1_center, 2, EC1_rot(pop,1));
% 
%         EC2.vertices = Rotation(EC2.vertices, EC2_center, 2, theta2(pop,1)+EC2_rot(pop,1));
%         EC2_center1 = Rotation(EC2_center1, EC2_center, 2, EC2_rot(pop,1));
%         EC2_center2 = Rotation(EC2_center2, EC2_center, 2, EC2_rot(pop,1));
% 
%         EC3.vertices = Rotation(EC3.vertices, EC3_center, 2, theta3(pop,1)+EC3_rot(pop,1));
%         EC3_center1 = Rotation(EC3_center1, EC3_center, 2, EC3_rot(pop,1));
% 
%         EC4.vertices = Rotation(EC4.vertices, EC4_center, 2, theta4(pop,1)+EC4_rot(pop,1));
%         EC4_center1 = Rotation(EC4_center1, EC4_center, 2, EC4_rot(pop,1));
% 
%         %Sort chromosomes based on their fitness values
%         [srt,I]=sort(d);
% 
% %         
%         for i=1:pop/2
%             bestgb(i,1) = gb_rot(I(i,1),1);
%             best1(i,1) = r1(I(i,1),:);
%             best2(i,1) = r2(I(i,1),:);
%             best3(i,1) = r3(I(i,1),:);
%             best4(i,1) = r4(I(i,1),:);
%             best1(i,2) = y1(I(i,1),:);
%             best2(i,2) = y2(I(i,1),:);
%             best3(i,2) = y3(I(i,1),:);
%             best4(i,2) = y4(I(i,1),:);
%             best1(i,3) = theta1(I(i,1),:);
%             best2(i,3) = theta2(I(i,1),:);
%             best3(i,3) = theta3(I(i,1),:);
%             best4(i,3) = theta4(I(i,1),:);
%             best_rot1(i,1) = EC1_rot(I(i,1),1);
%             best_rot2(i,1) = EC2_rot(I(i,1),1);
%             best_rot3(i,1) = EC3_rot(I(i,1),1);
%             best_rot4(i,1) = EC4_rot(I(i,1),1);
%             bestd(i,1) = d(I(i,1),1);
%             best_placement(i,1) = EC2_placement(I(i,1),1);
%         end
%         
%             if (g==iter)
%         %Prints best solution(of last generation)
%             %Gearbox
%             gearbox.vertices = Rotation(gearbox.vertices, Anchor, 2, gb_rot(I(1,1),1));
%             EM_start = Rotation(EM_start, Anchor, 2, gb_rot(I(1,1),1));
%             EM_end = EM_start + EM_vector;
% 
%             %IGBT
%             EC1.vertices = EC1.vertices - EC1_center + coords1(I(1,1),:);
%             EC1_center1 = EC1_center1 - EC1_center + coords1(I(1,1),:);
%             EC1_center2 = EC1_center2 - EC1_center + coords1(I(1,1),:);
%             EC1_center3 = EC1_center3 - EC1_center + coords1(I(1,1),:);
%             EC1_center = coords1(I(1,1),:);
%             EC1.vertices = Rotation(EC1.vertices, EC1_center, 2, -theta1(I(1,1),1) - EC1_rot(I(1,1),1));
%             EC1_center1 = Rotation(EC1_center1, EC1_center, 2, -theta1(I(1,1),1) - EC1_rot(I(1,1),1));
%             EC1_center2 = Rotation(EC1_center2, EC1_center, 2, -theta1(I(1,1),1) - EC1_rot(I(1,1),1));
%             EC1_center3 = Rotation(EC1_center3, EC1_center, 2, -theta1(I(1,1),1) - EC1_rot(I(1,1),1));
% 
%             %Capacitor
%             EC2.vertices = EC2.vertices - EC2_center + coords2(I(1,1),:);
%             EC2_center1 = EC2_center1 - EC2_center + coords2(I(1,1),:);
%             EC2_center2 = EC2_center2 - EC2_center + coords2(I(1,1),:);
%             EC2_center = coords2(I(1,1),:);
%             EC2.vertices = Rotation(EC2.vertices, EC2_center, 2, -theta2(I(1,1),1) - EC2_rot(I(1,1),1));
%             EC2_center1 = Rotation(EC2_center1, EC2_center, 2, -theta2(I(1,1),1) - EC2_rot(I(1,1),1));
%             EC2_center2 = Rotation(EC2_center2, EC2_center, 2, -theta2(I(1,1),1) - EC2_rot(I(1,1),1));
% 
%             %EMI-filter
%             EC3.vertices = EC3.vertices - EC3_center + coords3(I(1,1),:);
%             EC3_center1 = EC3_center1 - EC3_center + coords3(I(1,1),:);
%             EC3_center = coords3(I(1,1),:);
%             EC3.vertices = Rotation(EC3.vertices, EC3_center, 2, -theta3(I(1,1),1) - EC3_rot(I(1,1),1));
%             EC3_center1 = Rotation(EC3_center1, EC3_center, 2, -theta3(I(1,1),1) - EC3_rot(I(1,1),1));
% 
%             %Control board
%             EC4.vertices = EC4.vertices - EC4_center + coords4(I(1,1),:);
%             EC4_center1 = EC4_center1 - EC4_center + coords4(I(1,1),:);
%             EC4_center = coords4(I(1,1),:);
%             EC4.vertices = Rotation(EC4.vertices, EC4_center, 2, -theta4(I(1,1),1) - EC4_rot(I(1,1),1));
%             EC4_center1 = Rotation(EC4_center1, EC4_center, 2, -theta4(I(1,1),1) - EC4_rot(I(1,1),1));
%             
% 
%             hold on    
%             figure(pop+1)
%             BRaum = patch(bauraum,'FaceColor',       [0 0 1.0], ...
%                      'EdgeColor',       'none',        ...
%                      'FaceLighting',    'gouraud',     ...
%                      'AmbientStrength', 0.15);
%             BRaum.FaceVertexAlphaData = 0.2;    % Set constant transparency 
%             BRaum.FaceAlpha = 'flat' ;          % Interpolate to find face transparency
%             hold on
%             view([45 90 135]);
%             xlabel('X') 
%             ylabel('Y')
%             zlabel('Z')
%             
%             
%             hold on     
%             patch(gearbox,'FaceColor',       [0.8 0.8 1], ...
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
%             end
%        fprintf('Selection coefficient of Generation %d: %d \n',g-1,sel_coeff);      
%        fprintf('Fitness of Generation %d: %d \n' ,g,srt(1,1));
% 
%  end

tEnd = cputime - tStart