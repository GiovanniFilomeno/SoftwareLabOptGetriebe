function [r,y,theta,rotations,EC2_placement, Fitness] = NewGeneration(pop,iter,r,y,theta,rotations,EC2_placement, Fitness)
%Sort chromosomes based on their Fitness values
[srt,I]=sort(Fitness);

fprintf('Fitness of Generation %d: %d \n',iter-1,srt(1,1));

%Select the best half of chromosomes from this generation
best_r = zeros(pop/2,4);
best_y = zeros(pop/2,4);
best_theta = zeros(pop/2,4);
best_rotations = zeros(pop/2,5);
best_Fitness = zeros(pop/2,1);
best_placement = zeros(pop/2,1);

for i=1:pop/2
    best_r(i,:) = r(I(i,1),:);
    best_y(i,:) = y(I(i,1),:);
    best_theta(i,:) = theta(I(i,1),:);
    best_rotations(i,:) = rotations(I(i,1),:);
    best_Fitness(i,:) = Fitness(I(i,1),:);
    best_placement(i,:) = EC2_placement(I(i,1),:);
end    
    
for k = 1:2:pop/2

    %Selection
    sel_coeff = Fitness(I(pop/2,1),1);
    parents1 = zeros(2,3);
    parents2 = zeros(2,4);
    parents3 = zeros(2,3);
    parents4 = zeros(2,3);
    parents_rot = zeros(2,5);
    for j=1:2
        rand_sel = 1 + round(pop.*rand());
        for i = rand_sel:pop
            if Fitness(i,1)<sel_coeff
                parents1(j,1) = r(i,1);
                parents1(j,2) = y(i,1);
                parents1(j,3) = theta(i,1);
                parents2(j,1) = r(i,2);
                parents2(j,2) = y(i,2);
                parents2(j,3) = theta(i,2);
                parents2(j,4) = EC2_placement(i,1);
                parents3(j,1) = r(i,3);
                parents3(j,2) = y(i,3);
                parents3(j,3) = theta(i,3);
                parents4(j,1) = r(i,4);
                parents4(j,2) = y(i,4);
                parents4(j,3) = theta(i,4);
                parents_rot(j,1) = rotations(i,1);
                parents_rot(j,2) = rotations(i,2);
                parents_rot(j,3) = rotations(i,3);
                parents_rot(j,4) = rotations(i,4);
                parents_rot(j,5) = rotations(i,5);

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
    child1(1,1) = alpha.*parents1(1,1) + (1-alpha).*parents1(2,1);
    child1(2,1) = (1-alpha).*parents1(1,1) + alpha.*parents1(2,1);
    child2(1,1) = alpha.*parents2(1,1) + (1-alpha).*parents2(2,1);
    child2(2,1) = (1-alpha).*parents2(1,1) + alpha.*parents2(2,1);
    child3(1,1) = alpha.*parents3(1,1) + (1-alpha).*parents3(2,1);
    child3(2,1) = (1-alpha).*parents3(1,1) + alpha.*parents3(2,1);
    child4(1,1) = alpha.*parents4(1,1) + (1-alpha).*parents4(2,1);
    child4(2,1) = (1-alpha).*parents4(1,1) + alpha.*parents4(2,1);

    child1(1,2) = alpha.*parents1(1,2) + (1-alpha).*parents1(2,2);
    child1(2,2) = (1-alpha).*parents1(1,2) + alpha.*parents1(2,2);
    child2(1,2) = alpha.*parents2(1,2) + (1-alpha).*parents2(2,2);
    child2(2,2) = (1-alpha).*parents2(1,2) + alpha.*parents2(2,2);
    child3(1,2) = alpha.*parents3(1,2) + (1-alpha).*parents3(2,2);
    child3(2,2) = (1-alpha).*parents3(1,2) + alpha.*parents3(2,2);
    child4(1,2) = alpha.*parents4(1,2) + (1-alpha).*parents4(2,2);
    child4(2,2) = (1-alpha).*parents4(1,2) + alpha.*parents4(2,2);

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
    child_rot(2,:) = (1-alpha)*parents_rot(1,:) + alpha*parents_rot(2,:);

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
    r(k,1) = child1(1,1);
    y(k,1) = child1(1,2);
    theta(k,1) = child1(1,3);
    r(k+1,1) = child1(2,1);
    y(k+1,1) = child1(2,2);
    theta(k+1,1) = child1(2,3);
    r(k,2) = child2(1,1);
    y(k,2) = child2(1,2);
    theta(k,2) = child2(1,3);
    EC2_placement(k,1) = child2(1,4);
    r(k+1,2) = child2(2,1);
    y(k+1,2) = child2(2,2);
    theta(k+1,2) = child2(2,3);
    EC2_placement(k+1,1) = child2(2,4);
    r(k,3) = child3(1,1);
    y(k,3) = child3(1,2);
    theta(k,3) = child3(1,3);
    r(k+1,3) = child3(2,1);
    y(k+1,3) = child3(2,2);
    theta(k+1,3) = child3(2,3);
    r(k,4) = child4(1,1);
    y(k,4) = child4(1,2);
    theta(k,4) = child4(1,3);
    r(k+1,4) = child4(2,1);
    y(k+1,4) = child4(2,2);
    theta(k+1,4) = child4(2,3);
    rotations(k,1) = child_rot(1,1);
    rotations(k+1,1) = child_rot(2,1);
    rotations(k,2) = child_rot(1,2);
    rotations(k+1,2) = child_rot(2,2);
    rotations(k,3) = child_rot(1,3);
    rotations(k+1,3) = child_rot(2,3);
    rotations(k,4) = child_rot(1,4);
    rotations(k+1,4) = child_rot(2,4);
    rotations(k,5) = child_rot(1,5);
    rotations(k+1,5) = child_rot(2,5);

end

%Write other half of new population from the best of previous generation     
for i=(pop/2+1):pop
    r(i,1) = best_r(i-(pop/2),1);
    r(i,2) = best_r(i-(pop/2),2);
    r(i,3) = best_r(i-(pop/2),3);
    r(i,4) = best_r(i-(pop/2),4);
    y(i,1) = best_y(i-(pop/2),1);
    y(i,2) = best_y(i-(pop/2),2);
    y(i,3) = best_y(i-(pop/2),3);
    y(i,4) = best_y(i-(pop/2),4);
    theta(i,1) = best_theta(i-(pop/2),1);
    theta(i,1) = best_theta(i-(pop/2),2);
    theta(i,1) = best_theta(i-(pop/2),3);
    theta(i,1) = best_theta(i-(pop/2),4);
    EC2_placement(i,1) = best_placement(i-(pop/2),1); 
    rotations(i,1) = best_rotations(i-(pop/2),1);
    rotations(i,2) = best_rotations(i-(pop/2),2);
    rotations(i,3) = best_rotations(i-(pop/2),3);
    rotations(i,4) = best_rotations(i-(pop/2),4);
    rotations(i,5) = best_rotations(i-(pop/2),5);
    Fitness(i,1) = best_Fitness(i-(pop/2),1);
end
    
    
end