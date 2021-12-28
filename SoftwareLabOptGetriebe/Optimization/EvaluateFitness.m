function Fitness = EvaluateFitness(Mean_fitness,pop,bauraum,gearbox,EC1,EC1_width, EC1_center,EC1_SurfCenter,EC2,EC2_center,EC2_SurfCenter,EC3,EC3_center,EC3_SurfCenter,EC4,EC4_center,EC4_SurfCenter,Anchor,EM_start,EM_vector,r,y,theta,Rotations,EC2_placement, Fitness)
% Function that reads the possible rotation of the gearbox, as well as the
% dimensions of the Electric motor and Electric components and then creates
% random sets of values of coordinates and angles between these limits. The
% number of sets created is given by the size of population **pop**. The
% generation of each Electric Component coordinates and angles is connected
% to the values of the other neighboring geometries, as this code application
% assumes the best solution as the one where the ECs are wrapped around the
% Electric motor.
%
% :param double Mean_fitness: Average of the fitness of all the individuals
%                             in the population. Used to decide which
%                             individuals will be considered good enough to generate
%                             splines on them. For the first generation
%                             gets a value of 20000 based on performance of
%                             the first generations in many trials.
% :param int pop: Size of the population.
% :param struct bauraum: struct that represents the Bauraum geometry.
% :param struct gearbox: struct that represents the Gearbox geometry.
% :param struct EC1: struct that represents the IGBT geometry.
% :param double EC1_width: Width of the IGBT. Used for spline creation.
% :param array EC1_center: Vector of coordinates of the body center of the
%                          IGBT.
% :param array EC1_SurfCenter: Vector of coordinates of the centers of the faces of the
%                              IGBT.
% :param struct EC2: struct that represents the Capacitor geometry.
% :param array EC2_center: Vector of coordinates of the body center of the
%                          Capacitor.
% :param array EC2_SurfCenter: Vector of coordinates of the centers of the faces of the
%                              Capacitor.
% :param struct EC3: struct that represents the EMI-filter geometry.
% :param array EC3_center: Vector of coordinates of the body center of the
%                          	EMI-filter.
% :param array EC3_SurfCenter: Vector of coordinates of the centers of the faces of the
%                              EMI-filter.
% :param struct EC4: struct that represents the Control Board geometry.
% :param array EC4_center: Vector of coordinates of the body center of the
%                          Control board.
% :param array EC4_SurfCenter: Vector of coordinates of the centers of the faces of the
%                              Control board.
% :param array Anchor: Coordinates of a point at the center of the gearbox
%                      shaft, which will be used for rotations.
% :param array EM_start: Coordinates of the center of the electric motor.
% :param array EM_vector: Vector pointing the direction of the gearbox.  
% :param array r: Array of coordinates for the distance between each
%                 individual and the center of the electric motor.
% :param array y: Array of coordinates for the Y-position of each
%                 individual.
% :param array theta: Array of angles of the position of each individual's ECs
%                     with respect of the Electric Motor
% :param array Rotations: Array of angles that define the orientation of
%                         each individual's ECs.
% :param array EC2_placement: Array of logical values that define the position of each individual's 
%                             Capacitor with respect of IGBT of said individual.
% :param array Fitness: Current Fitness array, filled with zeros at the
%                       beginning or the fitness of the previous generation. To be overwritten
%                       with the values obtained from this function call.
%
% :return: 
%   **Fitness** : popx1 array with the fitness value of each individual in
%   the population. Initially filled with zeros, as this function is called each generation to evaluate the
%   new individuals this array is overwritten every generation with the
%   fitness values of the new individuals. Fitness calculation is obtained
%   here with inner :func:`OptimizationFunction` and the subsequent
%   computation of the splines via :func:`spline_connection` and :func:`spline_connection_gb`
%
% :rtype: double array
%
% **Example in Code**
%
% .. code-block::
%
%   Fitness = EvaluateFitness(Mean_fitness, iter, g,pop,bauraum,gearbox,EC1,EC1_width, EC1_center,EC1_SurfCenter,EC2,EC2_center,EC2_SurfCenter,EC3,EC3_center,EC3_SurfCenter,EC4,EC4_center,EC4_SurfCenter,Anchor,EM_start,EM_vector,r,y,theta,Rotations,EC2_placement, Fitness);
%

for i=1:pop
    
    %Move the components and associated variables to the position of that chromosome
    %Gearbox
    if i == 1
        gearbox.vertices = Rotation(gearbox.vertices, Anchor, 2, Rotations(i,1));
        EM_start = Rotation(EM_start, Anchor, 2, Rotations(i,1));
    else
        gearbox.vertices = Rotation(gearbox.vertices, Anchor, 2, Rotations(i,1)-Rotations(i-1,1));
        EM_start = Rotation(EM_start, Anchor, 2, Rotations(i,1)-Rotations(i-1,1));
        
    end
    EM_end = EM_start + EM_vector;
    
    coords1(i,:) = [EM_start(1)+ r(i,1).*sin(theta(i,1)) y(i,1) EM_start(3) + r(i,1).*cos(theta(i,1))];
    coords2(i,:) = [EM_start(1)+ r(i,2).*sin(theta(i,2)) y(i,2) EM_start(3) + r(i,2).*cos(theta(i,2))];
    coords3(i,:) = [EM_start(1)+ r(i,3).*sin(theta(i,3)) y(i,3) EM_start(3) + r(i,3).*cos(theta(i,3))];
    coords4(i,:) = [EM_start(1)+ r(i,4).*sin(theta(i,4)) y(i,4) EM_start(3) + r(i,4).*cos(theta(i,4))];
    
    %IGBT
    EC1_mov = zeros(pop,1);
    if i==1
        EC1_mov(i) = theta(i,1) + Rotations(i,2);
    else
        EC1_mov(i) = theta(i,1) + Rotations(i,2) - theta(i-1,1) - Rotations(i-1,2);
    end   
    EC1.vertices = EC1.vertices - EC1_center + coords1(i,:);
    EC1_SurfCenter(1,:) = EC1_SurfCenter(1,:) - EC1_center + coords1(i,:);
    EC1_SurfCenter(2,:) = EC1_SurfCenter(2,:) - EC1_center + coords1(i,:);
    EC1_SurfCenter(3,:) = EC1_SurfCenter(3,:) - EC1_center + coords1(i,:);
    EC1_SurfCenter(4,:) = EC1_SurfCenter(4,:) - EC1_center + coords1(i,:);
    EC1_center = coords1(i,:);
    EC1.vertices = Rotation(EC1.vertices, EC1_center, 2, -EC1_mov(i,1));
    EC1_SurfCenter(1,:) = Rotation(EC1_SurfCenter(1,:), EC1_center, 2, -EC1_mov(i,1));
    EC1_SurfCenter(2,:) = Rotation(EC1_SurfCenter(2,:), EC1_center, 2, -EC1_mov(i,1));
    EC1_SurfCenter(3,:) = Rotation(EC1_SurfCenter(3,:), EC1_center, 2, -EC1_mov(i,1));
    EC1_SurfCenter(4,:) = Rotation(EC1_SurfCenter(4,:), EC1_center, 2, -EC1_mov(i,1));
    
    
    
    %Capacitor
    EC2_mov = zeros(pop,1);
    if i==1
        EC2_mov(i) = theta(i,2) + Rotations(i,3);
    else
        EC2_mov(i) = theta(i,2) + Rotations(i,3) - theta(i-1,2) - Rotations(i-1,3);
    end   
    EC2.vertices = EC2.vertices - EC2_center + coords2(i,:);
    EC2_SurfCenter(1,:) = EC2_SurfCenter(1,:) - EC2_center + coords2(i,:);
    EC2_SurfCenter(2,:) = EC2_SurfCenter(2,:) - EC2_center + coords2(i,:);
    EC2_SurfCenter(3,:) = EC2_SurfCenter(3,:) - EC2_center + coords2(i,:);
    EC2_SurfCenter(4,:) = EC2_SurfCenter(4,:) - EC2_center + coords2(i,:);
    EC2_center = coords2(i,:);
    EC2.vertices = Rotation(EC2.vertices, EC2_center, 2, -EC2_mov(i,1));
    EC2_SurfCenter(1,:) = Rotation(EC2_SurfCenter(1,:), EC2_center, 2, -EC2_mov(i,1));
    EC2_SurfCenter(2,:) = Rotation(EC2_SurfCenter(2,:), EC2_center, 2, -EC2_mov(i,1));
    EC2_SurfCenter(3,:) = Rotation(EC2_SurfCenter(3,:), EC2_center, 2, -EC2_mov(i,1));
    EC2_SurfCenter(4,:) = Rotation(EC2_SurfCenter(4,:), EC2_center, 2, -EC2_mov(i,1));
    
    %EMI-filter
    EC3_mov = zeros(pop,1);
    if i==1
        EC3_mov(i) = theta(i,3) + Rotations(i,4);
    else
        EC3_mov(i) = theta(i,3) + Rotations(i,4) - theta(i-1,3) - Rotations(i-1,4);
    end   
    EC3.vertices = EC3.vertices - EC3_center + coords3(i,:);
    EC3_SurfCenter(1,:) = EC3_SurfCenter(1,:) - EC3_center + coords3(i,:);
    EC3_SurfCenter(2,:) = EC3_SurfCenter(2,:) - EC3_center + coords3(i,:);
    EC3_SurfCenter(3,:) = EC3_SurfCenter(3,:) - EC3_center + coords3(i,:);
    EC3_SurfCenter(4,:) = EC3_SurfCenter(4,:) - EC3_center + coords3(i,:);
    EC3_center = coords3(i,:);
    EC3.vertices = Rotation(EC3.vertices, EC3_center, 2, -EC3_mov(i,1));
    EC3_SurfCenter(1,:) = Rotation(EC3_SurfCenter(1,:), EC3_center, 2, -EC3_mov(i,1));
    EC3_SurfCenter(2,:) = Rotation(EC3_SurfCenter(2,:), EC3_center, 2, -EC3_mov(i,1));
    EC3_SurfCenter(3,:) = Rotation(EC3_SurfCenter(3,:), EC3_center, 2, -EC3_mov(i,1));
    EC3_SurfCenter(4,:) = Rotation(EC3_SurfCenter(4,:), EC3_center, 2, -EC3_mov(i,1));
    
    %Control board
    EC4_mov = zeros(pop,1);
    if i==1
        EC4_mov(i) = theta(i,4) + Rotations(i,5);
    else
        EC4_mov(i) = theta(i,4) + Rotations(i,5) - theta(i-1,4) - Rotations(i-1,5);
    end   
    EC4.vertices = EC4.vertices - EC4_center + coords4(i,:);
    EC4_SurfCenter(1,:) = EC4_SurfCenter(1,:) - EC4_center + coords4(i,:);
    EC4_SurfCenter(2,:) = EC4_SurfCenter(2,:) - EC4_center + coords4(i,:);
    EC4_SurfCenter(3,:) = EC4_SurfCenter(3,:) - EC4_center + coords4(i,:);
    EC4_SurfCenter(4,:) = EC4_SurfCenter(4,:) - EC4_center + coords4(i,:);
    EC4_center = coords4(i,:);
    EC4.vertices = Rotation(EC4.vertices, EC4_center, 2, -EC4_mov(i,1));
    EC4_SurfCenter(1,:) = Rotation(EC4_SurfCenter(1,:), EC4_center, 2, -EC4_mov(i,1));
    EC4_SurfCenter(2,:) = Rotation(EC4_SurfCenter(2,:), EC4_center, 2, -EC4_mov(i,1));
    EC4_SurfCenter(3,:) = Rotation(EC4_SurfCenter(3,:), EC4_center, 2, -EC4_mov(i,1));
    EC4_SurfCenter(4,:) = Rotation(EC4_SurfCenter(4,:), EC4_center, 2, -EC4_mov(i,1));
    
    cable_radius = 2.3;
    radius_of_curvature = 13.8; 
    EM_radius = 117;
    

    %Checking for collisions
    
    %Gearbox and bauraum
    IN = inpolyhedron(bauraum,gearbox.vertices);
    bauraum_gearbox = sum(IN)/length(IN);
    
    
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
    IN = inpolyhedron(gearbox,EC1.vertices);
    out1_gearbox = sum(IN)/length(IN);
    IN = inpolyhedron(gearbox,EC2.vertices);
    out2_gearbox = sum(IN)/length(IN);
    IN = inpolyhedron(gearbox,EC3.vertices);
    out3_gearbox = sum(IN)/length(IN);
    IN = inpolyhedron(gearbox,EC4.vertices);
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
    
    %Calculation of the distance and application of penalties from
    %collisions to generate initial fitness.
    Fitness(i,1) = OptimizationFunction(bauraum_gearbox, EC1_center, EC1_SurfCenter, EC2_SurfCenter, EC3_SurfCenter, EC4_SurfCenter, EM_start, EM_end, in1_bauraum, in2_bauraum, in3_bauraum, in4_bauraum, out1_gearbox,  out2_gearbox, out3_gearbox, out4_gearbox, collision12, collision13, collision14, collision23, collision24, collision34, EC2_placement(i,1));  

    Clearance = 5;
    %If fitness of an individual is below average/5, said individual
    %proceeds to spline generation where fitness can be applied a penalty
    %if spline generation is unsuccessful.
     if Fitness(i,1) < Mean_fitness/5
        if EC2_placement(i,1) == 1   %Counterclockwise        
            vector1 = EC1_SurfCenter(2,:) - EC1_SurfCenter(1,:);
            vector2 = EC2_SurfCenter(1,:) - EC2_SurfCenter(2,:);
            pt2 = EC1_SurfCenter(1,:) - (Clearance.*vector1)/norm(vector1);
            pt3 = EC2_SurfCenter(2,:) - (Clearance.*vector2)/norm(vector2);
            [splinelength,max_curvature, curvature, collisionwith1, collisionwith2, Penalty_spline, fn] = spline_connection (EC1_SurfCenter(1,:), pt2, pt3, EC2_SurfCenter(2,:),EC1_center, EC2_center, EC1, EC2, bauraum, cable_radius, radius_of_curvature);
            Fitness(i,1) = Fitness(i,1) + Penalty_spline;
            if Penalty_spline == 0
                Fitness(i,1) = Fitness(i,1) + splinelength/10;
            end
            vector1 = EC2_SurfCenter(2,:) - EC2_SurfCenter(1,:);
            vector2 = EC3_SurfCenter(1,:) - EC3_SurfCenter(2,:);
            pt2 = EC2_SurfCenter(1,:) - (Clearance.*vector1)/norm(vector1);
            pt3 = EC3_SurfCenter(2,:) - (Clearance.*vector2)/norm(vector2);
            [splinelength,max_curvature, curvature, collisionwith1, collisionwith2, Penalty_spline, fn] = spline_connection (EC2_SurfCenter(1,:), pt2, pt3, EC3_SurfCenter(2,:),EC2_center, EC3_center, EC2, EC3, bauraum, cable_radius, radius_of_curvature);
            Fitness(i,1) = Fitness(i,1) + Penalty_spline;
            if Penalty_spline == 0
                Fitness(i,1) = Fitness(i,1) + splinelength/10;
            end
            vector1 = EC1_SurfCenter(4,:) - EC1_SurfCenter(3,:);
            vector2 = EC4_SurfCenter(3,:) - EC4_SurfCenter(4,:);
            pt2 = EC1_SurfCenter(3,:) - (Clearance.*vector1)/norm(vector1);
            pt3 = EC4_SurfCenter(4,:) - (Clearance.*vector2)/norm(vector2);
            [splinelength,max_curvature, curvature, collisionwith1, collisionwith2, Penalty_spline, fn] = spline_connection (EC1_SurfCenter(3,:), pt2, pt3, EC4_SurfCenter(4,:),EC1_center, EC4_center, EC1, EC4, bauraum, cable_radius, radius_of_curvature);
            Fitness(i,1) = Fitness(i,1) + Penalty_spline;
            if Penalty_spline == 0
                Fitness(i,1) = Fitness(i,1) + splinelength/10;
            end

            EM_connection = EM_end;
            clearance_mid = 20;
            clearance_posy = 30;
            clearance_negy = 10;
            rad_gb = EM_radius;
            vector1 = EC1_SurfCenter(2,:) - EC1_SurfCenter(1,:);

            %pt3 mid spline
            vector_pt1 = [EC1_SurfCenter(2,1)+(-clearance_mid) EM_connection(2)+(-clearance_mid) EC1_SurfCenter(2,3)];
            vector_pt2 = EM_connection+[0 -clearance_mid 0]; 
            vector_mid = vector_pt1 - vector_pt2;
            vector_mid = vector_mid / norm(vector_mid);           
            pt3 = EM_connection+[0 -clearance_mid 0]+ rad_gb*vector_mid;
            pt6 = EM_connection + [cable_radius 0 0];
            pt11 = EM_connection + [-cable_radius 0 0];
            pt7 = EM_connection+[0 -clearance_mid 0] + [cable_radius 0 0];
            pt12 = EM_connection+[0 -clearance_mid 0] + [-cable_radius 0 0];
            pt10 = EC1_SurfCenter(2,:) + [0 EC1_width/4 0];
            pt15 = EC1_SurfCenter(2,:) + [0 -EC1_width/4 0];
            pt9 = pt10 + (clearance_posy.*vector1)/norm(vector1);
            pt14 = pt15 + (clearance_negy.*vector1)/norm(vector1);
            %pt8 pos y spline
            vector_pt3 = [pt9(1) EM_connection(2)+(-clearance_mid) EC1_SurfCenter(2,3)];
            vector_pt4 = pt6 +[0 -clearance_mid 0]; 
            vector_posy = vector_pt3 - vector_pt4;
            vector_posy = vector_posy / norm(vector_posy);
            pt8 = EM_connection+[0 -clearance_posy 0]+ rad_gb*vector_posy;
            %pt13 neg y spline
            vector_pt5 = [pt14(1) EM_connection(2)+(-clearance_mid) EC1_SurfCenter(2,3)];
            vector_pt6 = pt11 +[0 -clearance_mid 0]; 
            vector_negy = vector_pt5 - vector_pt6;
            vector_negy = vector_negy / norm(vector_negy);
            pt13 = EM_connection+[0 -clearance_negy 0]+ rad_gb*vector_negy;
            [splinelength,max_curvature, curvature, collisionwith1, collisionwithgb,fn_mid,fn_posy,fn_negy] = spline_connection_gb (EM_connection, EM_connection+[0 -clearance_mid 0], pt3, EC1_SurfCenter(2,:)+(clearance_mid.*vector1)/norm(vector1), EC1_SurfCenter(2,:), pt6, pt7, pt8, pt9, pt10, pt11, pt12, pt13, pt14, pt15,EC1_center, EM_start, EM_end, gearbox, EC1,  cable_radius, radius_of_curvature);
            Fitness(i,1) = Fitness(i,1) + splinelength/10;
        else                    %Clockwise

            vector1 = EC1_SurfCenter(1,:) - EC1_SurfCenter(2,:);
            vector2 = EC2_SurfCenter(2,:) - EC2_SurfCenter(1,:);
            pt2 = EC1_SurfCenter(2,:) - (Clearance.*vector1)/norm(vector1);
            pt3 = EC2_SurfCenter(1,:) - (Clearance.*vector2)/norm(vector2);        
            [splinelength,max_curvature, curvature, collisionwith1, collisionwith2, Penalty_spline, fn] = spline_connection (EC1_SurfCenter(2,:), pt2, pt3, EC2_SurfCenter(1,:),EC1_center, EC2_center, EC1, EC2, bauraum, cable_radius, radius_of_curvature);
            Fitness(i,1) = Fitness(i,1) + Penalty_spline;
            if Penalty_spline == 0
                Fitness(i,1) = Fitness(i,1) + splinelength/10;
            end
            vector1 = EC2_SurfCenter(1,:) - EC2_SurfCenter(2,:);
            vector2 = EC3_SurfCenter(2,:) - EC3_SurfCenter(1,:);
            pt2 = EC2_SurfCenter(2,:) - (Clearance.*vector1)/norm(vector1);
            pt3 = EC3_SurfCenter(1,:) - (Clearance.*vector2)/norm(vector2);        
            [splinelength,max_curvature, curvature, collisionwith1, collisionwith2, Penalty_spline, fn] = spline_connection (EC2_SurfCenter(2,:), pt2, pt3, EC3_SurfCenter(1,:),EC2_center, EC3_center, EC2, EC3, bauraum, cable_radius, radius_of_curvature);
            Fitness(i,1) = Fitness(i,1) + Penalty_spline;
            if Penalty_spline == 0
                Fitness(i,1) = Fitness(i,1) + splinelength/10;
            end
            vector1 = EC1_SurfCenter(4,:) - EC1_SurfCenter(3,:);
            vector2 = EC4_SurfCenter(3,:) - EC4_SurfCenter(4,:);
            pt2 = EC1_SurfCenter(3,:) - (Clearance.*vector1)/norm(vector1);
            pt3 = EC4_SurfCenter(4,:) - (Clearance.*vector2)/norm(vector2);        
            [splinelength,max_curvature, curvature, collisionwith1, collisionwith2, Penalty_spline, fn] = spline_connection (EC1_SurfCenter(3,:), pt2, pt3, EC4_SurfCenter(4,:),EC1_center, EC4_center, EC1, EC4, bauraum, cable_radius, radius_of_curvature);
            Fitness(i,1) = Fitness(i,1) + Penalty_spline;
            if Penalty_spline == 0
                Fitness(i,1) = Fitness(i,1) + splinelength/10;
            end

            EM_connection = EM_end;
            clearance_mid = 20;
            clearance_posy = 30;
            clearance_negy = 10;
            rad_gb = EM_radius;
            vector1 = EC1_SurfCenter(1,:) - EC1_SurfCenter(2,:);

            %pt3 mid spline
            vector_pt1 = [EC1_SurfCenter(1,1)+(clearance_mid) EM_connection(2)+(-clearance_mid) EC1_SurfCenter(1,3)];
            vector_pt2 = EM_connection+[0 -clearance_mid 0]; 
            vector_mid = vector_pt1 - vector_pt2;
            vector_mid = vector_mid / norm(vector_mid);           
            pt3 = EM_connection+[0 -clearance_mid 0]+ rad_gb*vector_mid;
            pt6 = EM_connection + [cable_radius 0 0];
            pt11 = EM_connection + [-cable_radius 0 0];
            pt7 = EM_connection+[0 -clearance_mid 0] + [cable_radius 0 0];
            pt12 = EM_connection+[0 -clearance_mid 0] + [-cable_radius 0 0];
            pt10 = EC1_SurfCenter(1,:) + [0 EC1_width/4 0];
            pt15 = EC1_SurfCenter(1,:) + [0 -EC1_width/4 0];
            pt9 = pt10 + (clearance_posy.*vector1)/norm(vector1);
            pt14 = pt15 + (clearance_negy.*vector1)/norm(vector1);
            %pt8 pos y spline
            vector_pt3 = [pt9(1) EM_connection(2)+(-clearance_mid) EC1_SurfCenter(1,3)];
            vector_pt4 = pt6 +[0 -clearance_mid 0]; 
            vector_posy = vector_pt3 - vector_pt4;
            vector_posy = vector_posy / norm(vector_posy);
            pt8 = EM_connection+[0 -clearance_posy 0]+ rad_gb*vector_posy;
            %pt13 neg y spline
            vector_pt5 = [pt14(1) EM_connection(2)+(-clearance_mid) EC1_SurfCenter(1,3)];
            vector_pt6 = pt11 +[0 -clearance_mid 0]; 
            vector_negy = vector_pt5 - vector_pt6;
            vector_negy = vector_negy / norm(vector_negy);
            pt13 = EM_connection+[0 -clearance_negy 0]+ rad_gb*vector_negy;
            [splinelength,max_curvature, curvature, collisionwith1, collisionwithgb,fn_mid,fn_posy,fn_negy] = spline_connection_gb (EM_connection, EM_connection+[0 -clearance_mid 0], pt3, EC1_SurfCenter(1,:)+(clearance_mid.*vector1)/norm(vector1), EC1_SurfCenter(1,:), pt6, pt7, pt8, pt9, pt10, pt11, pt12, pt13, pt14, pt15,EC1_center, EM_start, EM_end, gearbox, EC1, cable_radius, radius_of_curvature);
            Fitness(i,1) = Fitness(i,1) + splinelength/10;
        end
     end

    
    
end


end