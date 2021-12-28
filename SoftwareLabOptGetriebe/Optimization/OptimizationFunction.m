function fitness = OptimizationFunction(bauraum_gearbox, EC1_center, EC1_SurfCenter, EC2_SurfCenter, EC3_SurfCenter, EC4_SurfCenter, EM_start, EM_end, in1_bauraum, in2_bauraum, in3_bauraum, in4_bauraum, out1_gearbox,  out2_gearbox, out3_gearbox, out4_gearbox, collision12, collision13, collision14, collision23, collision24, collision34, EC2_placement) 
% Function that calculates the distance between IGBT and Electric Motor,
% as well as the distance between components that are connected using :func:`point_to_line`. 
% Then uses the collisions calculated with :func:`inpolyhedron` for the current individual
% being evaluated in :func:`EvaluateFitness` to apply penalties in case of
% collisions and obtain the fitness value for said individual
%
% :param double bauraum_gearbox: Interference between gearbox and bauraum.
% :param array EC1_center: Vector of coordinates of the body center of the
%                          IGBT.
% :param array EC1_SurfCenter: Vector of coordinates of the centers of the faces of the
%                              IGBT.
% :param array EC2_SurfCenter: Vector of coordinates of the centers of the faces of the
%                              Capacitor.
% :param array EC3_SurfCenter: Vector of coordinates of the centers of the faces of the
%                              EMI-filter.
% :param array EC4_SurfCenter: Vector of coordinates of the centers of the faces of the
%                              Control board.
% :param array EM_start: Coordinates of the center of the electric motor on one end.
% :param array EM_end: Coordinates of the end of the electric motor opposite to EM_start.
% :param double in1_bauraum: Checks if IGBT is fully inside bauraum.
% :param double in2_bauraum: Checks if Capacitor is fully inside bauraum.
% :param double in3_bauraum: Checks if EMI-filter is fully inside bauraum.
% :param double in4_bauraum: Checks if Control board is fully inside bauraum.
% :param double out1_gearbox: Checks if IGBT is colliding with gearbox.
% :param double out2_gearbox: Checks if Capacitor is colliding with gearbox.
% :param double out3_gearbox: Checks if EMI-filter is colliding with gearbox.
% :param double out4_gearbox: Checks if Control board is colliding with gearbox.
% :param double collision12: Checks if IGBT is colliding with Capacitor.
% :param double collision13: Checks if IGBT is colliding with EMI-filter.
% :param double collision14: Checks if IGBT is colliding with Control board.
% :param double collision23: Checks if Capacitor is colliding with EMI-filter.
% :param double collision24: Checks if Capacitor is colliding with Control board.
% :param double collision34: Checks if EMI-filter is colliding with Control board.
% :param array EC2_placement: Array of logical values that define the position of each individual's 
%                             Capacitor with respect of IGBT of said individual.
%
% :return: 
%   **Fitness** : Double value with the fitness of individual in
%   the population that is currently being evaluated in :func:`EvaluateFitness`. Depending on the obtained value from this function fitness can still be modified via :func:`spline_connection` and :func:`spline_connection_gb`
%
% :rtype: double
%
% **Example in Code**
%
% .. code-block::
%
%   Fitness(i,1) = OptimizationFunction(bauraum_gearbox, EC1_center, EC1_SurfCenter, EC2_SurfCenter, EC3_SurfCenter, EC4_SurfCenter, EM_start, EM_end, in1_bauraum, in2_bauraum, in3_bauraum, in4_bauraum, out1_gearbox,  out2_gearbox, out3_gearbox, out4_gearbox, collision12, collision13, collision14, collision23, collision24, collision34, EC2_placement(i,1));  
%    



    %Penalty factors depending on the situation
    Pen_gearbox = 5000*5;
    Pen_bauraum = 5000*5;
    Pen_collision = 5000*5;
    Pen_distance = 1000;
    min_distance = 30;
    
    %Checks if gearbox rotation caused it to go out of the bauraum
    if bauraum_gearbox ~= 0.9375
        Pen_baugb = 1000;
    else
        Pen_baugb = 0;
    end
        
%% Calculating distances, applying penalties and computing fitness

    %IGBT to Electric motor
    if EC1_center(1,2) > EM_start(1,2) || EC1_center(1,2) < EM_end(1,2)
         distance = point_to_line(EC1_center, EM_start,EM_end);
         fitness1 = distance + (1-in1_bauraum)*Pen_bauraum +(out1_gearbox)*Pen_gearbox;
     else
         distance = point_to_line(EC1_center, EM_start,EM_end);
         fitness1 = distance - 135.85 + (1-in1_bauraum)*Pen_bauraum +(out1_gearbox)*Pen_gearbox;
    end
 
   
    %Counterclockwise placement of Condensator
    if EC2_placement == 1       
        %IGBT to condensator
            distance = norm(EC1_SurfCenter(1,:) - EC2_SurfCenter(2,:));
            if distance < min_distance
                distance = distance+Pen_distance;
            end
            fitness2 = distance + (1-in2_bauraum)*Pen_bauraum +(out2_gearbox)*Pen_gearbox;
        %Condensator to EMI-filter    
            distance = norm(EC2_SurfCenter(1,:) - EC3_SurfCenter(2,:));
            if distance < min_distance
                distance = distance+Pen_distance;
            end
            fitness3 = distance + (1-in3_bauraum)*Pen_bauraum +(out3_gearbox)*Pen_gearbox;
            
    %Clockwise placement of Condensator
    else
        %IGBT to condensator
            distance = norm(EC1_SurfCenter(2,:) - EC2_SurfCenter(1,:));
            if distance < min_distance
                distance = distance+Pen_distance;
            end
            fitness2 = distance + (1-in2_bauraum)*Pen_bauraum +(out2_gearbox)*Pen_gearbox;
        %Condensator to EMI-filter    
            distance = norm(EC2_SurfCenter(2,:) - EC3_SurfCenter(1,:));
            if distance < min_distance
                distance = distance+Pen_distance;
            end
            fitness3 = distance + (1-in3_bauraum)*Pen_bauraum +(out3_gearbox)*Pen_gearbox;
    end
        
    %IGBT to Control board    
    distance = norm(EC1_SurfCenter(3,:) - EC4_SurfCenter(4,:));
    if distance < min_distance
                distance = distance+Pen_distance;
    end
    fitness4 = distance + (1-in4_bauraum)*Pen_bauraum +(out4_gearbox)*Pen_gearbox;
    
    
    fitness = Pen_baugb + fitness1 + fitness2 + fitness3 + fitness4 + Pen_collision*collision12 + Pen_collision*collision13 + Pen_collision*collision14 + Pen_collision*collision23 + Pen_collision*collision24 + Pen_collision*collision34;

 end