function fitness = OptimizationFunction(bauraum_gearbox, EC1_center, EC1_SurfCenter, EC2_SurfCenter, EC3_SurfCenter, EC4_SurfCenter, EM_start, EM_end, in1_bauraum, in2_bauraum, in3_bauraum, in4_bauraum, out1_gearbox,  out2_gearbox, out3_gearbox, out4_gearbox, collision12, collision13, collision14, collision23, collision24, collision34, EC2_placement) 
    Pen_gearbox = 5000*5;
    Pen_bauraum = 5000*5;
    Pen_collision = 5000*5;
    Pen_distance = 1000;
    min_distance = 30;
    
%     EC1_SurfCenter(1,1) = EC1_SurfCenter(1,1) + 10;
%     EC1_SurfCenter(2,1) = EC1_SurfCenter(2,1) - 10;
%     EC2_SurfCenter(1,1) = EC2_SurfCenter(1,1) + 10;
%     EC2_SurfCenter(2,1) = EC2_SurfCenter(2,1) - 10;
%     EC3_SurfCenter(1,1) = EC3_SurfCenter(1,1) + 10;
%     EC3_SurfCenter(2,1) = EC3_SurfCenter(2,1) - 10;

    if bauraum_gearbox ~= 0.9375
        Pen_baugb = 1000;
    else
        Pen_baugb = 0;
    end
        

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
                distance = distance*Pen_distance;
            end
            fitness2 = distance + (1-in2_bauraum)*Pen_bauraum +(out2_gearbox)*Pen_gearbox;
        %Condensator to EMI-filter    
            distance = norm(EC2_SurfCenter(1,:) - EC3_SurfCenter(2,:));
            if distance < min_distance
                distance = distance*Pen_distance;
            end
            fitness3 = distance + (1-in3_bauraum)*Pen_bauraum +(out3_gearbox)*Pen_gearbox;
            
    %Clockwise placement of Condensator
    else
        %IGBT to condensator
            distance = norm(EC1_SurfCenter(2,:) - EC2_SurfCenter(1,:));
            if distance < min_distance
                distance = distance*Pen_distance;
            end
            fitness2 = distance + (1-in2_bauraum)*Pen_bauraum +(out2_gearbox)*Pen_gearbox;
        %Condensator to EMI-filter    
            distance = norm(EC2_SurfCenter(2,:) - EC3_SurfCenter(1,:));
            if distance < min_distance
                distance = distance*Pen_distance;
            end
            fitness3 = distance + (1-in3_bauraum)*Pen_bauraum +(out3_gearbox)*Pen_gearbox;
    end
        
    %IGBT to Control board    
    distance = norm(EC1_SurfCenter(3,:) - EC4_SurfCenter(4,:));
    if distance < min_distance
                distance = distance*Pen_distance;
    end
    fitness4 = distance + (1-in4_bauraum)*Pen_bauraum +(out4_gearbox)*Pen_gearbox;
    
    
    fitness = Pen_baugb + fitness1 + fitness2 + fitness3 + fitness4 + Pen_collision*collision12 + Pen_collision*collision13 + Pen_collision*collision14 + Pen_collision*collision23 + Pen_collision*collision24 + Pen_collision*collision34;

 end