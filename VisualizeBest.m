function VisualizeBest(iter,bauraum, gearbox,Anchor, EC1,EC1_width, EC1_center, EC1_SurfCenter, EC2, EC2_center, EC2_SurfCenter, EC3, EC3_center, EC3_SurfCenter, EC4, EC4_center, EC4_SurfCenter, EM_start,EM_vector, Fitness,r,y,theta,Rotations,EC2_placement)

%Move the components and associated variables to the position of that chromosome
%Gearbox
[srt,I] = sort(Fitness);

gearbox.vertices = Rotation(gearbox.vertices, Anchor, 2, Rotations(I(1,1),1));
EM_start = Rotation(EM_start, Anchor, 2, Rotations(I(1,1),1));
EM_end = EM_start + EM_vector;

coords1(I(1,1),:) = [EM_start(1)+ r(I(1,1),1).*sin(theta(I(1,1),1)) y(I(1,1),1) EM_start(3) + r(I(1,1),1).*cos(theta(I(1,1),1))];
coords2(I(1,1),:) = [EM_start(1)+ r(I(1,1),2).*sin(theta(I(1,1),2)) y(I(1,1),2) EM_start(3) + r(I(1,1),2).*cos(theta(I(1,1),2))];
coords3(I(1,1),:) = [EM_start(1)+ r(I(1,1),3).*sin(theta(I(1,1),3)) y(I(1,1),3) EM_start(3) + r(I(1,1),3).*cos(theta(I(1,1),3))];
coords4(I(1,1),:) = [EM_start(1)+ r(I(1,1),4).*sin(theta(I(1,1),4)) y(I(1,1),4) EM_start(3) + r(I(1,1),4).*cos(theta(I(1,1),4))];

%IGBT

EC1_mov = theta(I(1,1),1) + Rotations(I(1,1),2);
EC1.vertices = EC1.vertices - EC1_center + coords1(I(1,1),:);
EC1_SurfCenter(1,:) = EC1_SurfCenter(1,:) - EC1_center + coords1(I(1,1),:);
EC1_SurfCenter(2,:) = EC1_SurfCenter(2,:) - EC1_center + coords1(I(1,1),:);
EC1_SurfCenter(3,:) = EC1_SurfCenter(3,:) - EC1_center + coords1(I(1,1),:);
EC1_SurfCenter(4,:) = EC1_SurfCenter(4,:) - EC1_center + coords1(I(1,1),:);
EC1_center = coords1(I(1,1),:);
EC1.vertices = Rotation(EC1.vertices, EC1_center, 2, -EC1_mov);
EC1_SurfCenter(1,:) = Rotation(EC1_SurfCenter(1,:), EC1_center, 2, -EC1_mov);
EC1_SurfCenter(2,:) = Rotation(EC1_SurfCenter(2,:), EC1_center, 2, -EC1_mov);
EC1_SurfCenter(3,:) = Rotation(EC1_SurfCenter(3,:), EC1_center, 2, -EC1_mov);
EC1_SurfCenter(4,:) = Rotation(EC1_SurfCenter(4,:), EC1_center, 2, -EC1_mov);


%Capacitor

EC2_mov = theta(I(1,1),2) + Rotations(I(1,1),3);
EC2.vertices = EC2.vertices - EC2_center + coords2(I(1,1),:);
EC2_SurfCenter(1,:) = EC2_SurfCenter(1,:) - EC2_center + coords2(I(1,1),:);
EC2_SurfCenter(2,:) = EC2_SurfCenter(2,:) - EC2_center + coords2(I(1,1),:);
EC2_SurfCenter(3,:) = EC2_SurfCenter(3,:) - EC2_center + coords2(I(1,1),:);
EC2_SurfCenter(4,:) = EC2_SurfCenter(4,:) - EC2_center + coords2(I(1,1),:);
EC2_center = coords2(I(1,1),:);
EC2.vertices = Rotation(EC2.vertices, EC2_center, 2, -EC2_mov);
EC2_SurfCenter(1,:) = Rotation(EC2_SurfCenter(1,:), EC2_center, 2, -EC2_mov);
EC2_SurfCenter(2,:) = Rotation(EC2_SurfCenter(2,:), EC2_center, 2, -EC2_mov);
EC2_SurfCenter(3,:) = Rotation(EC2_SurfCenter(3,:), EC2_center, 2, -EC2_mov);
EC2_SurfCenter(4,:) = Rotation(EC2_SurfCenter(4,:), EC2_center, 2, -EC2_mov);


%EMI-filter

EC3_mov = theta(I(1,1),3) + Rotations(I(1,1),4);
EC3.vertices = EC3.vertices - EC3_center + coords3(I(1,1),:);
EC3_SurfCenter(1,:) = EC3_SurfCenter(1,:) - EC3_center + coords3(I(1,1),:);
EC3_SurfCenter(2,:) = EC3_SurfCenter(2,:) - EC3_center + coords3(I(1,1),:);
EC3_SurfCenter(3,:) = EC3_SurfCenter(3,:) - EC3_center + coords3(I(1,1),:);
EC3_SurfCenter(4,:) = EC3_SurfCenter(4,:) - EC3_center + coords3(I(1,1),:);
EC3_center = coords3(I(1,1),:);
EC3.vertices = Rotation(EC3.vertices, EC3_center, 2, -EC3_mov);
EC3_SurfCenter(1,:) = Rotation(EC3_SurfCenter(1,:), EC3_center, 2, -EC3_mov);
EC3_SurfCenter(2,:) = Rotation(EC3_SurfCenter(2,:), EC3_center, 2, -EC3_mov);
EC3_SurfCenter(3,:) = Rotation(EC3_SurfCenter(3,:), EC3_center, 2, -EC3_mov);
EC3_SurfCenter(4,:) = Rotation(EC3_SurfCenter(4,:), EC3_center, 2, -EC3_mov);


%Control board
EC4_mov = theta(I(1,1),4) + Rotations(I(1,1),5);
EC4.vertices = EC4.vertices - EC4_center + coords4(I(1,1),:);
EC4_SurfCenter(1,:) = EC4_SurfCenter(1,:) - EC4_center + coords4(I(1,1),:);
EC4_SurfCenter(2,:) = EC4_SurfCenter(2,:) - EC4_center + coords4(I(1,1),:);
EC4_SurfCenter(3,:) = EC4_SurfCenter(3,:) - EC4_center + coords4(I(1,1),:);
EC4_SurfCenter(4,:) = EC4_SurfCenter(4,:) - EC4_center + coords4(I(1,1),:);
EC4_center = coords4(I(1,1),:);
EC4.vertices = Rotation(EC4.vertices, EC4_center, 2, -EC4_mov);
EC4_SurfCenter(1,:) = Rotation(EC4_SurfCenter(1,:), EC4_center, 2, -EC4_mov);
EC4_SurfCenter(2,:) = Rotation(EC4_SurfCenter(2,:), EC4_center, 2, -EC4_mov);
EC4_SurfCenter(3,:) = Rotation(EC4_SurfCenter(3,:), EC4_center, 2, -EC4_mov);
EC4_SurfCenter(4,:) = Rotation(EC4_SurfCenter(4,:), EC4_center, 2, -EC4_mov);

figure(iter)
BRaum = patch(bauraum,'FaceColor',       [0 0 1.0], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15); 
BRaum.FaceAlpha = 0.4 ;          % Interpolate to find face transparency
hold on
view([45 90 135]);
xlabel('X') 
ylabel('Y')
zlabel('Z')
camlight('headlight');
material('default');

patch(gearbox,'FaceColor',       [0.8 0.8 1], ...
     'FaceLighting',    'gouraud',     ...
     'AmbientStrength', 0.15);

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
Clearance = 5;

cable_radius = 2.3;
radius_of_curvature = 13.8; 
EM_radius = 117;

if iter > 10
%         figure(i+1)
%         Visualize(bauraum,gearbox,EC1,EC2,EC3,EC4);
%         hold on
    if EC2_placement == 1   %Counterclockwise        

        vector1 = EC1_SurfCenter(2,:) - EC1_SurfCenter(1,:);
        vector2 = EC2_SurfCenter(1,:) - EC2_SurfCenter(2,:);
        pt2 = EC1_SurfCenter(1,:) - (Clearance.*vector1)/norm(vector1);
        pt3 = EC2_SurfCenter(2,:) - (Clearance.*vector2)/norm(vector2);
        [splinelength,max_curvature, curvature, collisionwith1, collisionwith2, Penalty_spline, fn] = spline_connection (EC1_SurfCenter(1,:), pt2, pt3, EC2_SurfCenter(2,:),EC1_center, EC2_center, EC1, EC2, bauraum, cable_radius, radius_of_curvature);
        fnplt(fn,'r',2);
        hold on
        vector1 = EC2_SurfCenter(2,:) - EC2_SurfCenter(1,:);
        vector2 = EC3_SurfCenter(1,:) - EC3_SurfCenter(2,:);
        pt2 = EC2_SurfCenter(1,:) - (Clearance.*vector1)/norm(vector1);
        pt3 = EC3_SurfCenter(2,:) - (Clearance.*vector2)/norm(vector2);
        [splinelength,max_curvature, curvature, collisionwith1, collisionwith2, Penalty_spline, fn] = spline_connection (EC2_SurfCenter(1,:), pt2, pt3, EC3_SurfCenter(2,:),EC2_center, EC3_center, EC2, EC3, bauraum, cable_radius, radius_of_curvature);
        fnplt(fn,'r',2);
        hold on
        vector1 = EC1_SurfCenter(4,:) - EC1_SurfCenter(3,:);
        vector2 = EC4_SurfCenter(3,:) - EC4_SurfCenter(4,:);
        pt2 = EC1_SurfCenter(3,:) - (Clearance.*vector1)/norm(vector1);
        pt3 = EC4_SurfCenter(4,:) - (Clearance.*vector2)/norm(vector2);
        [splinelength,max_curvature, curvature, collisionwith1, collisionwith2, Penalty_spline, fn] = spline_connection (EC1_SurfCenter(3,:), pt2, pt3, EC4_SurfCenter(4,:),EC1_center, EC4_center, EC1, EC4, bauraum, cable_radius, radius_of_curvature);
        fnplt(fn,'r',2);
        hold on
        
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
        fnplt(fn_mid,'r',2);
        fnplt(fn_posy,'r',2);
        fnplt(fn_negy,'r',2);
        
    else                    %Clockwise

        vector1 = EC1_SurfCenter(1,:) - EC1_SurfCenter(2,:);
        vector2 = EC2_SurfCenter(2,:) - EC2_SurfCenter(1,:);
        pt2 = EC1_SurfCenter(2,:) - (Clearance.*vector1)/norm(vector1);
        pt3 = EC2_SurfCenter(1,:) - (Clearance.*vector2)/norm(vector2);        
        [splinelength,max_curvature, curvature, collisionwith1, collisionwith2, Penalty_spline, fn] = spline_connection (EC1_SurfCenter(2,:), pt2, pt3, EC2_SurfCenter(1,:),EC1_center, EC2_center, EC1, EC2, bauraum, cable_radius, radius_of_curvature);
        fnplt(fn,'r',2);
        hold on
        vector1 = EC2_SurfCenter(1,:) - EC2_SurfCenter(2,:);
        vector2 = EC3_SurfCenter(2,:) - EC3_SurfCenter(1,:);
        pt2 = EC2_SurfCenter(2,:) - (Clearance.*vector1)/norm(vector1);
        pt3 = EC3_SurfCenter(1,:) - (Clearance.*vector2)/norm(vector2);        
        [splinelength,max_curvature, curvature, collisionwith1, collisionwith2, Penalty_spline, fn] = spline_connection (EC2_SurfCenter(2,:), pt2, pt3, EC3_SurfCenter(1,:),EC2_center, EC3_center, EC2, EC3, bauraum, cable_radius, radius_of_curvature);
        fnplt(fn,'r',2);
        hold on
        vector1 = EC1_SurfCenter(4,:) - EC1_SurfCenter(3,:);
        vector2 = EC4_SurfCenter(3,:) - EC4_SurfCenter(4,:);
        pt2 = EC1_SurfCenter(3,:) - (Clearance.*vector1)/norm(vector1);
        pt3 = EC4_SurfCenter(4,:) - (Clearance.*vector2)/norm(vector2);        
        [splinelength,max_curvature, curvature, collisionwith1, collisionwith2, Penalty_spline, fn] = spline_connection (EC1_SurfCenter(3,:), pt2, pt3, EC4_SurfCenter(4,:),EC1_center, EC4_center, EC1, EC4, bauraum, cable_radius, radius_of_curvature);
        fnplt(fn,'r',2);
        hold on
        
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
        fnplt(fn_mid,'r',2);
        fnplt(fn_posy,'r',2);
        fnplt(fn_negy,'r',2);

    end
end

end