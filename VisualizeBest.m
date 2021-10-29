function VisualizeBest(iter,bauraum, gearbox,Anchor, EC1, EC1_center, EC2, EC2_center, EC3, EC3_center, EC4, EC4_center,EM_start,Fitness,r,y,theta,Rotations)

%Move the components and associated variables to the position of that chromosome
%Gearbox
[srt,I] = sort(Fitness);

gearbox.vertices = Rotation(gearbox.vertices, Anchor, 2, Rotations(I(1,1),1));
EM_start = Rotation(EM_start, Anchor, 2, Rotations(I(1,1),1));

coords1(I(1,1),:) = [EM_start(1)+ r(I(1,1),1).*sin(theta(I(1,1),1)) y(I(1,1),1) EM_start(3) + r(I(1,1),1).*cos(theta(I(1,1),1))];
coords2(I(1,1),:) = [EM_start(1)+ r(I(1,1),2).*sin(theta(I(1,1),2)) y(I(1,1),2) EM_start(3) + r(I(1,1),2).*cos(theta(I(1,1),2))];
coords3(I(1,1),:) = [EM_start(1)+ r(I(1,1),3).*sin(theta(I(1,1),3)) y(I(1,1),3) EM_start(3) + r(I(1,1),3).*cos(theta(I(1,1),3))];
coords4(I(1,1),:) = [EM_start(1)+ r(I(1,1),4).*sin(theta(I(1,1),4)) y(I(1,1),4) EM_start(3) + r(I(1,1),4).*cos(theta(I(1,1),4))];

%IGBT

EC1_mov = theta(I(1,1),1) + Rotations(I(1,1),2);
EC1.vertices = EC1.vertices - EC1_center + coords1(I(1,1),:);
EC1_center = coords1(I(1,1),:);
EC1.vertices = Rotation(EC1.vertices, EC1_center, 2, -EC1_mov);


%Capacitor

EC2_mov = theta(I(1,1),2) + Rotations(I(1,1),3);
EC2.vertices = EC2.vertices - EC2_center + coords2(I(1,1),:);
EC2_center = coords2(I(1,1),:);
EC2.vertices = Rotation(EC2.vertices, EC2_center, 2, -EC2_mov);


%EMI-filter

EC3_mov = theta(I(1,1),3) + Rotations(I(1,1),4);
EC3.vertices = EC3.vertices - EC3_center + coords3(I(1,1),:);
EC3_center = coords3(I(1,1),:);
EC3.vertices = Rotation(EC3.vertices, EC3_center, 2, -EC3_mov);


%Control board
EC4_mov = theta(I(1,1),4) + Rotations(I(1,1),5);
EC4.vertices = EC4.vertices - EC4_center + coords4(I(1,1),:);
EC4_center = coords4(I(1,1),:);
EC4.vertices = Rotation(EC4.vertices, EC4_center, 2, -EC4_mov);

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


end