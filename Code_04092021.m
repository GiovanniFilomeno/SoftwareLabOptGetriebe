clear all
clc

%Read Gearbox
FV=stlread("Gearbox1.STL");
V = FV.Points;
V1 = 1000.*V;
F = FV.ConnectivityList;
N = faceNormal(FV);
fv = struct('faces',F,'vertices',V1);

%Read Assembly space
BR=stlread("Bauraum.STL");
V = BR.Points;
F = BR.ConnectivityList;
N = faceNormal(BR);
bauraum = struct('faces',F,'vertices',V);

%Read electric component
% Cube=stlread("cube.STL");
% V = Cube.Points;
% V2 = 0.35.*V;
% F = Cube.ConnectivityList;
% N = faceNormal(Cube);
% EC1 = struct('faces',F,'vertices',V2);

Cube=stlread("cube.STL");
V = Cube.Points;
V2 = zeros(26,3);
V2(:,1) = 1.3477.*V(:,1);
V2(:,2) = 0.58428.*V(:,2);
V2(:,3) = 0.7778513.*V(:,3);
V2 = 0.25.*V2;
V3 = 0.25.*V;
F = Cube.ConnectivityList;
N = faceNormal(Cube);
EC1 = struct('faces',F,'vertices',V2);

EC2 = struct('faces',F,'vertices',V3);



%Plot Bauraum
BRaum = patch(bauraum,'FaceColor',       [0 0 1.0], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15);
BRaum.FaceVertexAlphaData = 0.2;    % Set constant transparency 
BRaum.FaceAlpha = 'flat' ;          % Interpolate to find face transparency
pause(1)
hold on
view([45 90 135]);
pause(1)
xlabel('X') 
ylabel('Y')
zlabel('Z')
%axis([0 10 0 10 0 10])
camlight('headlight');
material('default');


%Definition(manually) of the axis of the electric motor
EM_start = [-91.61 -2295 1];
EM_vector = [-210 0 0];
EM_end = EM_start + EM_vector;
zero = [0 0 0];

%Definition(manually of the axis of the gearbox shaft 
Anchor = [14 -2518 -34.015];
point = [2539 14 -40];
vector = [0 -1 0];
anglex = acos(vector(1,1)/norm(vector));
angley = acos(vector(1,2)/norm(vector));
anglez = acos(vector(1,3)/norm(vector));

%Positioning of the gearbox inside the Bauraum
fv.vertices = fv.vertices - Anchor + point;
EM_start = EM_start - Anchor + point;
Anchor = point;
shaftx = [2539 2539];
shafty = [14 -130];
shaftz = [-40 -40];
% shaftx = [fv.vertices(637,1) fv.vertices(637,1) + vector(1,1)];
% shafty = [fv.vertices(637,2) fv.vertices(637,2) + vector(1,2)];
% shaftz = [fv.vertices(637,3) fv.vertices(637,3) + vector(1,3)];
plot3(shaftx', shafty', shaftz','LineWidth',2,'Color','red')
pause(2)
fv.vertices = rotation(fv.vertices,Anchor, 1, anglex);
fv.vertices = rotation(fv.vertices,Anchor, 2, angley);
fv.vertices = rotation(fv.vertices,Anchor, 3, anglez);
EM_start = rotation(EM_start, Anchor, 1, anglex);
EM_start = rotation(EM_start, Anchor, 2, angley);
EM_start = rotation(EM_start, Anchor, 3, anglez);
EM_vector = rotation(EM_vector, zero, 1, anglex);
EM_vector = rotation(EM_vector, zero, 2, angley);
EM_vector = rotation(EM_vector, zero, 3, anglez);
EM_end = EM_start + EM_vector;
EM_shaftx = [EM_start(1,1) EM_start(1,1) + EM_vector(1,1)];
EM_shafty = [EM_start(1,2) EM_start(1,2) + EM_vector(1,2)];
EM_shaftz = [EM_start(1,3) EM_start(1,3) + EM_vector(1,3)];
hold on


%Defining the limits of the bauraum to generate the population inside.
xmax = max(bauraum.vertices(:,1));
xmin = min(bauraum.vertices(:,1));
ymax = max(bauraum.vertices(:,2));
ymin = min(bauraum.vertices(:,2));
zmax = max(bauraum.vertices(:,3));
zmin = min(bauraum.vertices(:,3));

%Rotation of the Gearbox 90 deg to be completely inside the bauraum
fv.vertices = rotation(fv.vertices,Anchor, 2, 13*pi/36);
EM_start = rotation(EM_start, Anchor, 2, 13*pi/36);
EM_vector = rotation(EM_vector, zero, 2, 13*pi/36);
EM_end = EM_start + EM_vector;
EM_shaftx = [EM_start(1,1) EM_start(1,1) + EM_vector(1,1)];
EM_shafty = [EM_start(1,2) EM_start(1,2) + EM_vector(1,2)];
EM_shaftz = [EM_start(1,3) EM_start(1,3) + EM_vector(1,3)];
hold on
plot3(EM_shaftx', EM_shafty', EM_shaftz','LineWidth',2,'Color','green')
fv.FaceVertexAlphaData = 0.2;    % Set constant transparency 
fv.FaceAlpha = 'flat' ;          % Interpolate to find face transparency
patch(fv,'FaceColor',       [0.8 0.8 1.0], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15);
pause(1)
IN = inpolyhedron(bauraum,fv.vertices);
out_gearbox = sum(IN)/length(IN)
%Defining the limits and center of the Electric component
EC1_xmax = max(EC1.vertices(:,1));
EC1_xmin = min(EC1.vertices(:,1));
EC1_ymax = max(EC1.vertices(:,2));
EC1_ymin = min(EC1.vertices(:,2));
EC1_zmax = max(EC1.vertices(:,3));
EC1_zmin = min(EC1.vertices(:,3));
EC1_center = [(EC1_xmax + EC1_xmin)/2 (EC1_ymax + EC1_ymin)/2 (EC1_zmax + EC1_zmin)/2];     

EC2_xmax = max(EC2.vertices(:,1));
EC2_xmin = min(EC2.vertices(:,1));
EC2_ymax = max(EC2.vertices(:,2));
EC2_ymin = min(EC2.vertices(:,2));
EC2_zmax = max(EC2.vertices(:,3));
EC2_zmin = min(EC2.vertices(:,3));
EC2_center = [(EC2_xmax + EC2_xmin)/2 (EC2_ymax + EC2_ymin)/2 (EC2_zmax + EC2_zmin)/2]; 

% EC1.vertices = rotation(EC1.vertices, EC1_center, 2, pi/4);

% point_EC1 = [2239 -190.5 295.2];
% EC1.vertices = EC1.vertices - EC1_center + point_EC1;
% EC1_center = point_EC1;
%      
% IN = inpolyhedron(bauraum,EC1.vertices);     
% br_ec1 = sum(IN)/length(IN);
% IN = inpolyhedron(fv,EC1.vertices);     
% fv_ec1 = sum(IN)/length(IN);
% patch(EC1,'FaceColor',       [0.8 0.8 1.0], ...
%          'EdgeColor',       'none',        ...
%          'FaceLighting',    'gouraud',     ...
%          'AmbientStrength', 0.15);
% pause(1)
     
     
     
%Genetic Algorithm implementation
pop = 100;
iter = 20;
Penalty = 2000;
%Create first sample population and evaluate fitness
x1 = xmin + (xmax-xmin).*rand(pop,1);
y1= ymin + (ymax-ymin).*rand(pop,1);
z1 = zmin + (zmax-zmin).*rand(pop,1);
x2 = xmin + (xmax-xmin).*rand(pop,1);
y2= ymin + (ymax-ymin).*rand(pop,1);
z2 = zmin + (zmax-zmin).*rand(pop,1);
gb_rot = 0 + (pi/2).*rand(pop,1);
coords1 = [x1,y1,z1];
coords2 = [x2,y2,z2];
d1 = zeros(pop,1);
d2 = zeros(pop,1);
d = zeros(pop,1);
for i=1:pop
    
    EC1.vertices = EC1.vertices - EC1_center + coords1(i,:);
    EC1_center = coords1(i,:);
    EC2.vertices = EC2.vertices - EC2_center + coords2(i,:);
    EC2_center = coords2(i,:);
    
    IN = inpolyhedron(bauraum,EC1.vertices);
    in1_bauraum = sum(IN)/length(IN);
    IN = inpolyhedron(fv,EC1.vertices);
    out1_gearbox = sum(IN)/length(IN);
    IN = inpolyhedron(bauraum,EC2.vertices);
    in2_bauraum = sum(IN)/length(IN);
    IN = inpolyhedron(fv,EC2.vertices);
    out2_gearbox = sum(IN)/length(IN);
    IN = inpolyhedron(EC1,EC2.vertices);
    collision12 = sum(IN)/length(IN);
%     d1(i,1) = optimization_function(EC1_center, EM_start, EM_end, in1_bauraum, out1_gearbox);
%     d2(i,1) = optimization_function(EC2_center, EM_start, EM_end, in2_bauraum, out2_gearbox);
    d(i,1) = optimization_function(EC1_center, EC2_center, EM_start, EM_end, in1_bauraum, out1_gearbox, in2_bauraum, out2_gearbox, collision12);

end
%  [srt1,I1]=sort(d1);
%  [srt2,I2]=sort(d2);
%  
% for i = 1:pop
%     EC1.vertices = EC1.vertices - EC1_center + coords1(i,:);
%     EC1_center = coords1(i,:);
%     EC2.vertices = EC2.vertices - EC2_center + coords2(i,:);
%     EC2_center = coords2(i,:);
%     EC1.vertices = EC1.vertices - EC1_center + coords1(I1(1,1),:);
%     EC1_center = coords1(I1(i,1),:);
%     EC2.vertices = EC2.vertices - EC2_center + coords2(I2(1,1),:);
%     EC2_center = coords2(I2(i,1),:);
%     
%     IN = inpolyhedron(EC1,EC2.vertices);
%     collision12 = sum(IN)/length(IN);
%     d(i,1) = srt1(i,1) + srt2(i,1) + Penalty*collision12;
%     
% end
    
[srt,I]=sort(d);

EC1.vertices = EC1.vertices - EC1_center + coords1(I(1,1),:);
EC1_center = coords1(I(1,1),:);
EC2.vertices = EC2.vertices - EC2_center + coords2(I(1,1),:);
EC2_center = coords2(I(1,1),:);

patch(EC1,'FaceColor',       [0.8 0.8 1.0], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15);
patch(EC2,'FaceColor',       [0 0.8 1.0], ...
     'EdgeColor',       'none',        ...
     'FaceLighting',    'gouraud',     ...
     'AmbientStrength', 0.15);
fprintf('Fitness of Generation 1: %d \n',srt(1,1));
pause(1)

best1 = zeros(pop/2,3);
best2 = zeros(pop/2,3);
for i=1:pop/2
    best1(i,:) = coords1(I(i,1),:);
    best2(i,:) = coords2(I(i,1),:);
end

 for g = 2:iter
      for k = 1: pop/2
         %Selection
        sel_coeff = mean(d);
        parents1 = zeros(2,3);
        d_parents1 = zeros(1,2);
        parents2 = zeros(2,3);
        d_parents2 = zeros(1,2);
        for j=1:2
            rand_sel = 1 + round(pop.*rand());
            for i = rand_sel:pop
                if d(i,1)<sel_coeff
                    parents1(j,1) = coords1(i,1);
                    parents1(j,2) = coords1(i,2);
                    parents1(j,3) = coords1(i,3);
                    parents2(j,1) = coords2(i,1);
                    parents2(j,2) = coords2(i,2);
                    parents2(j,3) = coords2(i,3);
                    d_parents1(1,j) = d(i,1);
                    d_parents2(1,j) = d(i,1);
                    
                    break
                end
            end
        end
 
        %Crossover
        child1 = zeros(2,3);
        child2 = zeros(2,3);
        alpha = 0.3;
        beta = 0.7;
        child1(1,:) = alpha.*parents1(1,:) + (1-alpha).*parents1(2,:);
        child1(2,:) = beta.*parents1(1,:) + (1-beta).*parents1(2,:);
        child2(1,:) = alpha.*parents2(1,:) + (1-alpha).*parents2(2,:);
        child2(2,:) = beta.*parents2(1,:) + (1-beta).*parents2(2,:);
%        child = [parents(1,1) parents(2,2) parents(3,3);parents(1,2) parents(2,3) parents(3,1);parents(1,3) parents(2,1) parents(3,2)];
%        child = [parents(1,1) parents(1,2) parents(1,3);parents(2,1) parents(2,2) parents(2,3);parents(3,1) parents(3,2) parents(3,3)];




%       Mutation
        for i = 1:2
             child1_start = EM_start - child1(i,:);
             child1_end = EM_end - child1(i,:);
             child2_start = EM_start - child2(i,:);
             child2_end = EM_end - child2(i,:);
             for j = 1:3
                mut_prob = rand();
                if mut_prob<=0.75
                    if child1(i,2) > EM_start(1,2)
                        child1(i,j) = child1(i,j) + 0.05*child1_start(1,j);
                    elseif child1(i,2) < EM_end(1,2)
                        child1(i,j) = child1(i,j) + 0.05*child1_end(1,j);
                    else
                        child1(i,j) = child1(i,j) + 0.05*child1_end(1,j);
                        child1(i,2) =child1(i,2) - 0.05*child1_end(1,2); 
                    end
                    
                    if child2(i,2) > EM_start(1,2)
                        child2(i,j) = child2(i,j) + 0.05*child2_start(1,j);
                    elseif child2(i,2) < EM_end(1,2)
                        child2(i,j) = child2(i,j) + 0.05*child2_end(1,j);
                    else
                        child2(i,j) = child2(i,j) + 0.05*child2_end(1,j);
                        child2(i,2) =child2(i,2) - 0.05*child2_end(1,2); 
                    end
                end
             end
        end
        
        %Write new generation 
        coords1(k,1) = child1(1,1);
        coords1(k,2) = child1(1,2);
        coords1(k,3) = child1(1,3);
        coords1(k+1,1) = child1(2,1);
        coords1(k+1,2) = child1(2,2);
        coords1(k+1,3) = child1(2,3);
        coords2(k,1) = child2(1,1);
        coords2(k,2) = child2(1,2);
        coords2(k,3) = child2(1,3);
        coords2(k+1,1) = child2(2,1);
        coords2(k+1,2) = child2(2,2);
        coords2(k+1,3) = child2(2,3);
      end
      
      for i=(pop/2+1):pop
          coords1(i,:) = best1(i-(pop/2),:);
          coords2(i,:) = best2(i-(pop/2),:);
      end
          
      
%         for i=1:pop
% 
%             EC1.vertices = EC1.vertices - EC1_center + coords1(i,:);
%             EC1_center = coords1(i,:);
%             EC2.vertices = EC2.vertices - EC2_center + coords2(i,:);
%             EC2_center = coords2(i,:);
% 
%             IN = inpolyhedron(bauraum,EC1.vertices);
%             in1_bauraum = sum(IN)/length(IN);
%             IN = inpolyhedron(fv,EC1.vertices);
%             out1_gearbox = sum(IN)/length(IN);
%             IN = inpolyhedron(bauraum,EC2.vertices);
%             in2_bauraum = sum(IN)/length(IN);
%             IN = inpolyhedron(fv,EC2.vertices);
%             out2_gearbox = sum(IN)/length(IN);
%             IN = inpolyhedron(EC1,EC2.vertices);
%             collision12 = sum(IN)/length(IN);
% %           d1(i,1) = optimization_function(EC1_center, EM_start, EM_end, in1_bauraum, out1_gearbox);
% %           d2(i,1) = optimization_function(EC2_center, EM_start, EM_end, in2_bauraum, out2_gearbox);
%             d(i,1) = optimization_function(EC1_center, EC2_center, EM_start, EM_end, in1_bauraum, out1_gearbox, in2_bauraum, out2_gearbox, collision12);


%         end
%          [srt1,I1]=sort(d1);
%          [srt2,I2]=sort(d2);
% 
%         for i = 1:pop
%             EC1.vertices = EC1.vertices - EC1_center + coords1(i,:);
%             EC1_center = coords1(i,:);
%             EC2.vertices = EC2.vertices - EC2_center + coords2(i,:);
%             EC2_center = coords2(i,:);
% 
%             IN = inpolyhedron(EC1,EC2.vertices);
%             collision12 = sum(IN)/length(IN);
%             d(i,1) = d1(i,1) + d2(i,1) + Penalty*collision12;
% 
%         end

%         [srt,I]=sort(d);
%         EC1.vertices = EC1.vertices - EC1_center + coords1(I1(1,1),:);
%         EC1_center = coords1(I1(1,1),:);
%         EC2.vertices = EC2.vertices - EC2_center + coords2(I2(1,1),:);
%         EC2_center = coords2(I2(1,1),:);


      for i=1:pop

            EC1.vertices = EC1.vertices - EC1_center + coords1(i,:);
            EC1_center = coords1(i,:);
            EC2.vertices = EC2.vertices - EC2_center + coords2(i,:);
            EC2_center = coords2(i,:);

            IN = inpolyhedron(bauraum,EC1.vertices);
            in1_bauraum = sum(IN)/length(IN);
            IN = inpolyhedron(fv,EC1.vertices);
            out1_gearbox = sum(IN)/length(IN);
            IN = inpolyhedron(bauraum,EC2.vertices);
            in2_bauraum = sum(IN)/length(IN);
            IN = inpolyhedron(fv,EC2.vertices);
            out2_gearbox = sum(IN)/length(IN);
            IN = inpolyhedron(EC1,EC2.vertices);
            collision12 = sum(IN)/length(IN);
            d(i,1) = optimization_function(EC1_center,EC2_center, EM_start, EM_end, in1_bauraum, out1_gearbox, in2_bauraum, out2_gearbox, collision12);

      end
        [srt,I]=sort(d);
        for i=1:pop/2
            best1(i,:) = coords1(I(i,1),:);
            best2(i,:) = coords2(I(i,1),:);
        end
        

        EC1.vertices = EC1.vertices - EC1_center + coords1(I(1,1),:);
        EC1_center = coords1(I(1,1),:);
        EC2.vertices = EC2.vertices - EC2_center + coords2(I(1,1),:);
        EC2_center = coords2(I(1,1),:);
        if g ==20
            
        patch(EC1,'FaceColor',       [0.8 0.8 1.0], ...
                 'EdgeColor',       'green',        ...
                 'FaceLighting',    'gouraud',     ...
                 'AmbientStrength', 0.15);
        patch(EC2,'FaceColor',       [0 0.8 1.0], ...
             'EdgeColor',       'red',        ...
             'FaceLighting',    'gouraud',     ...
             'AmbientStrength', 0.15);
        end 
       fprintf('Selection coefficient of Generation %d: %d \n',g-1,sel_coeff);      
       fprintf('Fitness of Generation %d: %d \n' ,g,srt(1,1));
       pause(1)
 end


% function fitness = optimization_function(EC_center, E_start, E_end, in_bauraum, out_gearbox) 
%     Pen = 2000;
% 
%     if EC_center(1,2) > E_start(1,2) || EC_center(1,2) < E_end(1,2)
%          distance = point_to_line(EC_center, E_start,E_end);
%          fitness = distance + (1-in_bauraum+out_gearbox)*Pen;
%      else
%          distance = point_to_line(EC_center, E_start,E_end);
%          fitness = distance - 135.85 + (1-in_bauraum+out_gearbox)*Pen;
%     end
% 
% end

function fitness = optimization_function(EC1_center,EC2_center, E_start, E_end, in1_bauraum, out1_gearbox,in2_bauraum, out2_gearbox, collision12) 
    Pen = 2000;

    if EC1_center(1,2) > E_start(1,2) || EC1_center(1,2) < E_end(1,2)
         distance = point_to_line(EC1_center, E_start,E_end);
         fitness1 = distance + (1-in1_bauraum+out1_gearbox)*Pen;
     else
         distance = point_to_line(EC1_center, E_start,E_end);
         fitness1 = distance - 135.85 + (1-in1_bauraum+out1_gearbox)*Pen;
    end
    
    if EC2_center(1,2) > E_start(1,2) || EC2_center(1,2) < E_end(1,2)
         distance = point_to_line(EC2_center, E_start,E_end);
         fitness2 = distance + (1-in2_bauraum+out2_gearbox)*Pen;
     else
         distance = point_to_line(EC2_center, E_start,E_end);
         fitness2 = distance - 135.85 + (1-in2_bauraum+out2_gearbox)*Pen;
    end    
    
    fitness = fitness1 + fitness2 + Pen*collision12;

 end

function d = point_to_line(pt,v1,v2)
      a = v1 - v2;
      b = pt - v2;

      d = norm(cross(a,b)) / norm(a);
end

function vertex = rotation(V, anchor,indice, angle)

    Rz = [ cos(angle), -sin(angle), 0 ;
          sin(angle), cos(angle), 0 ;
    0, 0, 1 ];
    Ry = [ cos(angle), 0, sin(angle) ;
    0, 1, 0 ;
          -sin(angle), 0, cos(angle) ];
    Rx = [ 1, 0, 0 ;
    0, cos(angle), -sin(angle);
    0, sin(angle), cos(angle) ];

    if(indice==1)
           vertex = (V-anchor)*Rx + anchor;
    end
    if(indice==2)
           vertex = (V-anchor)*Ry + anchor;
    end
    if(indice==3)
           vertex = (V-anchor)*Rz + anchor;
    end
    
end 