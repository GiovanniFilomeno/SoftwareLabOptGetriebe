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
Cube=stlread("cube.STL");
V = Cube.Points;
V2 = 0.35.*V;
F = Cube.ConnectivityList;
N = faceNormal(Cube);
EC1 = struct('faces',F,'vertices',V2);

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
EM_start = [-29.96 -2295 1];
EM_vector = [-271.7 0 0];
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
fv.vertices = rotation(fv.vertices,Anchor, 2, pi/2);
EM_start = rotation(EM_start, Anchor, 2, pi/2);
EM_vector = rotation(EM_vector, zero, 2, pi/2);
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

%Defining the limits and center of the Electric component
EC1_xmax = max(EC1.vertices(:,1));
EC1_xmin = min(EC1.vertices(:,1));
EC1_ymax = max(EC1.vertices(:,2));
EC1_ymin = min(EC1.vertices(:,2));
EC1_zmax = max(EC1.vertices(:,3));
EC1_zmin = min(EC1.vertices(:,3));
EC1_center = [(EC1_xmax + EC1_xmin)/2 (EC1_ymax + EC1_ymin)/2 (EC1_zmax + EC1_zmin)/2];     

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
iter = 5;
%Create first sample population and evaluate fitness
x = xmin + (xmax-xmin).*rand(pop,1);
y = ymin + (ymax-ymin).*rand(pop,1);
z = zmin + (zmax-zmin).*rand(pop,1);
gb_rot = 0 + (pi/2).*rand(pop,1);
coords = [x,y,z];
d = zeros(pop,1);
for i=1:pop
    
    EC1.vertices = EC1.vertices - EC1_center + coords(i,:);
    EC1_center = coords(i,:);
    
    IN = inpolyhedron(bauraum,EC1.vertices);
    in_bauraum = sum(IN)/length(IN);
    IN = inpolyhedron(fv,EC1.vertices);
    out_gearbox = sum(IN)/length(IN);
    d(i,1) = optimization_function(EC1_center, EM_start, EM_end, in_bauraum, out_gearbox);
    if d(i,1) < 100
        patch(EC1,'FaceColor',       [0.8 0.8 1.0], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15);
     i
     d(i,1)
     pause(0.1)
    end

end
%  [srt,I]=sort(d);
% 
% EC1.vertices = EC1.vertices - EC1_center + coords(I(1,1),:);
% EC1_center = coords(I(1,1),:);
% 
% patch(EC1,'FaceColor',       [0.8 0.8 1.0], ...
%          'EdgeColor',       'none',        ...
%          'FaceLighting',    'gouraud',     ...
%          'AmbientStrength', 0.15);
% srt(1,1)
% pause(1)

%  for g = 2:iter
%       for k = 1: 3: pop
%          %Selection
%         sel_coeff = mean(d);
%         parents = zeros(3,3);
%         for j=1:3
%             rand_sel = 1 + round(pop.*rand());
%             for i = rand_sel:pop
%                 if d(i,1)<sel_coeff
%                     parents(j,1) = coords(i,1);
%                     parents(j,2) = coords(i,2);
%                     parents(j,3) = coords(i,3);
%                     break
%                 end
%             end
%         end
%  
%         %Crossover
%         child = [parents(1,1) parents(2,2) parents(3,3);parents(1,2) parents(2,3) parents(3,1);parents(1,3) parents(2,1) parents(3,2)];
% %         child = [parents(1,1) parents(1,2) parents(1,3);parents(2,1) parents(2,2) parents(2,3);parents(3,1) parents(3,2) parents(3,3)];
% 
% 
% 
% 
%         %Mutation
% %         for i = 1:3
% %              for j = 1:3
% %                 mut_prob = rand();
% %                 if mut_prob<=0.1
% %                     d_child=point_to_line(child(i,1),EM_start,EM_end);
% %                     child(i,j) = child(i,j) + mut_prob*d_child; 
% %                 end
% %              end
% %         end
%         
%         %Write new generation 
%         coords(k,1) = child(1,1);
%         coords(k,2) = child(1,2);
%         coords(k,3) = child(1,3);
%         coords(k+1,1) = child(2,1);
%         coords(k+1,2) = child(2,2);
%         coords(k+1,3) = child(2,3);
%         coords(k+2,1) = child(2,1);
%         coords(k+2,2) = child(2,2);
%         coords(k+2,3) = child(2,3);
%       end
% 
%       for i=1:pop
% 
%             EC1.vertices = EC1.vertices - EC1_center + coords(i,:);
%             EC1_center = coords(i,:);
%             IN = inpolyhedron(bauraum,EC1.vertices);
%             in_bauraum = sum(IN)/length(IN);
%             IN = inpolyhedron(fv,EC1.vertices);
%             out_gearbox = sum(IN)/length(IN);
%             d(i,1) = optimization_function(EC1_center, EM_start, EM_end, in_bauraum, out_gearbox);
% 
%       end
%         [srt,I]=sort(d);
% 
%         EC1.vertices = EC1.vertices - EC1_center + coords(I(1,1),:);
%         EC1_center = coords(I(1,1),:);
% 
%         patch(EC1,'FaceColor',       [0.8 0.8 1.0], ...
%                  'EdgeColor',       'none',        ...
%                  'FaceLighting',    'gouraud',     ...
%                  'AmbientStrength', 0.15);
%        srt(1,1)
%        pause(1)
%  end
  
 
 
function fitness = optimization_function(EC_center, E_start, E_end, in_bauraum, out_gearbox)
    Pen = 2000;

    if EC_center(1,2) > E_start(1,2) || EC_center(1,2) < E_end(1,2)
         distance = point_to_line(EC_center, E_start,E_end);
         fitness = distance + (1-in_bauraum+out_gearbox)*Pen;
     else
         distance = point_to_line(EC_center, E_start,E_end);
         fitness = distance - 135.85 + (1-in_bauraum+out_gearbox)*Pen;
    end

     


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