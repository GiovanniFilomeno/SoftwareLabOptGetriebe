clc
clear all


FV=stlread("Gearbox1.STL");
V = FV.Points;
F = FV.ConnectivityList;
N = faceNormal(FV);

fv = struct('faces',F,'vertices',V);

BR=stlread("Bauraum.STL");
V = BR.Points;
V2 = 0.004.*V;
F = BR.ConnectivityList;
N = faceNormal(BR);

bauraum = struct('faces',F,'vertices',V2);

Cube=stlread("cube.STL");
V = Cube.Points;
V2 = .0005.*V;
F = Cube.ConnectivityList;
N = faceNormal(Cube);

cube = struct('faces',F,'vertices',V2);


BRaum = patch(bauraum,'FaceColor',       [0.8 0 1.0], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15);
BRaum.FaceVertexAlphaData = 0.2;    % Set constant transparency 
BRaum.FaceAlpha = 'flat' ;          % Interpolate to find face transparency
pause(1)
shaftx = [fv.vertices(637,1) fv.vertices(638,1)];
shafty = [fv.vertices(637,2) fv.vertices(638,2)];
shaftz = [fv.vertices(637,3) fv.vertices(638,3)];
hold on
view([45 90 135]);
pause(1)
xlabel('X') 
ylabel('Y')
zlabel('Z')
%axis([0 10 0 10 0 10])
camlight('headlight');
material('default');

Anchor = fv.vertices(637,:);

point = bauraum.vertices(36,:);
vector = [bauraum.vertices(43,1)-bauraum.vertices(38,1) bauraum.vertices(43,2)-bauraum.vertices(38,2) bauraum.vertices(43,3)-bauraum.vertices(38,3)];
anglex = acos(vector(1,1)/norm(vector));
angley = acos(vector(1,2)/norm(vector));
anglez = acos(vector(1,3)/norm(vector));

fv.vertices = fv.vertices - Anchor + point;
Anchor = fv.vertices(637,:);
shaftx = [fv.vertices(637,1) fv.vertices(637,1) + vector(1,1)];
shafty = [fv.vertices(637,2) fv.vertices(637,2) + vector(1,2)];
shaftz = [fv.vertices(637,3) fv.vertices(637,3) + vector(1,3)];
plot3(shaftx', shafty', shaftz','LineWidth',2,'Color','red')
pause(2)
fv.vertices = rotation(fv.vertices,Anchor, 1, anglex);
fv.vertices = rotation(fv.vertices,Anchor, 2, angley);
fv.vertices = rotation(fv.vertices,Anchor, 3, anglez);
fv.vertices = fv.vertices - [0,0,0.1];
patch(fv,'FaceColor',       [0.8 0.8 1.0], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15);
pause(1)
shaftx = [fv.vertices(637,1) fv.vertices(638,1)];
shafty = [fv.vertices(637,2) fv.vertices(638,2)];
shaftz = [fv.vertices(637,3) fv.vertices(638,3)];
hold on
plot3(shaftx', shafty', shaftz','LineWidth',2,'Color','red')
pause(1)

xmax = max(bauraum.vertices(:,1));
xmin = min(bauraum.vertices(:,1));
ymax = max(bauraum.vertices(:,2));
ymin = min(bauraum.vertices(:,2));
zmax = max(bauraum.vertices(:,3));
zmin = min(bauraum.vertices(:,3));

pop = 100;
iter = 5;
%Create first sample population and evaluate fitness
coords = xmin + (xmax-xmin).*rand(pop,1);
y = ymin + (ymax-ymin).*rand(pop,1);
z = zmin + (zmax-zmin).*rand(pop,1);
angles = (2*pi).*rand(pop,3);
coords = [coords,y,z];
d = zeros(pop,1);
in = zeros(pop,1);
colision = zeros(pop,1);
for i=1:pop
    d(i,1) = norm(Anchor-coords(i,:));
    cube.vertices = cube.vertices - (cube.vertices(1,:) - coords(i,:));
    cube.vertices = rotation(cube.vertices,cube.vertices(1,:), 1, angles(i,1));
    cube.vertices = rotation(cube.vertices,cube.vertices(1,:), 2, angles(i,2));
    cube.vertices = rotation(cube.vertices,cube.vertices(1,:), 3, angles(i,3));
    IN = inpolyhedron(bauraum,cube.vertices);
    in(i,1) = sum(IN)/length(IN);
    IN = inpolyhedron(fv,cube.vertices);
    colision(i,1) = sum(IN)/length(IN);
    cube.vertices = rotation(cube.vertices,cube.vertices(1,:), 1, -angles(i,1));
    cube.vertices = rotation(cube.vertices,cube.vertices(1,:), 2, -angles(i,2));
    cube.vertices = rotation(cube.vertices,cube.vertices(1,:), 3, -angles(i,3));
end
[srt,I]=sort(d);
for i=1:pop
    if in(I(i),1) == 1 && colision(I(i),1) == 0
    best = coords(I(i),:);
    best_angle = angles(I(i),:);
    break
    end
end
cube.vertices = cube.vertices - (cube.vertices(1,:) - best(1,:));
cube.vertices = rotation(cube.vertices,cube.vertices(1,:), 1, best_angle(1,1));
cube.vertices = rotation(cube.vertices,cube.vertices(1,:), 2, best_angle(1,2));
cube.vertices = rotation(cube.vertices,cube.vertices(1,:), 3, best_angle(1,3));

patch(cube,'FaceColor',       [0.8 0.8 1.0], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15);

cube.vertices = rotation(cube.vertices,cube.vertices(1,:), 1, -best_angle(1,1));
cube.vertices = rotation(cube.vertices,cube.vertices(1,:), 2, -best_angle(1,2));
cube.vertices = rotation(cube.vertices,cube.vertices(1,:), 3, -best_angle(1,3));
pause(1)

for g = 2:iter
     for k = 1: 3: pop
        %Selection
        sel_coeff = mean(d);
        parents = zeros(3,3);
        for j=1:3
            rand_sel = 1 + round(pop.*rand());
            for i = rand_sel:pop
                if d(i,1)<sel_coeff && in(i,1) == 1 && colision(i,1) == 0
                    parent(j,1) = coords(i,1);
                    parent(j,2) = coords(i,2);
                    parent(j,3) = coords(i,3);
                    break
                end
            end
        end

        %Crossover
        child = [parent(1,1) parent(2,2) parent(3,3);parent(1,2) parent(2,3) parent(3,1);parent(1,3) parent(2,1) parent(3,2)];



        %Mutation
        for i = 1:3
             for j = 1:3
                mut_prob = rand();
                if mut_prob<=0.1
                    child(i,j) = child(i,j) - (child(i,j)-Anchor(1,j))*(sel_coeff*0.1);                   
                end
             end
        end
%         
        %Write new generation 
        coords(k,1) = child(1,1);
        coords(k,2) = child(1,2);
        coords(k,3) = child(1,3);
        coords(k+1,1) = child(2,1);
        coords(k+1,2) = child(2,2);
        coords(k+1,3) = child(2,3);
        coords(k+2,1) = child(2,1);
        coords(k+2,2) = child(2,2);
        coords(k+2,3) = child(2,3);
     end

    for i=1:pop
        d(i,1) = norm(Anchor-coords(i,:));
        cube.vertices = cube.vertices - (cube.vertices(1,:) - coords(i,:));
        cube.vertices = rotation(cube.vertices,cube.vertices(1,:), 1, angles(i,1));
        cube.vertices = rotation(cube.vertices,cube.vertices(1,:), 2, angles(i,2));
        cube.vertices = rotation(cube.vertices,cube.vertices(1,:), 3, angles(i,3));
        IN = inpolyhedron(bauraum,cube.vertices);
        in(i,1) = sum(IN)/length(IN);
        IN = inpolyhedron(fv,cube.vertices);
        colision(i,1) = sum(IN)/length(IN);
        cube.vertices = rotation(cube.vertices,cube.vertices(1,:), 1, -angles(i,1));
        cube.vertices = rotation(cube.vertices,cube.vertices(1,:), 2, -angles(i,2));
        cube.vertices = rotation(cube.vertices,cube.vertices(1,:), 3, -angles(i,3));
    end
    [srt,I]=sort(d);
    for i=1:pop
        if in(I(i),1) == 1 && colision(I(i),1) == 0
        best = coords(I(i),:);
        best_angle = angles(I(i),:);
        break
        end
    end
    cube.vertices = cube.vertices - (cube.vertices(1,:) - best(1,:));
    cube.vertices = rotation(cube.vertices,cube.vertices(1,:), 1, best_angle(1,1));
    cube.vertices = rotation(cube.vertices,cube.vertices(1,:), 2, best_angle(1,2));
    cube.vertices = rotation(cube.vertices,cube.vertices(1,:), 3, best_angle(1,3));

    patch(cube,'FaceColor',       [0.8 0.8 1.0], ...
             'EdgeColor',       'none',        ...
             'FaceLighting',    'gouraud',     ...
             'AmbientStrength', 0.15);

    cube.vertices = rotation(cube.vertices,cube.vertices(1,:), 1, -best_angle(1,1));
    cube.vertices = rotation(cube.vertices,cube.vertices(1,:), 2, -best_angle(1,2));
    cube.vertices = rotation(cube.vertices,cube.vertices(1,:), 3, -best_angle(1,3));
    pause(1)
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