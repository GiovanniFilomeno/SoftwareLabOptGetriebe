clc
clear

%% Load STL mesh
% Stereolithography (STL) files are a common format for storing mesh data. STL
% meshes are simply a collection of triangular faces. This type of model is very
% suitable for use with MATLAB's PATCH graphics object.

% Import an STL mesh, returning a PATCH-compatible face-vertex structure
fv= stlread('Gearbox.stl');

%% Render
% The model is rendered with a PATCH graphics object. We also add some dynamic
% lighting, and adjust the material properties to change the specular
% highlighting.

% patch(fv,'FaceColor',       [0.8 0.8 1.0], ...
%          'EdgeColor',       'none',        ...
%          'FaceLighting',    'gouraud',     ...
%          'AmbientStrength', 0.15);

% Add a camera light, and tone down the specular highlighting
% camlight('headlight');
% material('dull');

% Fix the axes scaling, and set a nice view angle
% axis('image');
% view([45 90]);

% scale=5.0;
% fv.vertices=fv.vertices*scale;
% fv.vertices = rotation(fv.vertices, 2, pi/2);
%% Cylinders

% green zone cylinders
green_Cylcoor = {};
green_radius = {};
green_numberofcyl = 4;

figure

% cylinder 1
%Base at 
x0=0;y0=0;z0=0;
r = 0.3;
h=1;
[x,y,z]=cylinder(r,100);
x=x+x0;
y=y+y0;
z=z*h+z0;
g_vx = [0,0,0];
g_vy = [0,0,5];
green_radius = [green_radius r];
green_Cylcoor = [green_Cylcoor; [x0, y0, z0]];
surf(x,y,z, 'EdgeColor', 'g', 'FaceColor', 'g', 'FaceAlpha', 0.4)
hold on
xlabel('X')
ylabel('Y')
zlabel('Z')
view([45 90]);
drivex = [x0 x0];
drivey = [y0 y0];
drivez = [z0 h];
plot3(drivex',drivey',drivez', 'LineWidth',2)
hold on


% cylinder 2
hold on;
% Base at 
x0=-5;y0=0;z0=0;
r = 3.5;
h=5;
[x,y,z]=cylinder(r,100);
x=x+x0;
y=y+y0;
z=z*h+z0;
g_vx = [g_vx;[-5,0,0]];
g_vy = [g_vy;[-5,0,5]];
green_radius = [green_radius r];
green_Cylcoor = [green_Cylcoor; [x0, y0, z0]];
% surf(x,y,z, 'EdgeColor', 'g', 'FaceColor', 'g', 'FaceAlpha', 0.4)


% cylinder 3
hold on;
%Base at 
x0=-5;y0=-5;z0=0;
r = 3.5;
h=5;
[x,y,z]=cylinder(r,100);
x=x+x0;
y=y+y0;
z=z*h+z0;
g_vx = [g_vx;[-5,-5,0]];
g_vy = [g_vy;[-5,-5,5]];
green_radius = [green_radius r];
green_Cylcoor = [green_Cylcoor; [x0, y0, z0]];
% surf(x,y,z, 'EdgeColor', 'g', 'FaceColor', 'g', 'FaceAlpha', 0.4)


% cylinder 4
hold on;
%Base at 
x0=0;y0=-5;z0=0;
r = 3.5;
h=5;
[x,y,z]=cylinder(r,100);
x=x+x0;
y=y+y0;
z=z*h+z0;
g_vx = [g_vx;[0,-5,0]];
g_vy = [g_vy;[0,-5,5]];
green_radius = [green_radius r];
green_Cylcoor = [green_Cylcoor; [x0, y0, z0]];
% surf(x,y,z, 'EdgeColor', 'g', 'FaceColor', 'g', 'FaceAlpha', 0.4)

% transforming the data
green_Cylcoor = cell2mat(green_Cylcoor);
green_radius = cell2mat(green_radius);

%% red zone cyl
red_Cylcoor = {};
red_radius = {};
red_numberofcyl = 2;

hold on;
%Base at 
x0=0;y0=-3;z0=0;
r = 1.5;
h=5;
[x,y,z]=cylinder(r,100);
x=x+x0;
y=y+y0;
z=z*h+z0;
r_vx = [0,-3,0];
r_vy = [0,-3,5];
red_radius = [red_radius r];
red_Cylcoor =  [red_Cylcoor; [x0, y0, z0]];
% surf(x,y,z, 'EdgeColor', 'r', 'FaceColor', 'r', 'FaceAlpha', 0.4)


hold on;
%Base at 
x0=-3;y0=-3;z0=0;
r = 1;
h=5;
[x,y,z]=cylinder(r,100);
x=x+x0;
y=y+y0;
z=z*h+z0;
r_vx = [r_vx;[-3,-3,0]];
r_vy = [r_vy;[-3,-3,5]];
red_radius = [red_radius r];
red_Cylcoor = [red_Cylcoor; [x0, y0, z0]];
% surf(x,y,z, 'EdgeColor', 'r', 'FaceColor', 'r', 'FaceAlpha', 0.4)

% transforming the data
red_Cylcoor = cell2mat(red_Cylcoor);
red_radius = cell2mat(red_radius);

%% Creating obstacle
% hold on
% a = -pi : pi/2 : pi;                                % Define Corners
% ph = pi/4;                                          % Define Angular Orientation (‘Phase’)
% x = 0.1*[cos(a+ph); cos(a+ph)]/cos(ph);
% y = 0.1*[sin(a+ph); sin(a+ph)]/sin(ph);
% z = 0.1*[-ones(size(a)); ones(size(a))];                   % Plot Cube
% Block.x = x;
% Block.y = y;
% Block.z = z;

%% plotting the gearbox
Xmean = mean(fv.vertices(:,1));
Ymean = mean(fv.vertices(:,2));
Zmean = mean(fv.vertices(:,3));

hold on
[loc,fitness,best] = geneticalg(fv,100,5,5);
Xcoor = best(1,1)-Xmean+fv.vertices(:,1);
Ycoor = best(1,2)-Ymean+fv.vertices(:,2);
Zcoor = best(1,3)-Zmean+fv.vertices(:,3);
hold off;

%% in/out function -green zone

green_in_out_data = zeros([size(Xcoor),green_numberofcyl]);

for i=1:size(Xcoor)
    pt=[Xcoor(i) Ycoor(i) Zcoor(i)];
    for j= 1:green_numberofcyl
          r = green_radius(j);

          v1 = g_vx(j,:);
          v2 = g_vy(j,:);
          distance = point_to_line(pt,v1,v2);
        
        if distance < r
            green_in_out_data(i,j) = 1;
        
        elseif distance == r
            green_in_out_data(i,j) = 0;
        
        elseif distance > r
            green_in_out_data(i,j) = -1;
        end
    end
end
S = size(green_in_out_data);
A =reshape(green_in_out_data,S);
disp('The given object has ');
OutPts=out_percent(A);
disp(OutPts);
disp(' points outside the cylinders (green-zone) space');

%% in/out function -red zone

red_in_out_data = zeros([size(Xcoor),red_numberofcyl]);

for i=1:size(Xcoor)
    pt=[Xcoor(i) Ycoor(i) Zcoor(i)];
    for j= 1:red_numberofcyl
          r = red_radius(j);

          v1 = r_vx(j,:);
          v2 = r_vy(j,:);
          distance = point_to_line(pt,v1,v2);
        
        if distance < r
            red_in_out_data(i,j) = 1;
        
        elseif distance == r
            red_in_out_data(i,j) = 0;
        
        elseif distance > r
            red_in_out_data(i,j) = -1;
        end
    end
end
S = size(red_in_out_data);
B =reshape(red_in_out_data, S);
disp('The given object has ');
InPts=in_percent(B);
disp(InPts);
disp(' points inside the red-zone cylinders space');

%% functions

function d = point_to_line(pt,v1,v2)
      a = v1 - v2;
      b = pt - v2;

      d = norm(cross(a,b)) / norm(a);
end

function f = out_percent(data)

          j=0;
          k=0;
          l = 0;
          sz = size(data);

           for a=1:sz(1)

                  if any(data(a,:)==1)
                      j=j+1;             
%                   elseif data(a,b)==0
%                       k=k+1;              
                  elseif any(data(a,:)==-1)
                      l=l+1;
                  end               
           end          
          f = l;
end

function g = in_percent(data)

          j=0;
          k=0;
          l = 0;
          sz = size(data);

           for a=1:sz(1)
                  if any(data(a,:)==1)
                      j=j+1;             
%                   elseif data(a,b)==0
%                       k=k+1;              
                  elseif any(data(a,:)==-1)
                      l=l+1;
                  end
             
           end     
          g = j;
end          

function vertex = rotation(V, indice, angle)

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
           vertex = V*Rx;
    end
    if(indice==2)
           vertex = V*Ry;
    end
    if(indice==3)
           vertex = V*Rz;
    end
end 

function [x,d,best] = geneticalg(fv,pop,iter,ulim)
%Create first sample population and evaluate fitness
x = (ulim).*rand(pop,3);
x2 = (ulim).*rand(pop,3);

d = zeros(pop,1);
d2 = zeros(pop,1);
for i=1:pop
    d(i,1) = norm(x(i,:));
    d2(i,1) = norm(x2(i,:));
end
Xmean = mean(fv.vertices(:,1));
Ymean = mean(fv.vertices(:,2));
Zmean = mean(fv.vertices(:,3));

%Plot best of each generation
[srt,I]=sort(d);
best = x(I(1),:);
[srt,I]=sort(d2);
best2= x2(I(1),:);
Xcoor = best(1,1)-Xmean+fv.vertices(:,1);
Ycoor = best(1,2)-Ymean+fv.vertices(:,2);
Zcoor = best(1,3)-Zmean+fv.vertices(:,3);
shaftx = [max(Xcoor) min(Xcoor)];
shafty = [max(Ycoor) min(Ycoor)];
shaftz = [max(Zcoor) min(Zcoor)];
hold on
plot3(shaftx', shafty', shaftz','LineWidth',2)
plot3(Xcoor,Ycoor,Zcoor)
pause(0.5)

for g = 2:iter
    for k = 1: 3: pop
        %Selection
        sel_coeff = mean(d);
        sel_coeff2 = mean(d2);
        parent = zeros(3,3);
        parent2 = zeros(3,3);
        for j=1:3
            rand_sel = 1 + round(pop.*rand());
            for i = rand_sel:pop
                if d(i,1)<sel_coeff
                    parent(j,1) = x(i,1);
                    parent(j,2) = x(i,2);
                    parent(j,3) = x(i,3);
                    break
                end
            end
            rand_sel = 1 + round(pop.*rand());
            for i = rand_sel:pop
                if d2(i,1)<sel_coeff2
                    parent2(j,1) = x2(i,1);
                    parent2(j,2) = x2(i,2);
                    parent2(j,3) = x2(i,3);
                    break
                end
            end
        end

        %Crossover
        child = [parent(1,1) parent(2,2) parent(3,3);parent(1,2) parent(2,3) parent(3,1);parent(1,3) parent(2,1) parent(3,2)];
        child2 = [parent2(1,1) parent2(2,2) parent2(3,3);parent2(1,2) parent2(2,3) parent2(3,1);parent2(1,3) parent2(2,1) parent2(3,2)];
        %Mutation
        for i = 1:3
            for j = 1:3
                mut_prob = rand();
                if mut_prob<=0.05
                    child(i,j) = child(i,j) - (sel_coeff*0.50);
                    if child(i,j)<0
                        child(i,j) = abs(child(i,j));
                    end
                end
                mut_prob = rand();
                if mut_prob<=0.05
                    child2(i,j) = child2(i,j) - (sel_coeff*0.50);
                    if child2(i,j)<0
                        child2(i,j) = abs(child2(i,j));
                    end
                end
            end
        end
        %Write new generation
        %Gearbox
        x(k,1) = child(1,1);
        x(k,2) = child(1,2);
        x(k,3) = child(1,3);
        x(k+1,1) = child(2,1);
        x(k+1,2) = child(2,2);
        x(k+1,3) = child(2,3);
        x(k+2,1) = child(2,1);
        x(k+2,2) = child(2,2);
        x(k+2,3) = child(2,3);
        
        %Electric Component #1
        x2(k,1) = child2(1,1);
        x2(k,2) = child2(1,2);
        x2(k,3) = child2(1,3);
        x2(k+1,1) = child2(2,1);
        x2(k+1,2) = child2(2,2);
        x2(k+1,3) = child2(2,3);
        x2(k+2,1) = child2(2,1);
        x2(k+2,2) = child2(2,2);
        x2(k+2,3) = child2(2,3);
    end
    for i=1:pop
    d(i,1) = norm(x(i,:));
    end
    [srt,I]=sort(d);
    best = x(I(1),:);
    for i=1:pop
    d2(i,1) = norm(x2(i,:));
    end
    [srt,I]=sort(d2);
    best2 = x2(I(1),:);
    %Plot best of generation
    hold on
    Xcoor = best(1,1)-Xmean+fv.vertices(:,1);
    Ycoor = best(1,2)-Ymean+fv.vertices(:,2);
    Zcoor = best(1,3)-Zmean+fv.vertices(:,3);
    plot3(Xcoor,Ycoor,Zcoor)
    pause(0.5)

    
end
end