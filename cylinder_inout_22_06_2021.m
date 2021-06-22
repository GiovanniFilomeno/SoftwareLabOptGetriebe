clc
clear

%% Load STL mesh
% Stereolithography (STL) files are a common format for storing mesh data. STL
% meshes are simply a collection of triangular faces. This type of model is very
% suitable for use with MATLAB's PATCH graphics object.

% Import an STL mesh, returning a PATCH-compatible face-vertex structure
fv = stlread('Gearbox+EM.stl');

%% Render
% The model is rendered with a PATCH graphics object. We also add some dynamic
% lighting, and adjust the material properties to change the specular
% highlighting.

patch(fv,'FaceColor',       [0.8 0.8 1.0], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15);

% Add a camera light, and tone down the specular highlighting
camlight('headlight');
material('dull');

% Fix the axes scaling, and set a nice view angle
axis('image');
view([45 90]);

scale=5.0;
fv.vertices=fv.vertices*scale;
fv.vertices = rotation(fv.vertices, 2, pi/2);
%% Cylinders

% green zone cylinders
green_Cylcoor = {};
green_radius = {};
green_numberofcyl = 4;

figure

% cylinder 1
%Base at 
x0=0;y0=0;z0=0;
r = 2.5;
h=5;
[x,y,z]=cylinder(r,100);
x=x+x0;
y=y+y0;
z=z*h+z0;
g_vx = [0,0,0];
g_vy = [0,0,5];
green_radius = [green_radius r];
green_Cylcoor = [green_Cylcoor; [x0, y0, z0]];
surf(x,y,z, 'EdgeColor', 'g', 'FaceColor', 'g', 'FaceAlpha', 0.4)
xlabel('X')
ylabel('Y')
zlabel('Z')
view([45 90]);


% cylinder 2
hold on;
%Base at 
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
surf(x,y,z, 'EdgeColor', 'g', 'FaceColor', 'g', 'FaceAlpha', 0.4)


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
surf(x,y,z, 'EdgeColor', 'g', 'FaceColor', 'g', 'FaceAlpha', 0.4)


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
surf(x,y,z, 'EdgeColor', 'g', 'FaceColor', 'g', 'FaceAlpha', 0.4)

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
surf(x,y,z, 'EdgeColor', 'r', 'FaceColor', 'r', 'FaceAlpha', 0.4)


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
surf(x,y,z, 'EdgeColor', 'r', 'FaceColor', 'r', 'FaceAlpha', 0.4)

% transforming the data
red_Cylcoor = cell2mat(red_Cylcoor);
red_radius = cell2mat(red_radius);

%% plotting the gearbox

hold on
Xcoor = fv.vertices(:,1);
Ycoor = 4.5+fv.vertices(:,2);
Zcoor = fv.vertices(:,3);
plot3(Xcoor,Ycoor,Zcoor)

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
A =reshape(green_in_out_data, [S(1),S(3)]);
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
B =reshape(red_in_out_data, [S(1),S(3)]);
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
