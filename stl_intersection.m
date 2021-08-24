clc
clear

%% Load STL mesh
% Import an STL mesh, returning a PATCH-compatible face-vertex structure
FV1 = stlread('Gearbox+EM.stl');
V1 = FV1.Points;
F1 = FV1.ConnectivityList;
N1 = faceNormal(FV1);

fv1 = struct('faces',F1,'vertices',V1);
%% Render

% patch(fv1,'FaceColor',       [0.8 0.8 1.0], ...
%          'EdgeColor',       'none',        ...
%          'FaceLighting',    'gouraud',     ...
%          'AmbientStrength', 0.15);
% 
% % Add a camera light, and tone down the specular highlighting
% camlight('headlight');
% material('dull');
% 
% % Fix the axes scaling, and set a nice view angle
% axis('image');
% % view([-135 35]);
% view([45 35]);

scale=10.0;
fv1.vertices=fv1.vertices*scale;
fv1.vertices = rotation(fv1.vertices, 2, pi/2);

% Xcoor1 = fv1.vertices(:,1);
% Ycoor1 = 22+fv1.vertices(:,2);
% Zcoor1 = fv1.vertices(:,3);
fv1.vertices(:,1)= fv1.vertices(:,1);
fv1.vertices(:,2)= 22+fv1.vertices(:,2);
fv1.vertices(:,3) = fv1.vertices(:,3);
figure
% surf(Xcoor1,Ycoor1,Zcoor1);

% data = FV1;
TR = triangulation(fv1.faces,fv1.vertices) ;
trimesh(TR,'FaceColor','none','EdgeColor','blue')

hold on
% figure(2)
% patch(Xcoor,Ycoor,Zcoor);

%% %% Load STL mesh
% Import an STL mesh, returning a PATCH-compatible face-vertex structure

FV2 = stlread('Cube.stl');
V2 = FV2.Points;
F2 = FV2.ConnectivityList;
N2 = faceNormal(FV2);

fv2 = struct('faces',F2,'vertices',V2);
%% Render

% p2=patch(fv2,'FaceColor',       [0.8 0.8 1.0], ...
%          'EdgeColor',       'none',        ...
%          'FaceLighting',    'gouraud',     ...
%          'AmbientStrength', 0.15);
% 
% % Add a camera light, and tone down the specular highlighting
% camlight('headlight');
% material('dull');
% 
% % Fix the axes scaling, and set a nice view angle
% axis('image');
% % view([-135 35]);
% view([45 35]);
% 
% scale=0.2;
% fv2.vertices=fv2.vertices*scale;
scale=0.003;
fv2.vertices=fv2.vertices*scale;
hold on
% Xcoor2 =  17+fv2.vertices(:,1);
% Ycoor2 = 6+fv2.vertices(:,2);
% Zcoor2 = -2+fv2.vertices(:,3);
fv2.vertices(:,1) =  fv2.vertices(:,1);
fv2.vertices(:,2) = 2+fv2.vertices(:,2);
fv2.vertices(:,3) = -2+fv2.vertices(:,3);
TR2 = triangulation(fv2.faces,fv2.vertices) ;
trimesh(TR2,'FaceColor','none','EdgeColor','green')
hold off
% figure

% hold off
% figure(2)
% patch(Xcoor,Ycoor,Zcoor);
%% 

% figure(3)
% plot3(Xcoor1,Ycoor1,Zcoor1)
% hold on
% patch(Xcoor2,Ycoor2,Zcoor2)
% hold off


%% 
IN = inpolyhedron(fv1,fv2.vertices);
% disp(IN);
if any(IN(:)==1)
  disp("Objects are colliding");                         
else
  disp("Objects are free from collision");
end
%% Plot a cylinder 
% [T1, T2, T3] = ndgrid(X, Y, Z);
% A = [T1(:), T2(:), T3(:)];

% green zone cylinders
% green_Cylcoor = {};
% green_radius = {};
% centrelinesX = {};
% centrelinesY = {};
% green_numberofcyl = 4;
% 
% % cylinder 1
% figure
% 
% % r = 3.5;
% % [x,y,z] = cylinder(r,100);
% % h = 5;
% % Z = Z*h;
% % mesh(x,y,z);
% % v1 = [0,0,0];
% % v2 = [0,0,5];
% % radius = [radius, r];
% % Cylcoor = [Cylcoor, [x, y, z]];
% %surf(X,Y,Z)
% %plot(X,Y)
% 
% %Base at 
% x0=0;y0=0;z0=0;
% r = 2.5;
% h=5;
% [x,y,z]=cylinder(r,100);
% x=x+x0;
% y=y+y0;
% z=z*h+z0;
% % vx = [0,0,0];
% % vy = [0,0,5];
% g_vx = [0,0,0];
% g_vy = [0,0,5];
% % centrelines(:,:,1) = {v1; v2};
% % centrelinesX(1) = v1;
% % centrelinesY(2) = v2;
% green_radius = [green_radius r];
% green_Cylcoor = [green_Cylcoor; [x0, y0, z0]];
% surf(x,y,z, 'EdgeColor', 'g', 'FaceColor', 'g', 'FaceAlpha', 0.4)
% 
% % cylinder 2
% hold on;
% %Base at 
% x0=-5;y0=0;z0=0;
% r = 3.5;
% h=5;
% [x,y,z]=cylinder(r,100);
% x=x+x0;
% y=y+y0;
% z=z*h+z0;
% % vx{2} = [-5,0,0];
% % vy{2} = [-5,0,5];
% g_vx = [g_vx;[-5,0,0]];
% g_vy = [g_vy;[-5,0,5]];
% % centrelines(:,:,2) = {v1; v2};
% % centrelines(1,2) = v1;
% % centrelines(2,2) = v2;
% green_radius = [green_radius r];
% green_Cylcoor = [green_Cylcoor; [x0, y0, z0]];
% surf(x,y,z, 'EdgeColor', 'g', 'FaceColor', 'g', 'FaceAlpha', 0.4)
% 
% 
% % cylinder 3
% hold on;
% %Base at 
% x0=-5;y0=-5;z0=0;
% r = 3.5;
% h=5;
% [x,y,z]=cylinder(r,100);
% x=x+x0;
% y=y+y0;
% z=z*h+z0;
% % vx{3} = [-5,-5,0];
% % vy{3} = [-5,-5,5];
% g_vx = [g_vx;[-5,-5,0]];
% g_vy = [g_vy;[-5,-5,5]];
% % centrelines(:,:,3) = {v1; v2};
% % centrelines(1,3) = v1;
% % centrelines(2,3) = v2;
% green_radius = [green_radius r];
% green_Cylcoor = [green_Cylcoor; [x0, y0, z0]];
% surf(x,y,z, 'EdgeColor', 'g', 'FaceColor', 'g', 'FaceAlpha', 0.4)
% 
% % cylinder 4
% hold on;
% %Base at 
% x0=0;y0=-5;z0=0;
% r = 3.5;
% h=5;
% [x,y,z]=cylinder(r,100);
% x=x+x0;
% y=y+y0;
% z=z*h+z0;
% % vx{4} = [0,-5,0];
% % vy{4} = [0,-5,5];
% g_vx = [g_vx;[0,-5,0]];
% g_vy = [g_vy;[0,-5,5]];
% % centrelines(:,:,4) = {v1; v2};
% % centrelines(1,4) = v1;
% % centrelines(2,4) = v2;
% green_radius = [green_radius r];
% green_Cylcoor = [green_Cylcoor; [x0, y0, z0]];
% surf(x,y,z, 'EdgeColor', 'g', 'FaceColor', 'g', 'FaceAlpha', 0.4)
% 
% 
% green_Cylcoor = cell2mat(green_Cylcoor);
% green_radius = cell2mat(green_radius);
% % A = [T1(:), T2(:), T3(:)];
% % [x,y,z]=cylinder([0,0.5,0.5,0],50,[1 1 0]);
% % z([1,2],:)=0;
% % z([3,4],:)=5;
% % figure(1)
% % hm = mesh(x,y,z);
% % rotate(hm, [0 1 0], 90)
% % axis equal
% % figure(2)
% % hm = mesh(x,y,z);
% % rotate(hm, [0 1 0], 90)
% % axis equal
% % figure(3)
% % hm = mesh(x,y,z);
% % rotate(hm, [1 1 0], 45)
% % axis equal
% 
% % %% Translation
% % 
% % objf = double.empty;
% % n=min(Ycoor);
% % % m=min(Zcoor);
% % while n > -1.0
% % %     while m > 5
% %         
% %         Ycoor = Ycoor - 0.3;
% % %         Zcoor = Zcoor - 0.5;
% %         n=n-0.3;
% % %         m=m-0.5;
% %         figure
% %         hm = mesh(x,y,z);
% %         rotate(hm, [0 1 0], 90)
% %         axis equal
% %         patch('XData',Xcoor,'YData',Ycoor,'ZData',Zcoor);
% % %     end
% % end
% %% red zone cyl
% red_Cylcoor = {};
% red_radius = {};
% % centrelinesX = {};
% % centrelinesY = {};
% red_numberofcyl = 2;
% 
% hold on;
% %Base at 
% x0=0;y0=-3;z0=0;
% r = 1.5;
% h=5;
% [x,y,z]=cylinder(r,100);
% x=x+x0;
% y=y+y0;
% z=z*h+z0;
% % vx{4} = [0,-5,0];
% % vy{4} = [0,-5,5];
% % vx = [vx;[0,-5,0]];
% % vy = [vy;[0,-5,5]];
% r_vx = [0,-3,0];
% r_vy = [0,-3,5];
% % centrelines(:,:,4) = {v1; v2};
% % centrelines(1,4) = v1;
% % centrelines(2,4) = v2;
% red_radius = [red_radius r];
% red_Cylcoor =  [red_Cylcoor; [x0, y0, z0]];
% % mesh(x,y,z)
% surf(x,y,z, 'EdgeColor', 'r', 'FaceColor', 'r', 'FaceAlpha', 0.4)
% 
% hold on;
% %Base at 
% x0=-3;y0=-3;z0=0;
% r = 1;
% h=5;
% [x,y,z]=cylinder(r,100);
% x=x+x0;
% y=y+y0;
% z=z*h+z0;
% % vx{4} = [0,-5,0];
% % vy{4} = [0,-5,5];
% % vx = [vx;[0,-5,0]];
% % vy = [vy;[0,-5,5]];
% r_vx = [r_vx;[-3,-3,0]];
% r_vy = [r_vy;[-3,-3,5]];
% % centrelines(:,:,4) = {v1; v2};
% % centrelines(1,4) = v1;
% % centrelines(2,4) = v2;
% red_radius = [red_radius r];
% red_Cylcoor = [red_Cylcoor; [x0, y0, z0]];
% % mesh(x,y,z)
% surf(x,y,z, 'EdgeColor', 'r', 'FaceColor', 'r', 'FaceAlpha', 0.4)
% 
% red_Cylcoor = cell2mat(red_Cylcoor);
% red_radius = cell2mat(red_radius);
% 
% %% 
% 
% 
% hold on
% Xcoor1 = -2.+fv1.vertices(:,1);
% Ycoor1 = -1.+fv1.vertices(:,2);
% Zcoor1 = 4.+fv1.vertices(:,3);
% % polyshape([Xcoor Ycoor Zcoor]);
% plot3(Xcoor1,Ycoor1,Zcoor1)
% % plot3(polyshape)
% hold on;
% % hm = mesh(X,Y,Z);
% % rotate(hm, [0 1 0], 90)
% % axis equal
% hold on
% Xcoor2 = fv2.vertices(:,1);
% Ycoor2 = fv2.vertices(:,2);
% Zcoor2 = fv2.vertices(:,3);
% % polyshape([Xcoor Ycoor Zcoor]);
% plot3(Xcoor1,Ycoor1,Zcoor1)
% % plot3(polyshape)
% hold on;
% 
% %% in/out function
% % v1 = [0,0,0];
% % v2 = [0,0,5];
% % in_out_data = zeros(size(Xcoor));
% green_in_out_data = zeros([size(Xcoor1),green_numberofcyl]);
% 
% for i=1:size(Xcoor1)
%     pt=[Xcoor1(i) Ycoor1(i) Zcoor1(i)];
%     for j= 1:green_numberofcyl
%         r = green_radius(j);
% %         v1 = centrelines(1,j);
% %         v2 = centrelines(2,j);
%           v1 = g_vx(j,:);
%           v2 = g_vy(j,:);
%         distance = point_to_line(pt,v1,v2);
%         
%         if distance < r
%             green_in_out_data(i,j) = 1;
%         
%         elseif distance == r
%             green_in_out_data(i,j) = 0;
%         
%         elseif distance > r
%             green_in_out_data(i,j) = -1;
%         end
%     end
% end
% S = size(green_in_out_data);
% A =reshape(green_in_out_data, [S(1),S(3)]);
% disp('The given object has ');
% OutPts=out_percent(A);
% disp(OutPts);
% disp(' points outside the cylinders (green-zone) space');
% 
% %% 
% 
% red_in_out_data = zeros([size(Xcoor1),red_numberofcyl]);
% 
% for i=1:size(Xcoor1)
%     pt=[Xcoor1(i) Ycoor1(i) Zcoor1(i)];
%     for j= 1:red_numberofcyl
%         r = red_radius(j);
% %         v1 = centrelines(1,j);
% %         v2 = centrelines(2,j);
%           v1 = r_vx(j,:);
%           v2 = r_vy(j,:);
%         distance = point_to_line(pt,v1,v2);
%         
%         if distance < r
%             red_in_out_data(i,j) = 1;
%         
%         elseif distance == r
%             red_in_out_data(i,j) = 0;
%         
%         elseif distance > r
%             red_in_out_data(i,j) = -1;
%         end
%     end
% end
% S = size(red_in_out_data);
% B =reshape(red_in_out_data, [S(1),S(3)]);
% disp('The given object has ');
% InPts=in_percent(B);
% disp(InPts);
% disp(' points inside the red-zone cylinders space');
% 
% %% intersection 
% 
% % for i = 1:size(F1)
% %     vertex1=F1(i,1);
% %     vertex2=F1(i,1);
% %     vertex3=F1(i,1);
% %     
% %     
% %     
% %     
% %     shp = alphashape([V1(vertex1,1) V1(vertex2,1) V1(vertex3,1)],[V1(vertex1,2) V1(vertex2,2) V1(vertex3,2)],[V1(vertex1,3) V1(vertex2,3) V1(vertex3,3)]);
% 
% IN = inpolyhedron(fv1,fv2.vertices);
% 
% %     for j=1:size(V2)
% %         
% %         
% %         IN = inpolyhedron(fv1,V2(j,1),V2(j,2),V2(j,1));
% %         
% %     end
%      
% % end




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

% 
% function d = point_to_line(pt,v1,v2)
%       a = v1 - v2;
%       b = pt - v2;
% %         a = centrelines(1,1,j) - centrelines(2,1,j);
% %         b = pt - centrelines(2,1,j);
%       d = norm(cross(a,b)) / norm(a);
% end
% function f = out_percent(data)
% %       j=0;
% %       k=0;
% %       l = 0;
% %       sz = size(data);
% %       for a=1:sz(2)
% %           for i=1:sz(1)
% %               if data(i)==1
% %                   j=j+1;             
% %               elseif data(i)==0
% %                   k=k+1;              
% %               elseif data(i)==-1
% %                   l=l+1;
% %               end
% %           end
% %       end
% %       f = j;
% %       f= mrdivide(j,size(m));
%           j=0;
%           k=0;
%           l = 0;
%           sz = size(data);
% %           for a=1:sz(1)
% %               for b=1:sz(2)
% %                   if data(a,b)==1
% %                       j=j+1;             
% %                   elseif data(a,b)==0
% %                       k=k+1;              
% %                   elseif data(a,b)==-1
% %                       l=l+1;
% %                   end
% %               end
% %            end
%            for a=1:sz(1)
% %               for b=1:sz(2)
%                   if any(data(a,:)==1)
%                       j=j+1;             
% %                   elseif data(a,b)==0
% %                       k=k+1;              
%                   elseif any(data(a,:)==-1)
%                       l=l+1;
%                   end
% %               end
%            end
% %           
%           f = l;
% 
% end
% 
% function g = in_percent(data)
% %       j=0;
% %       k=0;
% %       l = 0;
% %       sz = size(data);
% %       for a=1:sz(2)
% %           for i=1:sz(1)
% %               if data(i)==1
% %                   j=j+1;             
% %               elseif data(i)==0
% %                   k=k+1;              
% %               elseif data(i)==-1
% %                   l=l+1;
% %               end
% %           end
% %       end
% %       f = j;
% %       f= mrdivide(j,size(m));
%           j=0;
%           k=0;
%           l = 0;
%           sz = size(data);
% %           for a=1:sz(1)
% %               for b=1:sz(2)
% %                   if data(a,b)==1
% %                       j=j+1;             
% %                   elseif data(a,b)==0
% %                       k=k+1;              
% %                   elseif data(a,b)==-1
% %                       l=l+1;
% %                   end
% %               end
% %            end
%            for a=1:sz(1)
% %               for b=1:sz(2)
%                   if any(data(a,:)==1)
%                       j=j+1;             
% %                   elseif data(a,b)==0
% %                       k=k+1;              
%                   elseif any(data(a,:)==-1)
%                       l=l+1;
%                   end
% %               end
%            end
% %           
%           g = j;
% 
% end          
% 
% 
% function vertex = rotation(V, indice, angle)
% 
%     Rz = [ cos(angle), -sin(angle), 0 ;
%           sin(angle), cos(angle), 0 ;
%     0, 0, 1 ];
%     Ry = [ cos(angle), 0, sin(angle) ;
%     0, 1, 0 ;
%           -sin(angle), 0, cos(angle) ];
%     Rx = [ 1, 0, 0 ;
%     0, cos(angle), -sin(angle);
%     0, sin(angle), cos(angle) ];
%     
%     if(indice==1)
%            vertex = V*Rx;
%     end
%     if(indice==2)
%            vertex = V*Ry;
%     end
%     if(indice==3)
%            vertex = V*Rz;
%     end
% end 
% %% 




