FV=stlread("F:\TUM\2nd_Semester\Software Lab\Matlab Files\Gearbox.stl");
V = FV.Points;
F = FV.ConnectivityList;
N = faceNormal(FV)
%N = patchnormals(FV)

%P(:,1) = P(:,1) + 0.5;     %add 5 to each vertex's x value
%P(:,2) = P(:,2) + 0.5;     %add 5 to each vertex's x value
%P(:,3) = P(:,3) + 0.5;     %add 5 to each vertex's x value
%V = rotation(V, 1, pi);

fv = struct('faces',F,'vertices',V);


%% Render
% The model is rendered with a PATCH graphics object. We also add some dynamic
% lighting, and adjust the material properties to change the specular
% highlighting.

patch(fv,'FaceColor',       [0.8 0.8 1.0], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15);

xlabel('X')
ylabel('Y')
zlabel('Z')

% Add a camera light, and tone down the specular highlighting
camlight('left');
material('default');

% Fix the axes scaling, and set a nice view angle
axis('image');
view([-135 35]);

hold on
%plot3(N(1:1216),N(1217:(1216*2)),N(2433:(1216*3)))
%plot3(N(1:640),N(641:(640*2)),N(1281:(640*3)))

%% 

ID = [];
for i = 1:length(V)
    ID(i) = i;
end

for i = 1:length(V)
    scatter3(V(i),V(i+length(V)),V(i+2*length(V)),5);
    text(V(i),V(i+length(V)),V(i+2*length(V)),string(ID(i)));
end 
    
trisurf(FV)











































%%
function vertex = rotation(V, indice, angle)
%     [V, F]=stlread("isopoison.stl"); 
%     angle=-pi/2;
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

function N = patchnormals(FV)
%Vertex normals of a triangulated mesh, area weighted, left-hand-rule
% N = patchnormals(FV) -struct with fields, faces Nx3 and vertices Mx3
%N: vertex normals as Mx3

%face corners index
A = FV.ConnectivityList(:,1);
B = FV.ConnectivityList(:,2);
C = FV.ConnectivityList(:,3);

%face normals
n = cross(FV.Points(A,:)-FV.Points(B,:),FV.Points(C,:)-FV.Points(A,:)); %area weighted

%vertice normals
N = zeros(size(FV.Points)); %init vertix normals
for i = 1:size(FV.ConnectivityList,1) %step through faces (a vertex can be reference any number of times)
    N(A(i),:) = N(A(i),:)+n(i,:); %sum face normals
    N(B(i),:) = N(B(i),:)+n(i,:);
    N(C(i),:) = N(C(i),:)+n(i,:);
end 
end
