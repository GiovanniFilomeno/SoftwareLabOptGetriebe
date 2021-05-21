%% Load STL mesh
% Stereolithography (STL) files are a common format for storing mesh data. STL
% meshes are simply a collection of triangular faces. This type of model is very
% suitable for use with MATLAB's PATCH graphics object.

% Import an STL mesh, returning a PATCH-compatible face-vertex structure
fv = stlread("Gearbox.stl");

P = fv.vertices;     %access the vertex data from triangulation              
C = fv.faces;     %access the connectivity data from triangulation
%P(:,1) = P(:,1) + 0.5;     %add 5 to each vertex's x value
%P(:,2) = P(:,2) + 0.5;     %add 5 to each vertex's x value
%P(:,3) = P(:,3) + 0.5;     %add 5 to each vertex's x value


P = rotation(P, 1, 1.5);


fv = struct('faces',C,'vertices',P);



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
view([-135 35]);

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
%fv = triangulation(C, P);     %Combine both components back into a triangulation variable