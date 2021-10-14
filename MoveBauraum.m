function  bauraum = MoveBauraum(bauraum)

%Defining the limits of the bauraum
% xmax = max(bauraum.vertices(:,1));
xmin = min(bauraum.vertices(:,1));
% ymax = max(bauraum.vertices(:,2));
ymin = min(bauraum.vertices(:,2));
% zmax = max(bauraum.vertices(:,3));
zmin = min(bauraum.vertices(:,3));

%Moving the bauraum to the origin
bauraum_mov = [xmin ymin zmin];
bauraum.vertices = bauraum.vertices - bauraum_mov;

% xmax = xmax-xmin;
% ymax = ymax-ymin;
% zmax = zmax-zmin;
% xmin=0;
% ymin=0;
% zmin=0;

end