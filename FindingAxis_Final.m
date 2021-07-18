FV=stlread("F:\TUM\2nd_Semester\Software Lab\Matlab Files\17.7\Gearbox1.STL");
V = FV.Points;
F = FV.ConnectivityList;
N = faceNormal(FV);

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

%% Obtain the exteme values in X,y and Z coorfinates

MaxX_Abs = -1e10; MinX_Abs = 1e10;
MaxY_Abs = -1e10; MinY_Abs = 1e10;
MaxZ_Abs = -1e10; MinZ_Abs = 1e10;


for i = 1:length(V)
    if V(i) > MaxX_Abs
        MaxX_Abs = V(i);
    end
    
    if V(i+length(V)) > MaxY_Abs
        MaxY_Abs = V(i+length(V));
    end
    
    if V(i + 2*length(V)) > MaxZ_Abs
        MaxZ_Abs = V(i + 2*length(V));
    end
    
    if V(i) < MinX_Abs
        MinX_Abs = V(i);
    end
    
    if V(i+length(V)) < MinY_Abs
        MinY_Abs = V(i+length(V));
    end
    
    if V(i + 2*length(V)) < MinZ_Abs
        MinZ_Abs = V(i + 2*length(V));
    end
end


%% Numbers all the vertices and prints the number of every vertix on the plot


ID = [];
for i = 1:length(V)
    ID(i) = i;
    %scatter3(V(i),V(i+length(V)),V(i+2*length(V)),5);
    %text(V(i),V(i+length(V)),V(i+2*length(V)),string(ID(i)));
end

%trisurf(FV)

%% Identifies the faces that are on the curved surface and stores their ID's to a vector "CounterN"
tolerance = 1e-4;

L_Normal = length(N);
L_Vertices = length(V);

CounterN = [];
j = 1;

for i = 1 : L_Normal
    %if ((N(i) == 0) && (N(i+L_Normal) ~= 0) && (N(i+2*L_Normal) ~= 0))||((N(i) ~= 0) && (N(i+L_Normal) ~= 0) && (N(i+2*L_Normal) == 0)) || ((N(i) ~= 0) && (N(i+L_Normal) ~= 0) && (N(i+2*L_Normal) ~= 0))
        
    if ((abs(N(i)) < tolerance) && (abs(N(i+L_Normal)) > tolerance) && (abs(N(i+2*L_Normal)) > tolerance))||((abs(N(i))  > tolerance) && (abs(N(i+L_Normal))  > tolerance) && (abs(N(i+2*L_Normal)) < tolerance)) || ((abs(N(i))  > tolerance) && (abs(N(i+L_Normal)) < tolerance) && (abs(N(i+2*L_Normal))  > tolerance))
        
        CounterN(j)= i;

        j = j+1;
    end
end

%%

counter1 = 0 ;
CylinderNr = 0;
LargestNodeID = 1;

Radius = [];
MaxX = []; MinX = [];
MaxY = []; MinY = [];
MaxZ = []; MinZ = [];
Dx =[]; Dy = []; Dz = [];
CenterX = [];
CenterY = [];
CenterZ = [];
CylinderHeight = [];

while (counter1 < L_Vertices)
    
    counter1 = counter1+1;
    check = 1;
    
    CylinderPoints = [];
    
    CylinderPoints(end+1) =  F(CounterN(counter1));
    CylinderPoints(end+1) =  F(CounterN(counter1) + L_Normal);
    CylinderPoints(end+1) =  F(CounterN(counter1) + 2*L_Normal);
   
    counter1 = counter1 + 2;
    
    while (check ~=0)
       counter2 = counter1; 

        for i = 1:length(CounterN)
               if ((ismember(F(CounterN(i)), CylinderPoints)) || (ismember(F(CounterN(i)+ L_Normal), CylinderPoints)) || (ismember(F(CounterN(i) + 2*L_Normal), CylinderPoints)))

                   if (ismember(F(CounterN(i)), CylinderPoints))
                       %do nothing;
                   else
                       CylinderPoints(end+1) =  F(CounterN(i));
                       counter1 = counter1+1;
                   end

                   if (ismember(F(CounterN(i)+ L_Normal), CylinderPoints))
                       %do nothing;

                   else
                       CylinderPoints(end+1) =  F(CounterN(i) + L_Normal);
                       counter1 = counter1+1;
                   end

                   if (ismember(F(CounterN(i) + 2*L_Normal), CylinderPoints))
                       %do nothing;
                   else
                       CylinderPoints(end+1) =  F(CounterN(i) + 2*L_Normal);
                       counter1 = counter1+1;
                   end
               end
        end 

        if (counter1 == counter2)
            check = 0;
            LargestNodeID = max(CylinderPoints);
            CylinderNr = CylinderNr+1;
        end      
    end

    CylindersX = [];
    CylindersY = [];
    CylindersZ = [];

    for i = 1:length(CylinderPoints)

        CylindersX(end+1) = V(CylinderPoints(i)); 
        CylindersY(end+1) = V(CylinderPoints(i)+ L_Vertices); 
        CylindersZ(end+1) = V(CylinderPoints(i)+ 2*L_Vertices);   

    end

    MaxX(end+1) = max(CylindersX);
    MinX(end+1) = min(CylindersX);
    MaxY(end+1) = max(CylindersY);
    MinY(end+1) = min(CylindersY);
    MaxZ(end+1) = max(CylindersZ);
    MinZ(end+1) = min(CylindersZ);

    Dx(end+1) = MaxX(CylinderNr)-MinX(CylinderNr);
    Dy(end+1) = MaxY(CylinderNr)-MinY(CylinderNr);
    Dz(end+1) = MaxZ(CylinderNr)-MinZ(CylinderNr);

    if (abs(Dx(CylinderNr)- Dy(CylinderNr)) < tolerance)
        Radius(end+1)= Dx(CylinderNr)/2;
        CenterX(CylinderNr) = MinX(CylinderNr) + Radius(CylinderNr);
        CenterY(CylinderNr) = MinY(CylinderNr) + Radius(CylinderNr);
        CylinderHeight(CylinderNr) = abs(MaxZ(CylinderNr) -  MinZ(CylinderNr));

        %zvector=linspace(MinZ(CylinderNr),MaxZ(CylinderNr));
        zvector=linspace(MinZ_Abs -0.02, MaxZ_Abs + 0.02);

        line(CenterX(CylinderNr)*ones(1,length(zvector)),CenterY(CylinderNr)*ones(1,length(zvector)),zvector);        

    elseif (abs(Dx(CylinderNr)- Dz(CylinderNr)) < tolerance)   
        Radius(end+1)= Dx(CylinderNr)/2;
        CenterX(CylinderNr) = MinX(CylinderNr) + Radius(CylinderNr);
        CenterZ(CylinderNr) = MinZ(CylinderNr) + Radius(CylinderNr);
        CylinderHeight(CylinderNr) = abs(MaxY(CylinderNr) -  MinY(CylinderNr));

        
        %xvector=linspace(MinX(CylinderNr),MaxX(CylinderNr));
        yvector=linspace(MinY_Abs -0.02, MaxY_Abs + 0.02);

        line(CenterX(CylinderNr)*ones(1,length(yvector)),yvector,CenterZ(CylinderNr)*ones(1,length(yvector)));        

    elseif (abs(Dy(CylinderNr)- Dz(CylinderNr)) < tolerance)
        Radius(end+1)= Dy(CylinderNr)/2;
        CenterY(CylinderNr) = MinY(CylinderNr) + Radius(CylinderNr);
        CenterZ(CylinderNr) = MinZ(CylinderNr) + Radius(CylinderNr);
        CylinderHeight(CylinderNr) = abs(MaxX(CylinderNr) -  MinX(CylinderNr));

        %xvector=linspace(MinX(CylinderNr),MaxX(CylinderNr));
        xvector=linspace(MinX_Abs -0.02, MaxX_Abs + 0.02);

        line(xvector,CenterY(CylinderNr)*ones(1,length(xvector)),CenterZ(CylinderNr)*ones(1,length(xvector)));
    end
end

%%
ElectricMotorRadius = max(Radius);
MotorIndex = find(Radius == ElectricMotorRadius);

axesX = [];
axesY = [];
axesZ = [];

distance = 0;

if isempty(CenterX)
   for i = 1:length(CenterY)
      if isempty(axesY)
         axesY(end+1) = CenterY(i); 
         axesZ(end+1) = CenterZ(i);
      else
        if ismembertol(CenterY(i),axesY,tolerance)
            %do nothing
        else
           axesY(end+1)=CenterY(i);
           axesZ(end+1)=CenterZ(i);
        end
      end
   end
    MotorCenterY = CenterY(MotorIndex);
	MotorCenterZ = CenterZ(MotorIndex);
      
    for i = 1:length(axesY)
        if (axesY(i) ~= MotorCenterY)
            if (sqrt(( MotorCenterY-axesY(i))^2 + (MotorCenterZ-axesZ(i))^2) > distance)
               distance =  sqrt(( MotorCenterY-axesY(i))^2 + (MotorCenterZ-axesZ(i))^2);
               j = i;
            end
        end
    end
    
    xvector=linspace(MinX_Abs -0.08, MaxX_Abs + 0.08);
    line(xvector,axesY(j)*ones(1,length(xvector)),axesZ(j)*ones(1,length(xvector)),'color','red');       
end
    
if isempty(CenterY)
   for i = 1:length(CenterX)
      if isempty(axesX)
         axesX(end+1) = CenterX(i); 
         axesZ(end+1) = CenterZ(i);
      else
        if ismembertol(CenterX(i),axesX,tolerance)
            %do nothing
        else
           axesX(end+1)=CenterX(i);
           axesZ(end+1)=CenterZ(i);
        end
      end       
   end
    MotorCenterX = CenterX(MotorIndex);
    MotorCenterZ = CenterZ(MotorIndex);
    
    for i = 1:length(axesX)
        if (axesX(i) ~= MotorCenterX)
            if (sqrt(( MotorCenterX-axesX(i))^2 + (MotorCenterZ-axesZ(i))^2) > distance)
               distance =  sqrt(( MotorCenterX-axesX(i))^2 + (MotorCenterZ-axesZ(i))^2);
               j = i;
            end
        end
    end
    
    yvector=linspace(MinY_Abs -0.08, MaxY_Abs + 0.08);
    line(axesX(j)*ones(1,length(yvector)),yvector,axesZ(j)*ones(1,length(yvector)),'color','red');  
end

if isempty(CenterZ)
   for i = 1:length(CenterY)
      if isempty(axesY)
         axesY(end+1) = CenterY(i); 
         axesX(end+1) = CenterX(i);
      else
        if ismembertol(CenterY(i),axesY,tolerance)
            %do nothing
        else
           axesY(end+1)=CenterY(i);
           axesX(end+1)=CenterX(i);
        end
      end       
   end
   MotorCenterX = CenterX(MotorIndex);
   MotorCenterY = CenterY(MotorIndex);
   
   for i = 1:length(axesX)
        if (axesX(i) ~= MotorCenterX)
            if (sqrt(( MotorCenterX-axesX(i))^2 + (MotorCenterY-axesY(i))^2) > distance)
               distance =  sqrt(( MotorCenterX-axesX(i))^2 + (MotorCenterY-axesY(i))^2);
               j = i;
            end
        end
    end
    
    zvector=linspace(MinZ_Abs -0.08, MaxZ_Abs + 0.08);
    line(axesX(j)*ones(1,length(zvector)),axesY(j)*ones(1,length(zvector)),zvector,'color','red');  
end







