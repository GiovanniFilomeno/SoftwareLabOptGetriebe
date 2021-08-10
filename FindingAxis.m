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
    %scatter3(V(i),V(i+length(V)),V(i+2*length(V)),10);
    %text(V(i),V(i+length(V)),V(i+2*length(V)),string(ID(i)));
end

%trisurf(FV)

%% Identifies the faces that are on the curved surface and stores their ID's to a vector "CounterN"
tolerance = 1e-5;

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
CylinderCenters_X_coordinates = [];
CylinderCenters_Y_coordinates = [];
CylinderCenters_Z_coordinates = [];
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
        CylinderCenters_X_coordinates(CylinderNr) = MinX(CylinderNr) + Radius(CylinderNr);
        CylinderCenters_Y_coordinates(CylinderNr) = MinY(CylinderNr) + Radius(CylinderNr);
        CylinderHeight(CylinderNr) = abs(MaxZ(CylinderNr) -  MinZ(CylinderNr));

        %zvector=linspace(MinZ(CylinderNr),MaxZ(CylinderNr));
        zvector=linspace(MinZ_Abs -0.02, MaxZ_Abs + 0.02);

        %line(CenterX(CylinderNr)*ones(1,length(zvector)),CenterY(CylinderNr)*ones(1,length(zvector)),zvector);        

    elseif (abs(Dx(CylinderNr)- Dz(CylinderNr)) < tolerance)   
        Radius(end+1)= Dx(CylinderNr)/2;
        CylinderCenters_X_coordinates(CylinderNr) = MinX(CylinderNr) + Radius(CylinderNr);
        CylinderCenters_Z_coordinates(CylinderNr) = MinZ(CylinderNr) + Radius(CylinderNr);
        CylinderHeight(CylinderNr) = abs(MaxY(CylinderNr) -  MinY(CylinderNr));

        
        %xvector=linspace(MinX(CylinderNr),MaxX(CylinderNr));
        yvector=linspace(MinY_Abs -0.02, MaxY_Abs + 0.02);

        %line(CenterX(CylinderNr)*ones(1,length(yvector)),yvector,CenterZ(CylinderNr)*ones(1,length(yvector)));        

    elseif (abs(Dy(CylinderNr)- Dz(CylinderNr)) < tolerance)
        Radius(end+1)= Dy(CylinderNr)/2;
        CylinderCenters_Y_coordinates(CylinderNr) = MinY(CylinderNr) + Radius(CylinderNr);
        CylinderCenters_Z_coordinates(CylinderNr) = MinZ(CylinderNr) + Radius(CylinderNr);
        CylinderHeight(CylinderNr) = abs(MaxX(CylinderNr) -  MinX(CylinderNr));

        %xvector=linspace(MinX(CylinderNr),MaxX(CylinderNr));
        xvector=linspace(MinX_Abs -0.02, MaxX_Abs + 0.02);

        %line(xvector,CenterY(CylinderNr)*ones(1,length(xvector)),CenterZ(CylinderNr)*ones(1,length(xvector)));
    end
end

%%
ElectricMotorRadius = max(Radius);
MotorIndex = find(Radius == ElectricMotorRadius);

FinalAxesX = [];
FinalAxesY = [];
FinalAxesZ = [];

distance = 0;

if isempty(CylinderCenters_X_coordinates)
   for i = 1:length(CylinderCenters_Y_coordinates)
      if isempty(FinalAxesY)
         FinalAxesY(end+1) = CylinderCenters_Y_coordinates(i); 
         FinalAxesZ(end+1) = CylinderCenters_Z_coordinates(i);
      else
        if ismembertol(CylinderCenters_Y_coordinates(i),FinalAxesY,tolerance)
            %do nothing
        else
           FinalAxesY(end+1)=CylinderCenters_Y_coordinates(i);
           FinalAxesZ(end+1)=CylinderCenters_Z_coordinates(i);
        end
      end
   end
    MotorCenterY = CylinderCenters_Y_coordinates(MotorIndex);
	MotorCenterZ = CylinderCenters_Z_coordinates(MotorIndex);   

    StartingAxisY = MotorCenterY;
    StartingAxisZ = MotorCenterZ;
    LastCylinderID = 0;

    for i = 1:length(FinalAxesY)
        for j = 1:length(CylinderCenters_Y_coordinates) %Loop over all cylinders and use following IF-STATEMENT to find the cylinders that lie on the same axis as the StartingAxis
            if ((abs(CylinderCenters_Y_coordinates(j) - StartingAxisY) < tolerance)&&(j ~= MotorIndex) && (j~= LastCylinderID))
                for k = 1:length(CylinderCenters_Y_coordinates) %Loop over all other cylinders lying on diferent axes and searching for meshing cylinders(gears)
                    if ((abs(CylinderCenters_Y_coordinates(k) - StartingAxisY)) > tolerance) && (k ~= MotorIndex) % && (j~=k)
                        DistanceBetweenAxes = sqrt(( StartingAxisY - CylinderCenters_Y_coordinates(k))^2 + (StartingAxisZ-CylinderCenters_Z_coordinates(k))^2);
                            if ((((Radius(j)+ Radius(k))- DistanceBetweenAxes ) > 0) && (abs(MinX(k)-MinX(j))<tolerance)&& (k~= LastCylinderID)) %&& (( DistanceBetweenAxes - (Radius(j)+ Radius(k))) > -0.01) 
                                fprintf('main cylinder is number %i and other cylinder is number %i \n',j,k)
                                StartingAxisY = CylinderCenters_Y_coordinates(k);
                                StartingAxisZ = CylinderCenters_Z_coordinates(k);
                                LastCylinderID = k;
                                 break;
                             end
                    end
                end
            end 
        end
    end

    xvector = linspace(MinX_Abs -0.08, MaxX_Abs + 0.08);
    line(xvector,CylinderCenters_Y_coordinates(LastCylinderID)*ones(1,length(xvector)),CylinderCenters_Z_coordinates(LastCylinderID)*ones(1,length(xvector)),'color','red');    
end    

if isempty(CylinderCenters_Y_coordinates)
   for i = 1:length(CylinderCenters_X_coordinates)
      if isempty(FinalAxesX)
         FinalAxesX(end+1) = CylinderCenters_X_coordinates(i); 
         FinalAxesZ(end+1) = CylinderCenters_Z_coordinates(i);
      else
        if ismembertol(CylinderCenters_X_coordinates(i),FinalAxesX,tolerance)
            %do nothing
        else
           FinalAxesX(end+1)=CylinderCenters_X_coordinates(i);
           FinalAxesZ(end+1)=CylinderCenters_Z_coordinates(i);
        end
      end       
   end
    MotorCenterX = CylinderCenters_X_coordinates(MotorIndex);
    MotorCenterZ = CylinderCenters_Z_coordinates(MotorIndex);
    
    StartingAxisX = MotorCenterX;
    StartingAxisZ = MotorCenterZ;
    LastCylinderID = 0;
    
    for i = 1:length(FinalAxesX)
        for j = 1:length(CylinderCenters_X_coordinates) %Loop over all cylinders and use following IF-STATEMENT to find the cylinders that lie on the same axis as the StartingAxis
            if ((abs(CylinderCenters_X_coordinates(j) - StartingAxisX) < tolerance)&&(j ~= MotorIndex) && (j~= LastCylinderID))
                for k = 1:length(CylinderCenters_X_coordinates) %Loop over all other cylinders lying on diferent axes and searching for meshing cylinders(gears)
                    if ((abs(CylinderCenters_X_coordinates(k) - StartingAxisX)) > tolerance) && (k ~= MotorIndex)
                        DistanceBetweenAxes = sqrt(( StartingAxisX - CylinderCenters_X_coordinates(k))^2 + (StartingAxisZ-CylinderCenters_Z_coordinates(k))^2);
                            if ((((Radius(j)+ Radius(k))- DistanceBetweenAxes ) > 0) && (abs(MinY(k)-MinY(j))<tolerance)&& (k~= LastCylinderID))
                                fprintf('main cylinder is number %i and other cylinder is number %i \n',j,k)
                                StartingAxisX = CylinderCenters_X_coordinates(k);
                                StartingAxisZ = CylinderCenters_Z_coordinates(k);
                                LastCylinderID = k;
                                break;
                             end
                    end
                end
            end 
        end
    end
    yvector = linspace(MinY_Abs -0.08, MaxY_Abs + 0.08);
    line(CylinderCenters_X_coordinates(LastCylinderID)*ones(1,length(yvector)),yvector,CylinderCenters_Z_coordinates(LastCylinderID)*ones(1,length(yvector)),'color','red');     
end

if isempty(CylinderCenters_Z_coordinates)
   for i = 1:length(CylinderCenters_Y_coordinates)
      if isempty(FinalAxesY)
         FinalAxesY(end+1) = CylinderCenters_Y_coordinates(i); 
         FinalAxesX(end+1) = CylinderCenters_X_coordinates(i);
      else
        if ismembertol(CylinderCenters_Y_coordinates(i),FinalAxesY,tolerance)
            %do nothing
        else
           FinalAxesY(end+1)=CylinderCenters_Y_coordinates(i);
           FinalAxesX(end+1)=CylinderCenters_X_coordinates(i);
        end
      end       
   end
   MotorCenterX = CylinderCenters_X_coordinates(MotorIndex);
   MotorCenterY = CylinderCenters_Y_coordinates(MotorIndex);

    StartingAxisX = MotorCenterX;
    StartingAxisY = MotorCenterY;
    LastCylinderID = 0;
    
for i = 1:length(FinalAxesX)
    for j = 1:length(CylinderCenters_X_coordinates) %Loop over all cylinders and use following IF-STATEMENT to find the cylinders that lie on the same axis as the StartingAxis
        if ((abs(CylinderCenters_X_coordinates(j) - StartingAxisX) < tolerance)&&(j ~= MotorIndex) && (j~= LastCylinderID))
            for k = 1:length(CylinderCenters_X_coordinates) %Loop over all other cylinders lying on diferent axes and searching for meshing cylinders(gears)
                if ((abs(CylinderCenters_X_coordinates(k) - StartingAxisX)) > tolerance) && (k ~= MotorIndex) 
                    DistanceBetweenAxes = sqrt(( StartingAxisX - CylinderCenters_X_coordinates(k))^2 + (StartingAxisY-CylinderCenters_Y_coordinates(k))^2);
                        if ((((Radius(j)+ Radius(k))- DistanceBetweenAxes ) > 0) && (abs(MinZ(k)-MinZ(j))<tolerance)&& (k~= LastCylinderID)) 
                            fprintf('main cylinder is number %i and other cylinder is number %i \n',j,k)
                            StartingAxisX = CylinderCenters_X_coordinates(k);
                            StartingAxisY = CylinderCenters_Y_coordinates(k);
                            LastCylinderID = k;
                            break;
                         end
                end
            end
        end 
    end
end
    zvector = linspace(MinZ_Abs -0.08, MaxZ_Abs + 0.08);
    line(CylinderCenters_X_coordinates(LastCylinderID)*ones(1,length(zvector)),CylinderCenters_Y_coordinates(LastCylinderID)*ones(1,length(zvector)),zvector,'color','red');       
end






