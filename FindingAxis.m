function [EM_start, EM_vector, Anchor,ElectricMotorRadius] = FindingAxis(file)


    if ~exist(file,'file')
        error(['File ''%s'' not found. If the file is not on MATLAB''s path' ...
               ', be sure to specify the full path to the file.'], file);
    end
    %%Reads STL file and obtains vertices, Faces and Normal Vectors within the geometry
    FV = stlread(file);
    %FV = stlread("C:\Users\ahmed\OneDrive\Desktop\Software Lab\Matlab Files\Current\Gearbox1.stl");
    %FV = stlread("C:\Users\ahmed\OneDrive\Desktop\Software Lab\Matlab Files\Current\Gearbox1(Rotated1).stl");
    %FV = stlread("C:\Users\ahmed\OneDrive\Desktop\Software Lab\Matlab Files\Current\Gearbox1(Rotated2).stl");

    V = 1000*FV.Points;  %%Obtains Vertices matrix
    F = FV.ConnectivityList;    %Obtains Faces Matrix
    N = faceNormal(FV); %Obtains Matrix of normal vectors

    fv = struct('faces',F,'vertices',V);    %Creates a structure from the obtaines faces and vertices

    %% Global Variables
    tolerance = 1e-2; %tolerance value used in the whole algorithm
    NumberOfVertices = length(V);
    NumberOfNormalVectors = length(N);

    %% Render
    % The model is rendered with a PATCH graphics object. We also add some dynamic
    % lighting, and adjust the material properties to change the specular
    % highlighting.

%     patch(fv,'FaceColor',       [0.8 0.8 1.0], ...
%              'EdgeColor',       'none',        ...
%              'FaceLighting',    'gouraud',     ...
%              'AmbientStrength', 0.15);
% 
%     xlabel('X')
%     ylabel('Y')
%     zlabel('Z')
% 
%     % Add a camera light, and tone down the specular highlighting
%     camlight('left');
%     material('default');
% 
%     % Fix the axes scaling, and set a nice view angle
%     axis('image');
%     view([-135 35]);

    %% Obtain the extreme values in X,Y and Z coordinates for the whole geometry by looping over all vertices

    MaxX_Abs = -1e10; MinX_Abs = 1e10;  %Maximum X value for whole geometry
    MaxY_Abs = -1e10; MinY_Abs = 1e10;  %Maximum X value for whole geometry
    MaxZ_Abs = -1e10; MinZ_Abs = 1e10;  %Maximum X value for whole geometry

    for i = 1:NumberOfVertices
        if V(i) > MaxX_Abs
            MaxX_Abs = V(i);
        end

        if V(i+NumberOfVertices) > MaxY_Abs
            MaxY_Abs = V(i+NumberOfVertices);
        end

        if V(i + 2*NumberOfVertices) > MaxZ_Abs
            MaxZ_Abs = V(i + 2*NumberOfVertices);
        end

        if V(i) < MinX_Abs
            MinX_Abs = V(i);
        end

        if V(i+NumberOfVertices) < MinY_Abs
            MinY_Abs = V(i+NumberOfVertices);
        end

        if V(i + 2*NumberOfVertices) < MinZ_Abs
            MinZ_Abs = V(i + 2*NumberOfVertices);
        end
    end

    %% Assignes an ID number to all the vertices of the geometry
    % ID = [];
    % for i = 1:NumberOfVertices
    %     ID(i) = i;
    % end

    %% Prints Scatter points and ID Numbers at the locations of all the vertices 
    % for i = 1:NumberOfVertices
    %     scatter3(V(i),V(i+NumberOfVertices),V(i+2*NumberOfVertices),5); %Plots scatter points
    %     text(V(i),V(i+NumberOfVertices),V(i+2*NumberOfVertices),string(ID(i))); %Prints ID numbers next at vertex locations
    % end

    %% Plots the Fces constructing the geometry
    % trisurf(FV)

    %% Identifies the faces that are on the curved surface and stores their ID's to a vector "ID_CurvedSurfaces". This is done by looping over all normal vectors and obtaining the ones that have only one normal component equal to zero. Faces with two normal components equal to zero are flat surfaces and are not wanted. 
    ID_CurvedSurfaces = []; % Lists the ID's of the faces on the curved surfaces of the cylindrical shapes in the geometry
    Counter_CurvedSurfaces = 1; %Counts the number of curved faces on the surfaces

    for i = 1 : NumberOfNormalVectors
        if ((abs(N(i)) < tolerance) && (abs(N(i+NumberOfNormalVectors)) > tolerance) && (abs(N(i+2*NumberOfNormalVectors)) > tolerance))||((abs(N(i))  > tolerance) && (abs(N(i+NumberOfNormalVectors))  > tolerance) && (abs(N(i+2*NumberOfNormalVectors)) < tolerance)) || ((abs(N(i))  > tolerance) && (abs(N(i+NumberOfNormalVectors)) < tolerance) && (abs(N(i+2*NumberOfNormalVectors)) > tolerance))
            ID_CurvedSurfaces(Counter_CurvedSurfaces)= i;
            Counter_CurvedSurfaces = Counter_CurvedSurfaces+1;
        end

    end

    %%
    counter1 = 0 ;
    CylinderNr = 0; %counts the number fo cylinders in the gearbox
    LargestNodeID = 1;

    Radius = []; % Contains the radius values of all the cylinders
    MaxX = []; MinX = []; % minimum and maximum X values of the cylinders (boundary values)
    MaxY = []; MinY = []; % minimum and maximum Y values of the cylinders (boundary values)
    MaxZ = []; MinZ = []; % minimum and maximum Z values of the cylinders (boundary values)
    Dx =[]; Dy = []; Dz = []; %Diameter values in X, Y and,Z directions
    CylinderCenters_X_coordinates = []; % Contains the X-coordinate of the centers of all cylinders
    CylinderCenters_Y_coordinates = []; % Contains the Y-coordinate of the centers of all cylinders
    CylinderCenters_Z_coordinates = []; % Contains the Y-coordinate of the centers of all cylinders
    CylinderHeight = []; %contains the length of all cylinders ordered in their respective ID's

    while (counter1 < NumberOfVertices)
        counter1 = counter1+1;
        check = 1;

        CylinderPoints = [];

        CylinderPoints(end+1) =  F(ID_CurvedSurfaces(counter1));
        CylinderPoints(end+1) =  F(ID_CurvedSurfaces(counter1) + NumberOfNormalVectors);
        CylinderPoints(end+1) =  F(ID_CurvedSurfaces(counter1) + 2*NumberOfNormalVectors);

        counter1 = counter1 + 2;

        while (check ~=0)
           counter2 = counter1; 

            for i = 1:length(ID_CurvedSurfaces)
                   if ((ismember(F(ID_CurvedSurfaces(i)), CylinderPoints)) || (ismember(F(ID_CurvedSurfaces(i)+ NumberOfNormalVectors), CylinderPoints)) || (ismember(F(ID_CurvedSurfaces(i) + 2*NumberOfNormalVectors), CylinderPoints)))

                       if (ismember(F(ID_CurvedSurfaces(i)), CylinderPoints))
                           %do nothing;
                       else
                           CylinderPoints(end+1) =  F(ID_CurvedSurfaces(i));
                           counter1 = counter1+1;
                       end

                       if (ismember(F(ID_CurvedSurfaces(i)+ NumberOfNormalVectors), CylinderPoints))
                           %do nothing;

                       else
                           CylinderPoints(end+1) =  F(ID_CurvedSurfaces(i) + NumberOfNormalVectors);
                           counter1 = counter1+1;
                       end

                       if (ismember(F(ID_CurvedSurfaces(i) + 2*NumberOfNormalVectors), CylinderPoints))
                           %do nothing;
                       else
                           CylinderPoints(end+1) =  F(ID_CurvedSurfaces(i) + 2*NumberOfNormalVectors);
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
            CylindersY(end+1) = V(CylinderPoints(i)+ NumberOfVertices); 
            CylindersZ(end+1) = V(CylinderPoints(i)+ 2*NumberOfVertices);   
        end

        MaxX(end+1) = max(CylindersX);
        MinX(end+1) = min(CylindersX);
        MaxY(end+1) = max(CylindersY);
        MinY(end+1) = min(CylindersY);
        MaxZ(end+1) = max(CylindersZ);
        MinZ(end+1) = min(CylindersZ);

        Dx(end+1) = MaxX(CylinderNr) - MinX(CylinderNr);
        Dy(end+1) = MaxY(CylinderNr) - MinY(CylinderNr);
        Dz(end+1) = MaxZ(CylinderNr) - MinZ(CylinderNr);

        if (abs(Dx(CylinderNr)- Dy(CylinderNr)) < tolerance)
            Radius(end+1)= Dx(CylinderNr)/2;
            CylinderCenters_X_coordinates(CylinderNr) = MinX(CylinderNr) + Radius(CylinderNr);
            CylinderCenters_Y_coordinates(CylinderNr) = MinY(CylinderNr) + Radius(CylinderNr);
            CylinderHeight(CylinderNr) = abs(MaxZ(CylinderNr) -  MinZ(CylinderNr));

            zvector=linspace(MinZ(CylinderNr),MaxZ(CylinderNr));
            zvector=linspace(MinZ_Abs -2, MaxZ_Abs + 2);

            %line(CylinderCenters_X_coordinates(CylinderNr)*ones(1,length(zvector)),CylinderCenters_Y_coordinates(CylinderNr)*ones(1,length(zvector)),zvector);        

        elseif (abs(Dx(CylinderNr)- Dz(CylinderNr)) < tolerance)   
            Radius(end+1)= Dx(CylinderNr)/2;
            CylinderCenters_X_coordinates(CylinderNr) = MinX(CylinderNr) + Radius(CylinderNr);
            CylinderCenters_Z_coordinates(CylinderNr) = MinZ(CylinderNr) + Radius(CylinderNr);
            CylinderHeight(CylinderNr) = abs(MaxY(CylinderNr) -  MinY(CylinderNr));


            xvector=linspace(MinX(CylinderNr),MaxX(CylinderNr));
            yvector=linspace(MinY_Abs -2, MaxY_Abs +2);

            %line(CylinderCenters_X_coordinates(CylinderNr)*ones(1,length(yvector)),yvector,CylinderCenters_Z_coordinates(CylinderNr)*ones(1,length(yvector)));        

        elseif (abs(Dy(CylinderNr)- Dz(CylinderNr)) < tolerance)
            Radius(end+1)= Dy(CylinderNr)/2;
            CylinderCenters_Y_coordinates(CylinderNr) = MinY(CylinderNr) + Radius(CylinderNr);
            CylinderCenters_Z_coordinates(CylinderNr) = MinZ(CylinderNr) + Radius(CylinderNr);
            CylinderHeight(CylinderNr) = abs(MaxX(CylinderNr) -  MinX(CylinderNr));

            xvector=linspace(MinX(CylinderNr),MaxX(CylinderNr));
            xvector=linspace(MinX_Abs -2, MaxX_Abs +2);

            %line(xvector,CylinderCenters_Y_coordinates(CylinderNr)*ones(1,length(xvector)),CylinderCenters_Z_coordinates(CylinderNr)*ones(1,length(xvector)));
        end
    end

    %
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
            for Counter_CurvedSurfaces = 1:length(CylinderCenters_Y_coordinates) %Loop over all cylinders and use following IF-STATEMENT to find the cylinders that lie on the same axis as the StartingAxis
                if ((abs(CylinderCenters_Y_coordinates(Counter_CurvedSurfaces) - StartingAxisY) < tolerance)&&(Counter_CurvedSurfaces ~= MotorIndex) && (Counter_CurvedSurfaces~= LastCylinderID))
                    for k = 1:length(CylinderCenters_Y_coordinates) %Loop over all other cylinders lying on diferent axes and searching for meshing cylinders(gears)
                        if ((abs(CylinderCenters_Y_coordinates(k) - StartingAxisY)) > tolerance) && (k ~= MotorIndex) % && (j~=k)
                            DistanceBetweenAxes = sqrt(( StartingAxisY - CylinderCenters_Y_coordinates(k))^2 + (StartingAxisZ-CylinderCenters_Z_coordinates(k))^2);
                                if ((((Radius(Counter_CurvedSurfaces)+ Radius(k))- DistanceBetweenAxes ) > 0) && (abs(MinX(k)-MinX(Counter_CurvedSurfaces))<tolerance)&& (k~= LastCylinderID)) %&& (( DistanceBetweenAxes - (Radius(j)+ Radius(k))) > -0.01) 
%                                     fprintf('main cylinder is number %i and other cylinder is number %i \n',Counter_CurvedSurfaces,k)
                                    StartingAxisY = CylinderCenters_Y_coordinates(k);
                                    StartingAxisZ = CylinderCenters_Z_coordinates(k);
                                    LastCylinderID = k;
                                 end
                        end
                    end
                end 
            end
        end

        xvector = linspace(MinX_Abs - 50, MaxX_Abs + 50);
        %line(xvector,CylinderCenters_Y_coordinates(LastCylinderID)*ones(1,length(xvector)),CylinderCenters_Z_coordinates(LastCylinderID)*ones(1,length(xvector)),'color','red'); 
        OutputShaftAxisX = [CylinderCenters_Y_coordinates(LastCylinderID),CylinderCenters_Z_coordinates(LastCylinderID)];
        EM_start = [MaxX(MotorIndex),MotorCenterY, MotorCenterZ];
        EM_vector = [(MinX(MotorIndex)-MaxX(MotorIndex)),0,0];
        Anchor = [MaxX_Abs, OutputShaftAxisX(1), OutputShaftAxisX(2)];
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
            for Counter_CurvedSurfaces = 1:length(CylinderCenters_X_coordinates) %Loop over all cylinders and use following IF-STATEMENT to find the cylinders that lie on the same axis as the StartingAxis
                if ((abs(CylinderCenters_X_coordinates(Counter_CurvedSurfaces) - StartingAxisX) < tolerance)&&(Counter_CurvedSurfaces ~= MotorIndex) && (Counter_CurvedSurfaces~= LastCylinderID))
                    for k = 1:length(CylinderCenters_X_coordinates) %Loop over all other cylinders lying on diferent axes and searching for meshing cylinders(gears)
                        if ((abs(CylinderCenters_X_coordinates(k) - StartingAxisX)) > tolerance) && (k ~= MotorIndex)
                            DistanceBetweenAxes = sqrt(( StartingAxisX - CylinderCenters_X_coordinates(k))^2 + (StartingAxisZ-CylinderCenters_Z_coordinates(k))^2);
                                if ((((Radius(Counter_CurvedSurfaces)+ Radius(k))- DistanceBetweenAxes ) > 0) && (abs(MinY(k)-MinY(Counter_CurvedSurfaces))<tolerance)&& (k~= LastCylinderID))
%                                     fprintf('main cylinder is number %i and other cylinder is number %i \n',Counter_CurvedSurfaces,k)
                                    StartingAxisX = CylinderCenters_X_coordinates(k);
                                    StartingAxisZ = CylinderCenters_Z_coordinates(k);
                                    LastCylinderID = k;
                                 end
                        end
                    end
                end 
            end
        end
        yvector = linspace(MinY_Abs -50, MaxY_Abs + 50);
        %line(CylinderCenters_X_coordinates(LastCylinderID)*ones(1,length(yvector)),yvector,CylinderCenters_Z_coordinates(LastCylinderID)*ones(1,length(yvector)),'color','red');
        OutputShaftAxisY = [CylinderCenters_X_coordinates(LastCylinderID),CylinderCenters_Z_coordinates(LastCylinderID)];
        EM_start = [MotorCenterX, MaxY(MotorIndex), MotorCenterZ];
        EM_vector = [0,(MinY(MotorIndex)-MaxY(MotorIndex)),0];
        Anchor = [OutputShaftAxisY(1), MaxY_Abs,OutputShaftAxisY(2) ];

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
        for Counter_CurvedSurfaces = 1:length(CylinderCenters_X_coordinates) %Loop over all cylinders and use following IF-STATEMENT to find the cylinders that lie on the same axis as the StartingAxis
            if ((abs(CylinderCenters_X_coordinates(Counter_CurvedSurfaces) - StartingAxisX) < tolerance)&&(Counter_CurvedSurfaces ~= MotorIndex) && (Counter_CurvedSurfaces~= LastCylinderID))
                for k = 1:length(CylinderCenters_X_coordinates) %Loop over all other cylinders lying on diferent axes and searching for meshing cylinders(gears)
                    if ((abs(CylinderCenters_X_coordinates(k) - StartingAxisX)) > tolerance) && (k ~= MotorIndex) 
                        DistanceBetweenAxes = sqrt(( StartingAxisX - CylinderCenters_X_coordinates(k))^2 + (StartingAxisY-CylinderCenters_Y_coordinates(k))^2);
                            if ((((Radius(Counter_CurvedSurfaces)+ Radius(k))- DistanceBetweenAxes ) > 0) && (abs(MinZ(k)-MinZ(Counter_CurvedSurfaces))<tolerance)&& (k~= LastCylinderID)) 
%                                 fprintf('main cylinder is number %i and other cylinder is number %i \n',Counter_CurvedSurfaces,k)
                                StartingAxisX = CylinderCenters_X_coordinates(k);
                                StartingAxisY = CylinderCenters_Y_coordinates(k);
                                LastCylinderID = k;
                             end
                    end
                end
            end 
        end
    end
        zvector = linspace(MinZ_Abs -50, MaxZ_Abs +50);
        %line(CylinderCenters_X_coordinates(LastCylinderID)*ones(1,length(zvector)),CylinderCenters_Y_coordinates(LastCylinderID)*ones(1,length(zvector)),zvector,'color','red');
        OutputShaftAxisZ = [CylinderCenters_X_coordinates(LastCylinderID),CylinderCenters_Y_coordinates(LastCylinderID)];
        EM_start = [MotorCenterX,MotorCenterY,MinZ(MotorIndex) ];
        EM_vector = [0,0,(MaxZ(MotorIndex)-MinZ(MotorIndex))];
        Anchor = [OutputShaftAxisZ(1), OutputShaftAxisZ(2), MinZ_Abs];
    end
end





