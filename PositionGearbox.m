function [gearbox, Anchor, EM_start, EM_vector, EM_end, rad_range] = PositionGearbox(bauraum,gearbox,EM_start, EM_vector, Anchor, point, zero)

%setting rotation angle values
if (EM_vector(1) ~= 0)
    if (EM_vector(1) > 0)
        anglex = 0;
        angley = 0;
        anglez = pi/2;
    else
        anglex = 0;
        angley = 0;
        anglez = -pi/2;
    end
end
    
if (EM_vector(2) ~= 0)
    if (EM_vector(2) > 0)
        anglex = 0;
        angley = 0;
        anglez = 0;
    else
        anglex = 0;
        angley = 0;
        anglez = 0;
    end
end
        
if (EM_vector(3) ~= 0)
    if (EM_vector(3) > 0)
        anglex = 0;
        angley = pi/2;
        anglez = -pi/2;
    else
        anglex = 0;
        angley = -pi/2;
        anglez = pi/2;
    end
end

%Positioning of the gearbox inside he Bauraum
gearbox.vertices = gearbox.vertices - Anchor + point;
EM_start = EM_start - Anchor + point;
Anchor = point;
gearbox.vertices = Rotation(gearbox.vertices,Anchor, 1, anglex);
gearbox.vertices = Rotation(gearbox.vertices,Anchor, 2, angley);
gearbox.vertices = Rotation(gearbox.vertices,Anchor, 3, anglez);
EM_start = Rotation(EM_start, Anchor, 1, anglex);
EM_start = Rotation(EM_start, Anchor, 2, angley);
EM_start = Rotation(EM_start, Anchor, 3, anglez);
EM_vector = Rotation(EM_vector, zero, 1, anglex);
EM_vector = Rotation(EM_vector, zero, 2, angley);
EM_vector = Rotation(EM_vector, zero, 3, anglez);
EM_end = EM_start + EM_vector;

out_gearbox = [];
for i= 1:360
    
    gearbox.vertices = Rotation(gearbox.vertices,Anchor, 2, pi/180);
    EM_start = Rotation(EM_start, Anchor, 2, pi/180);
    EM_vector = Rotation(EM_vector, zero, 2, pi/180);
    EM_end = EM_start + EM_vector;
    EM_shaftx = [EM_start(1,1) EM_start(1,1) + EM_vector(1,1)];
    EM_shafty = [EM_start(1,2) EM_start(1,2) + EM_vector(1,2)];
    EM_shaftz = [EM_start(1,3) EM_start(1,3) + EM_vector(1,3)];
    
    IN = inpolyhedron(bauraum,gearbox.vertices);
    out_gearbox(i) = sum(IN)/length(IN); 
end

valid_angles = [];
max_IN = max(out_gearbox)


for i= 1:360
    if (out_gearbox(i) == max_IN)
        valid_angles(end+1) = i;
        valid_radians = valid_angles.*(pi/180);
    end
end
out_gearbox = [];
angle_range = valid_angles(length(valid_angles)) - valid_angles(1);
rad_range = (length(valid_angles))*pi/180;

if angle_range == length(valid_angles) - 1
    first_angle = valid_radians(1);
else
    for i=1:length(valid_radians)
        if i==1
            valid_diff(i) = valid_radians(1) - valid_radians(length(valid_radians)); 
        else
            valid_diff(i) = valid_radians(i) - valid_radians(i-1);
        end
    end
    [srt,I] = sort(valid_diff, 'descend');
    first_angle = valid_radians(I(1,1));
end

gearbox.vertices = Rotation(gearbox.vertices,Anchor, 2, first_angle);
EM_start = Rotation(EM_start, Anchor, 2, first_angle);
EM_end = EM_start + EM_vector;




end