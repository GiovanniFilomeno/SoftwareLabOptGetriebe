function [r, y, theta, rotations, EC2_placement ]=InitialPopulation(rad_range,pop,EM_start,EM_vector,EM_radius,EC1_height,EC1_width,EC2_height,EC2_width,EC3_height,EC3_width)

gb_rot = rad_range.*rand(pop,1);


%IGBT
r1 = EM_radius + (EC1_height/2) + (EC1_height).*rand(pop,1);
y1= EM_start(2) + (EM_vector(2)).*rand(pop,1);
theta1 = -2*pi + 2*pi.*rand(pop,1);
EC1_rot = -pi/12 + (pi/6).*rand(pop,1);



%Condensator
r2 = EM_radius + (EC2_height/2) + (EC2_height).*rand(pop,1);
y2= EM_start(2) + (EM_vector(2)).*rand(pop,1);
angle_range2 = (EC1_width + EC2_width)/(2*(EM_radius)); 
theta2_direction = rand(pop,1);
EC2_placement = round(theta2_direction);
theta2 = zeros(pop,1);
for i=1:pop
    if EC2_placement(i) == 1
        theta2(i) = theta1(i) + angle_range2 + (pi/15)*rand(); %Component counterclockwise from 1
    else
        theta2(i) = theta1(i) - angle_range2 - (pi/15)*rand(); %Component clockwise from 1
    end
end   
EC2_rot = -pi/12 + (pi/6).*rand(pop,1);


%EMI-Filter
r3 = EM_radius + (EC3_height/2) + (EC3_height).*rand(pop,1);
y3= EM_start(2) + (EM_vector(2)).*rand(pop,1);
angle_range3 = (EC2_width + EC3_width)/(2*(EM_radius)); 
theta3 = zeros(pop,1);
for i=1:pop
    if theta2(i) > theta1(i)
        theta3(i) = theta2(i) + angle_range3 + (pi/15)*rand();
    else
        theta3(i) = theta2(i) - angle_range3 - (pi/15)*rand();
    end
end   
EC3_rot = -pi/12 + (pi/6).*rand(pop,1);


%Control board
r4 = r1 + (EC1_height) + (EC1_height).*rand(pop,1);
y4= EM_start(2) + (EM_vector(2)).*rand(pop,1);
theta4 = theta1;
EC4_rot = EC1_rot;

r = [r1,r2,r3,r4];

y = [y1,y2,y3,y4];

theta = [theta1,theta2,theta3,theta4];

rotations = [gb_rot, EC1_rot, EC2_rot, EC3_rot, EC4_rot];


 end