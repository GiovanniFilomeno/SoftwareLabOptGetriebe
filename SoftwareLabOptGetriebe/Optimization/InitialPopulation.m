function [r, y, theta, rotations, EC2_placement ]=InitialPopulation(rad_range,pop,EM_start,EM_vector,EM_radius,EC1_height,EC1_width,EC2_height,EC2_width,EC3_height,EC3_width)
% Function that reads the possible rotation of the gearbox, as well as the
% dimensions of the Electric motor and Electric components and then creates
% random sets of values of coordinates and angles between these limits. The
% number of sets created is given by the size of population **pop**. The
% generation of each Electric Component coordinates and angles is connected
% to the values of the other neighboring geometries, as this code application
% assumes the best solution as the one where the ECs are wrapped around the
% Electric motor.
%
% :param double rad_range: Range of rotation of the gearbox without going
%                          out of the bauraum. Obtained from :func:`PositionGearbox`.
% :param int pop: Size of the population for the genetic algorithm. User
%                 defined value.
% :param array EM_start: Vector of the Electric motor center coordinates.
% :param array EM_vector: Vector of Gearbox and Electric motor direction. 
% :param array EM_radius: Radius of the electric motor.
% :param double EC1_height: Height of the IGBT.
% :param double EC1_width: Width of the IGBT.
% :param double EC2_height: Height of the Capacitor.
% :param double EC2_width: Width of the Capacitor
% :param double EC3_height: Height of the EMI-filter
% :param double EC3_width: Width of the EMI- filter
%
% :return: 
%   *[r, y, theta, rotations, EC2_placement]*
%       - **r** : popx4(1 column per EC) array with random values that define the perpendicular distance
%         between the center of the electric motor and the respective ECs.
%       - **y** : popx4 array with random values that define the y-axis
%         location of the ECs. These random values are limited between the
%         minimum and maximum Y coordinate of the Electric motor.
%       - **theta** : popx4 array with random values between 0 and 2*pi that
%         define the position of the ECs around the circular Electric Motor.
%       - **rotations** : popx5(4 ECs and gearbox) array with random values
%         that define the rotations of the gearbox along its drive axis, as
%         well as the rotation of the ECs along its own center.
%       - **EC2_placement** : popx1 array with random boolean values.
%         Whether the value is 1 or 0 defines if the Capacitor will be
%         positioned to the left or right of the IGBT.
%
% :rtype: [ double array, double array, double array, double array, double array] 
%
% **Example in Code**
%
% .. code-block::
%
%   [r, y, theta, rotations, EC2_placement ]=InitialPopulation(rad_range,pop,EM_start,EM_vector,EM_radius,EC1_height,EC1_width,EC2_height,EC2_width,EC3_height,EC3_width);
%

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