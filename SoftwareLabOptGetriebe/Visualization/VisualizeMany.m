function VisualizeMany(n, iter,bauraum, gearbox,Anchor, EC1,EC1_width, EC1_center, EC1_SurfCenter, EC2, EC2_center, EC2_SurfCenter, EC3, EC3_center, EC3_SurfCenter, EC4, EC4_center, EC4_SurfCenter, EM_start,EM_vector, Fitness,r,y,theta,Rotations,EC2_placement)
%Parent function of :func:`Visualize` that uses the same inputs except for the
%counter **i**. This counter is created in VisualizeMany so Visualize
%can be called multiple times and print many solutions. Can be used to
%print the desired amount of best solutions for a defined generation.
%
%
%
% .. code-block:: 
%
%   VisualizeMany(n, iter,bauraum, gearbox,Anchor, EC1,EC1_width, EC1_center, EC1_SurfCenter, EC2, EC2_center, EC2_SurfCenter, EC3, EC3_center, EC3_SurfCenter, EC4, EC4_center, EC4_SurfCenter, EM_start,EM_vector, Fitness,r,y,theta,Rotations,EC2_placement)
%

for i = 1:n
    
    Visualize(i,n, iter,bauraum, gearbox,Anchor, EC1,EC1_width, EC1_center, EC1_SurfCenter, EC2, EC2_center, EC2_SurfCenter, EC3, EC3_center, EC3_SurfCenter, EC4, EC4_center, EC4_SurfCenter, EM_start,EM_vector, Fitness,r,y,theta,Rotations,EC2_placement)
    
end

end