pop = 100;
iter = 10;
ulim = 5;
%Create first sample population and evaluate fitness
x = (ulim).*rand(pop,3);
angle = (2*pi).*rand(pop,1);
X = zeros(pop,2);
Y = zeros(pop,2);
Z = zeros(pop,2);
for i=1:pop
X(i,1) = x(i,1);
X(i,2) = x(i,1)+sin(angle(i,1));
Y(i,1) = x(i,2);
Y(i,2) = x(i,2)+cos(angle(i,1));
Z(i,1) = x(i,3);
Z(i,2) = x(i,3);
end
drive = [1,0,0];
Xdrive = [0 5];
Ydrive = [0 0];
Zdrive = [0 0];
d = zeros(pop,1);
theta = zeros(pop,1);
for i=1:pop
    d(i,1) = norm(x(i,:));
    vector = [X(i,2)-X(i,1) Y(i,2)-Y(i,1) Z(i,2)-Z(i,1)];
    theta(i,1) = atan2(norm(cross(drive,vector)), dot(drive,vector));
end

%Plot each generation
    figure(1)
    plot3(Xdrive',Ydrive',Zdrive', 'LineWidth', 3)
    hold on
    plot3(X',Y',Z')
    axis([0 ulim 0 ulim 0 ulim])
    hold on
    pause(1)

for g = 2:iter
    for k = 1: 3: pop
        %Selection
        sel_coeff = mean(d);
        sel_coeff2 = mean(theta);
        parents = zeros(3,4);
        for j=1:3
            rand_sel = 1 + round(pop.*rand());
            for i = rand_sel:pop
                if d(i,1)<sel_coeff && theta(i,1)<sel_coeff2
                    parent(j,1) = x(i,1);
                    parent(j,2) = x(i,2);
                    parent(j,3) = x(i,3);
                    parent(j,4) = angle(i);
                    break
                end
            end
        end

        %Crossover
        child = [parent(1,1) parent(2,2) parent(3,3);parent(1,2) parent(2,3) parent(3,1);parent(1,3) parent(2,1) parent(3,2)];
        child_angle = [parent(1,4) parent(2,4) parent(3,4)];
        %Mutation
        for i = 1:3
            for j = 1:3
                mut_prob = rand();
                if mut_prob<=0.05
                    child(i,j) = child(i,j) - (sel_coeff*0.50);                    
                    if child(i,j)<0
                        child(i,j) = abs(child(i,j));
                    end
                end
            end
            if mut_prob<0.5
                child_angle(i) = child_angle(i) - (sel_coeff2*0.20);
                if child_angle(i)<0
                        child_angle(i) = abs(child_angle(i));
                end
            end
        end
        
        %Write new generation 
        x(k,1) = child(1,1);
        x(k,2) = child(1,2);
        x(k,3) = child(1,3);
        x(k+1,1) = child(2,1);
        x(k+1,2) = child(2,2);
        x(k+1,3) = child(2,3);
        x(k+2,1) = child(2,1);
        x(k+2,2) = child(2,2);
        x(k+2,3) = child(2,3);
        angle(k,1) =  child_angle(1);
        angle(k+1,1) = child_angle(2);
        angle(k+2,1) = child_angle(3);
    end
    for i=1:pop
        X(i,1) = x(i,1);
        X(i,2) = x(i,1)+cos(angle(i,1));
        Y(i,1) = x(i,2);
        Y(i,2) = x(i,2)+sin(angle(i,1));
        Z(i,1) = x(i,3);
        Z(i,2) = x(i,3);
        end
    for i=1:pop
        d(i,1) = norm(x(i,:));
        vector = [X(i,2)-X(i,1) Y(i,2)-Y(i,1) Z(i,2)-Z(i,1)];
        theta(i,1) = atan2(norm(cross(drive,vector)), dot(drive,vector));
    end
    [srt,I]=sort(d);
    best = x(I(1),:);
    %Plot each generation
    figure(g)
    plot3(Xdrive',Ydrive',Zdrive', 'LineWidth', 2)
    hold on
    plot3(X',Y',Z')
    axis([0 ulim 0 ulim 0 ulim])
    hold on
    pause(1)
end
        


