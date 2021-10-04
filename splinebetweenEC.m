clear
clc

x= [-1 0 2 3];
y= [2 2 4 4];
z= [9 9 -3 -3];
% cs = spline(x,[0 y 0], [0 z 0]);
% xx = linspace(-1,3,51);
% plot(x,y,z,'o',xx,ppval(cs,xx),'-');

plot3(x,y, z,'ro','LineWidth',2);
text(x,y,z,[repmat('  ',4,1), num2str((1:4)')])
ax = gca;
ax.XTick = [];
ax.YTick = [];
ax.ZTick = [];
xyz = [x ;y ;z];
hold on
fnplt(cscvn(xyz),'r',2)
hold off