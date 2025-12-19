t = linspace(0, 10, 100);
r = [t; t.^2; sin(t)];  % 3D vector function
plot3(r(1,:), r(2,:), r(3,:));
