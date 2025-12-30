%axis equal;
hold on

dist = @(v) sqrt(sum(v'*v));

r2 = zeros(size(r));

r3 = zeros(size(r));
r4 = zeros(size(r));

l = 1;
d = 50;
bias = -[ri;0;0];

r(:,500)+bias

for i=1:size(r,2)
  m = 1000;
  for j=i:1:min(i+d,size(r,2))
    dst = dist(r1(:,j)+bias)
    if dst < m
      m = dst;
    endif
  endfor
  m;
  if dist(r(:,i)+bias) < m
    r2(:,l) = r(:,i);
    l++;
  endif

  r3(:,i) = dist(r(:,i)+bias);
  r4(:,i) = dist(r1(:,i)+bias);

endfor

#r2 = resize(r2, 3, l-1);

plot3(r(1,:),r(2,:),r(3,:), "r", "LineWidth", 1, "marker", "none");
plot3(r1(1,:),r1(2,:),r1(3,:), "b", "LineWidth", 1, "marker", "none");

plot3(r2(1,:),r2(2,:),r2(3,:), "g", "LineWidth", 3, "marker", "o");


hold off
%{
hold on
plot(1:1000, r3);
plot(1:1000, r4);
hold off
%}
