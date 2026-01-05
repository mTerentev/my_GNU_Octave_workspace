R = 1;
ro = 1.2;
ri = 0.8;
alp = 20*pi/180;
n = 5;
rc = 0.01;

tf = @(x) (0.25-x)/tan(alp);

cf = @(x,y,t,s) s*sqrt(rc^2-(t-x).^2) + y;

tup=(ro-R)*n/(2*pi*R);
tdown=(ri-R)*n/(2*pi*R);
x1 = 0.25-(tup-rc)*tan(alp)-rc/cos(alp);
x2 = x1+rc*cos(alp);
x4 = 0.25-(tdown+rc)*tan(alp)+rc/cos(alp);
x3 = x4-rc*cos(alp);

ftooth = @(t) ...
      (t<x1) * tup + ...
      (x1<=t & t<x2) * cf(x1,tup-rc,t,1) + ...
      (x2<=t & t<x3) * tf(t) + ...
      (x3<=t & t<x4) * cf(x4,tdown+rc,t,-1) + ...
      (x4<=t) * tdown;

Rot = @(alp) [
    [cos(alp),-sin(alp), 0];
    [sin(alp), cos(alp), 0];
    [       0,        0, 1];
  ];



global df_dn
df_dn = @(f, n, u, v, d=10^(-6)) ...
      ( f(u + d*(n==1)/2,v + d*(n==2)/2) - f(u - d*(n==1)/2,v - d*(n==2)/2) ) / d;

rack = @(t) [ftooth(abs(t)); t; 0];

global y F;
y = @(u,v) Rot(u)*(rack(v)+[R; -R*u; 0]);
y1 = @(u,v) Rot(u)*(rack(v)+[R+0.01; -R*u; 0]);

F = @(u, v) cross(df_dn(y, 1, u, v), df_dn(y, 2, u, v))(3);
F1 = @(u, v) cross(df_dn(y1, 1, u, v), df_dn(y1, 2, u, v))(3);

function sols = Newton_Raphson(f, inits, iters = 10)

  global df_dn;
  n = 2;
  m = 2;
  sols = inits;
  for i=1:iters
    J = zeros(n,m,size(inits,2),size(inits,3));
    for k=1:m
      for i=1:size(inits(2))
        for j=1:size(inits(3))
          J(:,k,i,j) = df_dn(f,k,sols(1,i,j), sols(2,i,j));
        endfor
      endfor
    endfor
    del = zeros(2,size(inits,2),size(inits,3));
    for i=1:size(inits(2))
        for j=1:size(inits(3))
          del(:,i,j) = J(:,:,i,j)^(-1)*f(sols(1,i,j),sols(2,i,j));
        endfor
    endfor
    sols = sols - del;
  endfor
end

function sols = NR_all_sols(f, u, v)
  m = size(u,2);
  n = size(v,2);

  [U,V] = meshgrid(u,v);
  inits = [reshape(U,1,m,n);reshape(V,1,m,n)];
  sols = zeros(2,5);
  l=0;

  sols = Newton_Raphson(f, inits);
  sols = cluster(sols)

endfunction


u_res = 1000;
v_res = 1000;

global su sv;
su = linspace(-pi/n, pi/n, u_res);
sv = linspace(-pi/n, pi/n, v_res);

hu=0;
hv=0;

axis equal;
hold on;
warning("off","all");

for u = su(1):0.1:su(end)
  f = @(v) y(u, v);
  r = cat(2,arrayfun(f, sv, "UniformOutput", false){:,:});
  plot3(r(1,:),r(2,:),r(3,:), "k");
end

U = zeros(size(sv));
Ux = zeros(size(sv));
U1 = zeros(size(sv));

points = zeros(3,size(sv,2));

for i = 1:size(sv,2)
  solu = fsolve(@(u) F(u,sv(i)), u);
  points(:,i) = y(solu,sv(i));

endfor
plot3(points(1,:),points(2,:),points(3,:), "r", "LineWidth", 3, "marker", "none");

filtered_points = zeros(3,size(sv,2));

l=1;
for i = 1:size(sv,2)

  sols = NR_all_sols(@(u,v) (y(u,v)-points(i))(1:2), -1:0.05:1, -1:0.05:1);

  if size(sols,2) == 1
    filtered_points(l) = point(i);
  endif
end

plot3(filtered_points(1,:),filtered_points(2,:),filtered_points(3,:), "g", "LineWidth", 3, "marker", "none");

%r1 = cat(2,arrayfun(y1, U1, sv, "UniformOutput", false){:,:});
%plot6 = plot3(r1(1,:),r1(2,:),r1(3,:), "b", "LineWidth", 1, "marker", "none");

%{
rh = cat(2,arrayfun(@(v) y(0,v), sv, "UniformOutput", false){:,:});
global plot1 plot2 plot4 plot5 plot6
plot2 = plot3(rh(1,:),rh(2,:),rh(3,:), "r", "LineWidth", 3);

ph = y(0, 0);
plot1 = plot3(ph(1,:),ph(2,:),ph(3,:), "o", "LineWidth", 5);

du = ph + df_dn(y, 1, 0, 0);
dv = ph + df_dn(y, 2, 0, 0);

plot4 = plot3([ph(1),du(1)],[ph(2),du(2)],[ph(3),du(3)], "b", "LineWidth", 3);
plot5 = plot3([ph(1),dv(1)],[ph(2),dv(2)],[ph(3),dv(3)], "g", "LineWidth", 3);
%}
hold off;



