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

function sol = Newton_Raphson(f, init=[rand()-0.5;rand()-0.5], iters = 10)

  global df_dn;
  n = 2;
  m = 2;
  sol = init;
  for i=1:iters
    J = zeros(n,m);
    for k=1:m
      J(:,k) = df_dn(f,k,sol(1), sol(2));
    endfor
    sol = sol - J^(-1)*f(sol(1),sol(2));
    if abs(sol) > 2
      sol = Newton_Raphson(f);
      return
    endif
    if abs(f(sol(1),sol(2))) < 10^(-5)
      return
    endif
  endfor
end


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

l=1
for i = 1:size(sv,2)
  solu = fsolve(@(u) F(u,sv(i)), u);
  point = y(solu,sv(i));

  sol2 = Newton_Raphson(@(u,v) (y(u,v)-point)(1:2));
  %Q = y(su(1),Qs);
  if abs([solu;sv(i)]-sol2)<0.01
    U(:,l) = solu;
    Ux(:,l) = sv(i);
    l++;
  endif
end
U=resize(U,1,l-1);
Ux=resize(Ux,1,l-1);
r = cat(2,arrayfun(y, U, Ux, "UniformOutput", false){:,:});
plot6 = plot3(r(1,:),r(2,:),r(3,:), "r", "LineWidth", 3, "marker", "none");

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



