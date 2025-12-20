
R = 1;
ro = 1.2;
ri = 0.8;
alp = 20*pi/180;
n = 5;
rc = 0.01;

tup=(ro-R)*n/(2*pi*R);
tdown=(ri-R)*n/(2*pi*R);
x1 = 0.25-(tup-rc)*tan(alp)-rc/cos(alp);
x2 = x1+rc*cos(alp);
x4 = 0.25-(tdown+rc)*tan(alp)+rc/cos(alp);
x3 = x4-rc*cos(alp);

ftooth = @(t) piecewise(
      t<x1, tup,
      t<x2, cf(x1,tup-rc,t,1,rc),
      t<x3, tf(t,alp),
      t<x4, cf(x4,tdown+rc,t,-1,rc),
      t>=x4, tdown
)

%pretty(sym(ftooth))

Rot = @(alp) [
    cos(alp),-sin(alp),     0;
    sin(alp), cos(alp),     0;
       0*alp,        0,     1;
  ]

%pretty(sym(ftooth))

rack = @(t) [ftooth(abs(t),R,ro,ri,alp,n,rc); t; 0]

y = @(u,v) Rot(u)*(rack(v)+[R; -R*u; 0])

%dy_du = @(u, v) diff(y, u)
%dy_dv = @(u, v) diff(y, v)


u_res = 10
v_res = 100

su = linspace(-pi/n, pi/n, u_res)
sv = linspace(-pi/n, pi/n, v_res),

axis equal
hold on
for i = su
  f = matlabFunction(sym(y))
  r = cat(2,arrayfun(@(v) f(i,v), sv, "UniformOutput", false){:,:})
  plot3(r(1,:),r(2,:),r(3,:), "k")
end
hold off
