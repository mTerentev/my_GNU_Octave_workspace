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

global y;
y = @(u,v) Rot(u)*(rack(v)+[R; -R*u; 0]);

F = @(u, v) cross(df_dn(y, 1, u, v), df_dn(y, 2, u, v))(3);

function sol = Newton_Raphson(f, init, iters = 10)
  df_dx = @(f, x, d=10^(-6)) ...
      ( f(x + d/2) - f(x - d/2) ) / d;
  sol = init;
  for i=1:iters
    sol = sol - f(sol)/df_dx(f, sol);
    if f(sol) < 10^(-5)
      break
    endif
  endfor
  if f(sol) >= 2*10^(-5)
    sol = Newton_Raphson(f, init + 1);
  endif
end


global hu hv;
uslider = uicontrol (                    ...
         'style', 'slider',                ...
         'Units', 'normalized',            ...
         'position', [0, 0.05, 1, 0.05], ...
         'min', -pi/n,                         ...
         'max', pi/n,                        ...
         'value', 0,                      ...
         'callback', {@getU}          ...
       );
function getU (h, event) global hu; hu = get(h, 'value'); update(); end

vslider = uicontrol (                    ...
         'style', 'slider',                ...
         'Units', 'normalized',            ...
         'position', [0, 0, 1, 0.05], ...
         'min', -pi/n,                         ...
         'max', pi/n,                        ...
         'value', 0,                      ...
         'callback', {@getV}          ...
       );
function getV (h, event) global hv; hv = get(h, 'value'); update(); end



function update()
  global plot1 plot2 plot4 plot5 y sv hu hv df_dn;
  hu;
  hv;
  ph = y(hu, hv);
  set(plot1, "xdata",ph(1,:));
  set(plot1, "ydata",ph(2,:));
  set(plot1, "zdata",ph(3,:));

  rh = cat(2,arrayfun(@(v) y(hu,v), sv, "UniformOutput", false){:,:});
  set(plot2, "xdata",rh(1,:));
  set(plot2, "ydata",rh(2,:));
  set(plot2, "zdata",rh(3,:));

  du = ph + normalize(df_dn(y, 1, hu, hv))/5;
  dv = ph + normalize(df_dn(y, 2, hu, hv))/5;

  set(plot4, "xdata",[ph(1),du(1)]);
  set(plot4, "ydata",[ph(2),du(2)]);
  set(plot4, "zdata",[ph(3),du(3)]);

  set(plot5, "xdata",[ph(1),dv(1)]);
  set(plot5, "ydata",[ph(2),dv(2)]);
  set(plot5, "zdata",[ph(3),dv(3)]);

  refresh();
end

u_res = 1000;
v_res = 100;

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

V = zeros(size(su));

for i = 1:size(su,2)
  V(:,i) = Newton_Raphson(@(v) F(su(i),v), -5);
end

r1 = cat(2,arrayfun(y, su, V, "UniformOutput", false){:,:});
plot3(r1(1,:),r1(2,:),r1(3,:), "r", "LineWidth", 1, "marker", "o");


rh = cat(2,arrayfun(@(v) y(0,v), sv, "UniformOutput", false){:,:});
global plot1 plot2 plot4 plot5
plot2 = plot3(rh(1,:),rh(2,:),rh(3,:), "r", "LineWidth", 3);

ph = y(0, 0);
plot1 = plot3(ph(1,:),ph(2,:),ph(3,:), "o", "LineWidth", 5);

du = ph + df_dn(y, 1, 0, 0);
dv = ph + df_dn(y, 2, 0, 0);

plot4 = plot3([ph(1),du(1)],[ph(2),du(2)],[ph(3),du(3)], "b", "LineWidth", 3);
plot5 = plot3([ph(1),dv(1)],[ph(2),dv(2)],[ph(3),dv(3)], "g", "LineWidth", 3);

hold off;



