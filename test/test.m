R = 1;
ro = 1.2;
ri = 0.8;
alp = 20*pi/180;
n = 5;
rc = 0.1;

% Then define functions below...
function y = ftooth(t,R,ro,ri,alp,n,rc)

  tup=(ro-R)*n/(2*pi*R);
  tdown=(ri-R)*n/(2*pi*R);
  x1 = 0.25-(tup-rc)*tan(alp)-rc/cos(alp);
  x2 = x1+rc*cos(alp);
  x4 = 0.25-(tdown+rc)*tan(alp)+rc/cos(alp);
  x3 = x4-rc*cos(alp);

  y = (t<x1).*tup + ...
      (x1<=t & t<x2).*cf(x1,tup-rc,t,1,rc) + ...
      (x2<=t & t<x3).*tf(t,alp) + ...
      (x3<=t & t<x4).*cf(x4,tdown+rc,t,-1,rc) + ...
      (t>=x4).*tdown;
end

function R = Rot(alp)
  R = [
    cos(alp),-sin(alp),     0;
    sin(alp), cos(alp),     0;
           0,        0,     1;
  ]
end

rack = @(t) [t; ftooth(abs(t),R,ro,ri,alp,n,rc); 0]

y = @(u,v) Rot(u)*rack(v)

%u = linspace(0, 0.5, 10)
%v = linspace(0, 0.5, 10),


u_res = 100
v_res = 10
[u,v] = meshgrid(linspace(-0.5, 0.5, u_res), linspace(-0.5, 0.5, v_res));

r = reshape(cat(2,arrayfun(y, u, v, "UniformOutput", false){:,:}),3,u_res,v_res)
hold on
for i = 1:v_res
  plot3(r(1,:,i),r(2,:,i),r(3,:,i))
end
hold off
