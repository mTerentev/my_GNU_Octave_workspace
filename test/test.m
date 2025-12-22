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

tf = @(x) (0.25-x)/tan(alp);
cf = @(x,y,t,s) s*sqrt(rc^2-(t-x).^2) + y;

ftooth = @(t) (t<x1).*tup + ...
      (x1<=t & t<x2).*cf(x1,tup-rc,t,1) + ...
      (x2<=t & t<x3).*tf(t) + ...
      (x3<=t & t<x4).*cf(x4,tdown+rc,t,-1) + ...
      (x4<=t).*tdown

s = @(x) size(x,2)

Rot = @(alp) reshape([
    reshape(cos(alp),1,1,s(alp),s(alp)) ,reshape(-sin(alp),1,1,s(alp),s(alp)), reshape(alp*0,1,1,s(alp),s(alp));
    reshape(sin(alp),1,1,s(alp),s(alp)), reshape(cos(alp),1,1,s(alp),s(alp)), reshape(alp*0,1,1,s(alp),s(alp));
       reshape(alp*0,1,1,s(alp),s(alp)),    reshape(alp*0,1,1,s(alp),s(alp)), reshape(ones(size(alp*0)),1,1,s(alp),s(alp));
  ],3,3,s(alp),s(alp))

rack = @(t) [ftooth(abs(t)); t; 0*t]

y = @(u,v) num2cell(reshape(rack(v).+[R*ones(size(u)); -R*u; 0*u],3*s(u),s(v)),[1])

[u,v] = meshgrid(-0.5:0.1:0.5, -0.5:0.1:0.5)

M=Rot(u)

M1=reshape(rack(v).+[R*ones(size(u)); -R*u; 0*u],3,1,s(u),s(v))
% Reshape and use arrayfun
M_reshaped = reshape(M, 3, 3, []);
V_reshaped = reshape(M1, 3, 1, []);
R = zeros(3, 1, size(M_reshaped,3));

for k = 1:size(M_reshaped,3)
    R(:,:,k) = M_reshaped(:,:,k) * V_reshaped(:,:,k);
end

R = reshape(R, 3, size(u), size(v));

%r=num2cell(Rot(u),[1,2])
%r=y(u,v)
plot3(R(1,:,5), R(2,:,5), R(3,:,5))

%{
%disp(sym(y))

dy_du = @(u, v) diff(y, u)
dy_dv = @(u, v) diff(y, v)

F = @(u,v) (cross(dy_du(u,v), dy_dv(u,v)))

syms u v
disp(solve(F(u,v), u))

disp(F(1,2))
%}
function plot_family(f)

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
end

axis equal
