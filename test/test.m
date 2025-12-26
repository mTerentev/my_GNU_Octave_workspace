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
      (x4<=t).*tdown;

s = @(x) size(x,2);

Rot = @(alp) reshape([
    reshape(cos(alp),1,1,s(alp)) ,reshape(-sin(alp),1,1,s(alp));
    reshape(sin(alp),1,1,s(alp)), reshape(cos(alp),1,1,s(alp));
  ],2,2,s(alp));

rack = @(t) [reshape(ftooth(abs(t)),1,s(t)); reshape(t,1,s(t))];

vrshp = @(v) [reshape(v(1),1,s(v)); reshape(v(2),1,s(v))];

h = 0.01;

U = -pi/n:h:pi/n;
V = -pi/n:h:pi/n;


Y = Rot(U);

[u,v] = meshgrid(1:size(U,2), 1:size(V,2));

X = rack(V);
X1 = [reshape(R*ones(s(U),s(V)),1,s(U)*s(V)); reshape(-R*U(u),1,s(U)*s(V))];

M=Y(:,:,u);

M1=X(:,v);

R = zeros(2, 1, size(M,3));



for k = 1:size(M,3)
    R(:,k) = M(:,:,k) * (M1(:,k)+X1(:,k));
end

R = reshape(R, 2, size(u), size(v));

%r=num2cell(Rot(u),[1,2])
%r=y(u,v)
hold on
for k = 1:size(R,3)
  plot(R(1,:,k), R(2,:,k), "k");
end
hold off
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
