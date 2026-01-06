R = 1;
ro = 1.2;
ri = 0.8;
alp = 20*pi/180;
n = 5;
rc = 0.05;

tf = @(x) (0.25-x)/tan(alp);

cf = @(x,y,t,s) s*sqrt(rc^2-(t-x).^2) + y;

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

u_res = 200;
v_res = 2000;

su = linspace(-pi/n, pi/n, u_res);
sv = linspace(-pi/n, pi/n, v_res);

Rot_field = Rot(su);

[u,v] = meshgrid(1:u_res, 1:v_res);

[U, V] = meshgrid(su,sv);

X = rack(sv);
X1 = [reshape(R*ones(u_res,v_res),1,u_res*v_res); reshape(-R*U,1,u_res*v_res)];

MRot=Rot_field(:,:,u);

M1=X(:,v);

Y = zeros(2, u_res, v_res);



for i=1:u_res
  for j=1:v_res
    k = (i-1)*v_res + j;
    Y(:,i,j) = MRot(:,:,k) * (M1(:,k)+X1(:,k));
  end
end

df_dn = @(f, n, u, v, d=10^(-6)) ...
      ( f(u + d*(n==1)/2,v + d*(n==2)/2) - f(u - d*(n==1)/2,v - d*(n==2)/2) ) / d;

rack = @(t) [ftooth(abs(t)); t];

y = @(u,v) resize((Rot(u)*(rack(v)+[R; -R*u])),3,1);

F = @(u, v) cross(df_dn(y, 1, u, v), df_dn(y, 2, u, v))(3);

axis equal;
hold on;
warning("off","all");

for i = 1:u_res
    plot(Y(1,i,:),Y(2,i,:), "k");
end


points = zeros(2,v_res);

for i = 1:size(sv,2)
  u=0;
  solu = fsolve(@(u) F(u,sv(i)), u);
  points(:,i) = y(solu,sv(i))(1:2);
endfor
plot(points(1,:),points(2,:), "r", "LineWidth", 3, "marker", "none");



filtered_points = zeros(2,size(points,2));
l=1;
for i=1:size(points,2)

  field = reshape(vecnorm(Y - points(:,i), 2, 1), u_res, v_res);
  threshold = 0.001;
  zero_mask = field < threshold;
  [labels, num_regions] = bwlabel(zero_mask, 8);
  num_regions
  if num_regions < 2
    filtered_points(:,l) = points(:,i);
    l++;
  endif
end
filtered_points = resize(filtered_points, 2, l-1);
plot(filtered_points(1,:),filtered_points(2,:), "g", "LineWidth", 3, "marker", "none");

hold off;
