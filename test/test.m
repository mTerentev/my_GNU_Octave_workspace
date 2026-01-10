R = 1;
ro = 1.4;
ri = 0.6;
alp = 20*pi/180;
n = 5;
rc = 0.1;

tup=(ro-R)*n/(pi*R);
tdown=(ri-R)*n/(pi*R);
x1 = 0.5-(tup-rc)*tan(alp)-rc/cos(alp);
x2 = x1+rc*cos(alp);
x4 = 0.5-(tdown+rc)*tan(alp)+rc/cos(alp);
x3 = x4-rc*cos(alp);

tf = @(x) (0.5-x)/tan(alp);
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

rack = @(t) [reshape(ftooth(abs(t)),1,s(t)); reshape(t,1,s(t))]*pi*R/n;

u_res = 5000;
v_res = 5000;

su = linspace(-2*pi/n, 2*pi/n, u_res);
sv = linspace(-1, 1, v_res);
[u,v] = meshgrid(1:u_res, 1:v_res);
[U, V] = meshgrid(su,sv);


d = @(s) (s(end)-s(1))/size(s,2);
Rot_field = Rot(su);

X = rack(sv);
X1 = reshape([reshape(R*ones(u_res,v_res),1,u_res*v_res); reshape(-R*U,1,u_res*v_res)],2,1,v_res,u_res);

MRot=reshape(Rot_field(:,:,u),2,2,v_res,u_res);

M1=reshape(X(:,v),2,1,v_res,u_res);

function result = batchMTimesV(A,B)
  result = zeros(2, 1, size(A,3), size(A,4));
  result(1,1,:,:) = A(1,1,:,:).*B(1,1,:,:) + A(1,2,:,:).*B(2,1,:,:);
  result(2,1,:,:) = A(2,1,:,:).*B(1,1,:,:) + A(2,2,:,:).*B(2,1,:,:);
endfunction

function result = batchCross(u,v)
  result = zeros(1, 1, size(u,3), size(u,4));
  result(1,1,:,:) = u(1,1,:,:).*v(2,1,:,:) - u(2,1,:,:).*v(1,1,:,:);
endfunction

Y = batchMTimesV(MRot, M1 + X1);

dY_du = resize(diff(Y,1,4)/d(su),2,1,v_res-1,u_res-1);
dY_dv = resize(diff(Y,1,3)/d(sv),2,1,v_res-1,u_res-1);
F = batchCross(dY_du, dY_dv);

axis equal;
hold on;

for i=1:30:u_res
  %plot(Y(1,1,:,i),Y(2,1,:,i),"k")
endfor

points = zeros(2,1,v_res-1);
for i = 1:v_res-1
  [sol, ind] = min(abs(F(1,1,i,:)));
  points(:,:,i) = Y(:,:,i,ind);
endfor

filtered_points = zeros(2,1,v_res-1);
l=1;
k=1;
while true

  [sol, ind] = find(vecnorm(points-points(:,:,k)) < 0.01, 1, "last");
  if ind - k > 5
    k += ind-k;
  endif
  filtered_points(:,:,l) = points(:,:,k);
  l++;
  k++;
  if k > v_res-1 break; endif
endwhile
filtered_points = resize(filtered_points, 2, 1, l-1);


%plot(filtered_points(1,1,:), filtered_points(2,1,:), "g", "LineWidth", 3);

gear = zeros(2,1,size(filtered_points,3),n);
for i=1:n
  gear(:,:,:,i) = Rot(2*pi/n*i)*filtered_points(:,:,:);
endfor

gear = cat(1,gear);
plot(gear(1,:,:), gear(2,:,:));
hold off;
