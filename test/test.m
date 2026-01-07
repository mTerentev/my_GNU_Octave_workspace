


R = 1;
ro = 1.4;
ri = 0.6;
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

rack = @(t) [reshape(ftooth(abs(t)),1,s(t)); reshape(t,1,s(t))]*2*pi*R/n;

vrshp = @(v) [reshape(v(1),1,s(v)); reshape(v(2),1,s(v))];

u_res = 3000;
v_res = 3000;

su = linspace(-pi/n*2, pi/n*2, u_res);
sv = linspace(-pi/n, pi/n, v_res);
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
  plot(Y(1,1,:,i),Y(2,1,:,i),"k")
endfor

sv1 = sv;
points = zeros(2,1,v_res-1);

for i = 1:v_res-1
  [sol, ind] = min(abs(F(1,1,i,:)));
  points(:,:,i) = Y(:,:,i,ind);

##  field = zeros(size(Y));
##  field(1,1,:,:) = Y(1,1,:,:) - points(1,1,i);
##  field(2,1,:,:) = Y(2,1,:,:) - points(2,1,i);
##  threshold = 0.01;
##  zero_mask = field < threshold;

##  [x,y] = find(zero_mask,1);
##  dist = abs(x-i)+abs(y-ind)
##  if dist < 100
##    filtered_points(:,:,l) = points(:,:,i);
##    l++;
##  endif
endfor
disp("filtering");

[unique_points, ia, ic] = uniquetol(reshape(Y,2,v_res*u_res), 1e-6, 'ByRows', true);

% Step 3: Count multiplicity
counts = accumarray(ic, 1);

% Step 4: Classify points
boundary_idx = find(counts == 1);  % Unique mappings
interior_idx = find(counts > 1);   % Multiple mappings

boundaryPoints = unique_points(boundary_idx, :);
interiorPoints = unique_points(interior_idx, :);

% Step 5: Refine boundary using α-shape
if ~isempty(boundaryPoints)
    filtered_points = alphaShape(boundaryPoints(:,1), boundaryPoints(:,2), 0.1);
endif

%{
filtered_points = zeros(2,1,size(points,3));
l=1;

uniquetol()

for i=1:v_res-1
  if vecnorm(points(:,1,i),2,1) < ro+0.0001
    filtered_points(:,1,l) = points(:,1,i);
    l++;
  endif
  field = Y-points(:,1,i);
endfor
filtered_points = resize(filtered_points,2,1,l-1);
%}
gear = zeros(2,1,size(filtered_points,3),n);

for i=1:n
  gear(:,:,:,i) = Rot(2*pi/n*i)*filtered_points(:,:,:);
endfor

gear = cat(1,gear);



plot(gear(1,1,:),gear(2,1,:), "g", "LineWidth", 3, "marker", "none");
plot(points(1,1,:),points(2,1,:), "r", "LineWidth", 3, "marker", "none");

hold off;
