function gear = Gear(n, R, rack, u_res, v_res)
global Rot

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

%Solve Envelope
points = zeros(2,1,v_res-1);
for i = 1:v_res-1
  [sol, ind] = min(abs(F(1,1,i,:)));
  points(:,:,i) = Y(:,:,i,ind);
endfor

%Sparse points evenly
TotalL = 0;
for i = 1: v_res-2
  TotalL += vecnorm(points(:,1,i+1)-points(:,1,i));
endfor

points1 = zeros(2,1,v_res);
l=1;
st_point = points(:,:,1);
for i = 1:1000
  path=TotalL/v_res;
  L = vecnorm(points(:,1,l+1)-st_point);
  while path > L
    l++;
    path -= L;
    st_point = points(:,:,l);
    L = vecnorm(points(:,1,l+1)-st_point);
  endwhile
  st_point = (points(:,1,l+1)-st_point)/L*path + st_point;
  points1(:,:,i) = st_point;
  % if vecnorm(st_point-points(:,:,v_res-1)) < 0.01 break endif
endfor

% plot(points(1,:,:), points(2,:,:));
% plot(points1(1,:,:), points1(2,:,:), "linestyle", "-", "marker", "o");
% waitfor(gcf);
% drawnow;

points = points1;

function res = crs(x,y)
  res = x(1)*y(2)-y(1)*x(2);
endfunction

function bool = isIntersecting(a1,b1,a2,b2)
  % bool = crs(a2-a1, b2-a1);
  bool = (crs(a2-a1, b2-a1) * crs(a2-b1, b2-b1)) < 0 && (crs(a1-a2, b1-a2) * crs(a1-b2, b1-b2)) < 0;
endfunction

function res = Intersection(a1,b1,a2,b2)
  t = inv([(a1-b1),(a2-b2)])*(b2-b1);
  res = a1*t(1)+b1*(1-t(1));
endfunction

%Filter self-intersections out
filtered_points = zeros(2,1,v_res-1);
l=1;
k=1;
while true
  for j = k+10 : k+200
    if j > v_res-2 || j <= 0 break; endif  
    if isIntersecting(points(:,1,k),points(:,1,k+1),points(:,1,j),points(:,1,j+1))
      filtered_points(:,:,l) = points(:,:,k);
      l++;
      filtered_points(:,1,l) = Intersection(points(:,1,k),points(:,1,k+1),points(:,1,j),points(:,1,j+1));
      l++;
      k=j+1;
      break;
    endif
  endfor
  filtered_points(:,:,l) = points(:,:,k);
  l++;
  k++;
  if k > v_res-2 break; endif  
endwhile

filtered_points = resize(filtered_points, 2, 1, l-1);
plot(points(1,:,:), points(2,:,:), "linestyle", "-", "marker", "o");
plot(filtered_points(1,:,:), filtered_points(2,:,:), "color", "red", "linewidth", 3);
waitfor(gcf);
drawnow;

gear = [];
for i=1:n
  gear = [gear, Rot(2*pi/n*i)*filtered_points(:,:,:)];
endfor
endfunction

