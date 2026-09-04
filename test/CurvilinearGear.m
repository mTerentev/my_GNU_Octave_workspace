function gear = CurvilinearGear(n, R, tooth, transform_func, su, sv)
global Rot

axis equal;
hold on;

"Solve Envelope"

u_res = size(su,2);
v_res = size(sv,2);

[u,v] = meshgrid(1:u_res, 1:v_res);
[U, V] = meshgrid(su,sv);

d = @(s) (s(end)-s(1))/size(s,2);

Rot_field = Rot(transform_func{3}(su));

X = tooth;

MRot=reshape(Rot_field(:,:,u),2,2,v_res,u_res);
X1 = reshape([transform_func{1}(U)(:)'; transform_func{2}(U)(:)'],[2,1,size(U)]);
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

"Solve Y"
Y = batchMTimesV(MRot, M1);
Y += X1;


for i = 1:10:u_res
  plot(Y(1,1,:,i), Y(2,1,:,i), "black");
endfor

"Solve F"
dY_du = resize(diff(Y,1,4)/d(su),2,1,v_res-1,u_res-1);
dY_dv = resize(diff(Y,1,3)/d(sv),2,1,v_res-1,u_res-1);
F = batchCross(dY_du, dY_dv);

"Solve envelope"
points = zeros(2,1,v_res-1);
prev_ind = 0;
for i = 1:v_res-1
  arr = abs(F(1,1,i,:));
  do
    [sol, ind] = min(arr);
    arr(ind) = 200;

  until abs(ind - prev_ind) < v_res/100 || prev_ind == 0;
  prev_ind = ind;
  points(:,:,i) = Y(:,:,i,ind);
endfor

plot(points(1,:,:), points(2,:,:), "linestyle", "-", "marker", "o", "color", "blue");

"Sparse points evenly"

Np = 1000;

TotalL = 0;
for i = 1: v_res-2
  TotalL += vecnorm(points(:,1,i+1)-points(:,1,i));
endfor

points1 = zeros(2,1,v_res);
l=1;
q=1;
st_point = points(:,:,1);
for i = 1:v_res
  path=TotalL/Np;
  L = vecnorm(points(:,1,l+1)-st_point);
  while path > L && l < v_res-2
    l++;
    path -= L;
    st_point = points(:,:,l);
    L = vecnorm(points(:,1,l+1)-st_point);
  endwhile
  st_point = (points(:,1,l+1)-st_point)/L*path + st_point;
  points1(:,:,i) = st_point;
  q = i;
  if vecnorm(st_point-points(:,:,v_res-1)) < 0.01; break endif
endfor

points1 = resize(points1, 2, 1, q-1);

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

"Filter self-intersections out"

filtered_points = zeros(2,1,q-1);
l=1;
k=1;
while true
  mask = reshape(abs(points(1,:,:)-points(1,:,k)) < TotalL/Np*2  & abs(points(2,:,:)-points(2,:,k)) < TotalL/Np*2, q-1, 1);

  ind = find(mask, 1, "last");

  if ind - k > 20

  for j = ind-20 : ind
    if j > q-2 || j <= 0 break; endif  
    if isIntersecting(points(:,1,k),points(:,1,k+1),points(:,1,j),points(:,1,j+1))
      filtered_points(:,:,l) = points(:,:,k);
      l++;
      filtered_points(:,1,l) = Intersection(points(:,1,k),points(:,1,k+1),points(:,1,j),points(:,1,j+1));
      l++;
      k=j+1;
      break;
    endif
  endfor

  endif
  filtered_points(:,:,l) = points(:,:,k);
  l++;
  k++;
  if k > q-2 break; endif  
endwhile

filtered_points = resize(filtered_points, 2, 1, l-1);


% plot(points(1,:,:), points(2,:,:), "linestyle", "-", "marker", "o", "color", "red");
plot(filtered_points(1,:,:), filtered_points(2,:,:), "color", "red", "linewidth", 3);


waitfor(gcf);
drawnow;

gear = [];
for i=1:n
  gear = [gear, Rot(2*pi/n*i)*filtered_points(:,:,:)];
endfor

gear = [gear, gear(:,1)];

endfunction
