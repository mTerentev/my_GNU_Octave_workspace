function gear = CurvilinearGear(n, tooth, transform_func, su, sv)
global Rot

axis equal;
hold on;

"Solve Envelope"

u_res = size(su,2)
v_res = size(sv,2)

d = @(s) (s(end)-s(1))/size(s,2);

Rot_field = Rot(transform_func{3}(-su));
size(Rot_field);
X = tooth;

MRot=Rot_field;
X1 = reshape([transform_func{1}(su)(:)'; transform_func{2}(su)(:)'],[2,u_res]);
M1=reshape(X(:,1:v_res),2,v_res);

function C = batchMTimesV(A, v)

  A = oclArray(A);
  v = oclArray(v);

  A = reshape(A,[2,2*u_res]);
  C = A'*v;

  C = single(C);

  C = reshape(C,[2,u_res, v_res]);

end

function result = batchCross(u,v)
  u1 = oclArray(reshape(u(1,:,:),[u_res-1,v_res-1]));
  v2 = oclArray(reshape(v(2,:,:),[u_res-1,v_res-1]));
  result = single(u1.*v2);
  u2 = oclArray(reshape(u(2,:,:),[u_res-1,v_res-1]));
  v1 = oclArray(reshape(v(1,:,:),[u_res-1,v_res-1]));
  result -= single(u2.*v1);
endfunction

"Solve Y"
tic
Y = batchMTimesV(MRot, M1);
Y += X1;
toc


for i = 1:floor(u_res/30):u_res
  plot(Y(1,i,:), Y(2,i,:), "black");
endfor

%todo: accelerate diff()

"Solve F"
Y1 = oclArray(Y);
dY_du = resize(diff(Y1,1,3)/d(su),2,u_res-1,v_res-1);
dY_dv = resize(diff(Y1,1,2)/d(sv),2,u_res-1,v_res-1);
F = batchCross(dY_du, dY_dv);

"Solve envelope"
points = zeros(2,v_res-1);
prev_ind = 0;
for i = 1:v_res-1
  arr = abs(F(:,i));
  do
    [sol, ind] = min(arr);
    arr(ind) = 200;

  until abs(ind - prev_ind) < u_res/10 || prev_ind == 0;
  prev_ind = ind;
  points(:,i) = Y(:,ind,i);
endfor

plot(points(1,:,:), points(2,:,:), "linestyle", "-", "marker", "o", "color", "blue");

"Sparse points evenly"

Np = 1000;

TotalL = 0;
for i = 1: v_res-2
  TotalL += vecnorm(points(:,i+1)-points(:,i));
endfor

points1 = zeros(2,v_res);
l=1;
q=1;
st_point = points(:,1);
for i = 1:v_res
  path=TotalL/Np;
  L = vecnorm(points(:,l+1)-st_point);
  while path > L && l < v_res-2
    l++;
    path -= L;
    st_point = points(:,l);
    L = vecnorm(points(:,l+1)-st_point);
  endwhile
  st_point = (points(:,l+1)-st_point)/L*path + st_point;
  points1(:,i) = st_point;
  q = i;
  if vecnorm(st_point-points(:,v_res-1)) < 0.01; break endif
endfor

points1 = resize(points1, 2, q-1);

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

filtered_points = zeros(2,q-1);
l=1;
k=1;
while true
  mask = reshape(abs(points(1,:)-points(1,k)) < TotalL/Np*2  & abs(points(2,:)-points(2,k)) < TotalL/Np*2, q-1, 1);

  ind = find(mask, 1, "last");

  if ind - k > 20

  for j = ind-20 : ind
    if j > q-2 || j <= 0 break; endif  
    if isIntersecting(points(:,k),points(:,k+1),points(:,j),points(:,j+1))
      filtered_points(:,l) = points(:,k);
      l++;
      filtered_points(:,l) = Intersection(points(:,k),points(:,k+1),points(:,j),points(:,j+1));
      l++;
      k=j+1;
      break;
    endif
  endfor

  endif
  filtered_points(:,l) = points(:,k);
  l++;
  k++;
  if k > q-2 break; endif  
endwhile

filtered_points = resize(filtered_points, 2, l-1);


% plot(points(1,:,:), points(2,:,:), "linestyle", "-", "marker", "o", "color", "red");
plot(filtered_points(1,:), filtered_points(2,:), "color", "red", "linewidth", 3);


waitfor(gcf);
drawnow;

gear = [];
for i=1:n
  gear = [gear, Rot(2*pi/n*i)*filtered_points(:,:,:)];
endfor

gear = [gear, gear(:,1)];

endfunction
