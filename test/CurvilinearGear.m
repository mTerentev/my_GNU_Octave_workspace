function gear = CurvilinearGear(n, R, rack_func, transform_func, su, sv)
global Rot

u_res = size(su,2);
v_res = size(sv,2);

[u,v] = meshgrid(1:u_res, 1:v_res);
[U, V] = meshgrid(su,sv);

d = @(s) (s(end)-s(1))/size(s,2);

Rot_field = Rot(transform_func{3}(su));

X = rack_func(sv);

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

Y = batchMTimesV(MRot, M1 + X1);

for i=1:u_res
  %plot(Y(1,1,:,i),Y(2,1,:,i), "k");
endfor

dY_du = resize(diff(Y,1,4)/d(su),2,1,v_res-1,u_res-1);
dY_dv = resize(diff(Y,1,3)/d(sv),2,1,v_res-1,u_res-1);
F = batchCross(dY_du, dY_dv);

axis equal;
hold on;


points = zeros(2,1,v_res-1);
for i = 1:v_res-1
  [sol, ind] = min(abs(F(1,1,i,:)));
  points(:,:,i) = Y(:,:,i,ind);
endfor

%plot(points(1,1,:),points(2,1,:), "r", "linewidth", 3);

filtered_points = zeros(2,1,v_res-1);
l=1;
k=1;
while true

  [sol, ind] = find(vecnorm(points-points(:,:,k)) < 10*max(d(sv),d(su)), 1, "last");
  if ind - k > 3
    k += ind-k;
  endif
  filtered_points(:,:,l) = points(:,:,k);
  l++;
  k++;
  if k > v_res-1 break; endif
endwhile
filtered_points = resize(filtered_points, 2, 1, l-1);

gear = filtered_points;
%{
gear = [];
for i=1:n
  gear = [gear, Rot(2*pi/n*i)*filtered_points(:,:,:)];
endfor
%}
endfunction
