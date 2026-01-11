function gear = CurvilinearGear(n, R, rack, u_res, v_res)
global Rot

su = linspace(-pi, pi, u_res);
sv = linspace(-0.7*n, 1.3*n, v_res);
[u,v] = meshgrid(1:u_res, 1:v_res);
[U, V] = meshgrid(su,sv);


d = @(s) (s(end)-s(1))/size(s,2);
Rot_field = Rot(su);

X = rack(sv);

X1 = reshape([reshape(R*ones(v_res,u_res)+U/5,1,u_res*v_res); reshape(-(R.*U+R.*U.*U./10),1,u_res*v_res)],2,1,v_res,u_res);

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
