% graphics_toolkit("qt5");

pkg load ocl;

R = 1;
ro = 1.25;
ri = 0.75;
alp = 20*pi/180;
n = 8;
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

global Rot
Rot = @(alp) reshape([
    reshape(cos(alp),1,1,s(alp)) ,reshape(-sin(alp),1,1,s(alp));
    reshape(sin(alp),1,1,s(alp)), reshape(cos(alp),1,1,s(alp));
  ],2,2,s(alp));

rack = @(t) [reshape(ftooth(abs(t)),1,s(t)); reshape(t,1,s(t))]*pi*R/n;



u_res = 5000;
v_res = 7000;

su = linspace(-2*pi/n, 2*pi/n, u_res);
sv = linspace(-1,1,v_res);

tr_x = @(t) R*cos(t)-(-R.*t).*sin(t);
tr_y = @(t) R*sin(t)+(-R.*t).*cos(t);
tr_rot = @(t) t;

transform = {tr_x, tr_y, tr_rot};

global gear11
gear11 = CurvilinearGear(n, rack(sv), transform, su, sv);

global arr_size
arr_size = size(gear11, 2)

"Interpolate"
function y = gear_func(t)
  global arr_size
  global gear11
  i = floor(arr_size.*t)+1;
  q = ones(2,size(t,1))*(arr_size*t - i + 1);
  m1 = (gear11(:,i) + gear11(:,i+1)) / 2;
  m2 = (gear11(:,i+1) + gear11(:,i+2)) / 2;
  p1 = m1.*(1-q) + gear11(:,i+1).*q;
  p2 = gear11(:,i+1).*(1-q) + m2.*q;
  y = p1.*(1-q) + p2.*q;
endfunction

tls = 0:0.0001:0.999;
hold on;
plot(gear_func(tls)(1,:),gear_func(tls)(2,:), "linestyle", "-", "marker", ".");
drawnow;
waitfor(gcf);

% u_res = 3000;
% v_res = 10000;

su = linspace(-pi/n, pi/n, u_res);
sv = linspace(1.5/n,0.5/n,v_res);

tr_x = @(t) 3*R*cos(t);
tr_y = @(t) 3*R*sin(t);
tr_rot = @(t) 3*t - ones(size(t))*(pi+2*pi/n/2 - 2*2*pi/n);
transform = {tr_x, tr_y, tr_rot};

gear21 = CurvilinearGear(2*n, gear_func(sv), transform, su, sv);


rack11 = [];
for i=-3:3
  rack11 = [rack11, rack(-1:0.01:1)+[0;i*2*pi*R/n]];
endfor
rack11 += [2*R; 0];
rack11 = [rack11, [3*R;3.5*2*pi*R/n]];
rack11 = [rack11, [3*R;-3.5*2*pi*R/n]];
rack11 = [rack11, rack11(:,1)];

function rgb = hex2rgb(hex)
    hex = strrep(hex, '#', '');
    rgb = sscanf(hex, '%2x%2x%2x', [1, 3]) / 255;
end

axis equal;
hold on;
axis([-6*R 6*R -3*R 3*R]);

p1 = plot(gear11(1,:,:), gear11(2,:,:), hex2rgb("#29D69A"));
p2 = plot(gear21(1,:,:), gear21(2,:,:), hex2rgb("#D69A29"));
p3 = plot(rack11(1,:),rack11(2,:), hex2rgb("#9A29D6"));

fps = 120;
tme = 2;
rec1 = zeros(fps,2,size(gear11,2));
rec2 = zeros(fps,2,size(gear21,2));
rec3 = zeros(fps,2,size(rack11,2));

for i=1:(fps*tme)
    alp1 = i/fps/tme*2*pi/n;
    alp2 = -0.5*alp1 + pi/n/2 + pi/2;

    gear1 = Rot(alp1)*gear11;
    gear2 = Rot(alp2)*gear21;

    gear1 += [R; 0];
    gear2 += [-2*R; 0];
    rack1 = rack11 + [0; alp1*R];

    rec1(i,:,:) = gear1;
    rec2(i,:,:) = gear2;
    rec3(i,:,:) = rack1;
  endfor

while true
  for i=1:(fps*tme)
    set(p1, "XData", rec1(i,1,:)(:), "YData", rec1(i,2,:)(:));
    set(p2, "XData", rec2(i,1,:)(:), "YData", rec2(i,2,:)(:));
    set(p3, "XData", rec3(i,1,:)(:), "YData", rec3(i,2,:)(:));
    drawnow;
    pause(1/fps);
    clc;

  endfor
endwhile
