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

global Rot
Rot = @(alp) reshape([
    reshape(cos(alp),1,1,s(alp)) ,reshape(-sin(alp),1,1,s(alp));
    reshape(sin(alp),1,1,s(alp)), reshape(cos(alp),1,1,s(alp));
  ],2,2,s(alp));

rack = @(t) [reshape(ftooth(abs(t)),1,s(t)); reshape(t,1,s(t))]*pi*R/n;

u_res = 2000;
v_res = 2000;

gear11 = Gear(n, R, rack, u_res, v_res);
gear21 = Gear(2*n, 2*R, rack, u_res, v_res);
rack11 = [];
for i=-3:3
  rack11 = [rack11, rack(-1:0.01:1)+[0;i*2*pi*R/n]];
endfor
rack11 += [2*R; 0];
rack11 = [rack11, [3*R;3.5*2*pi*R/n]];
rack11 = [rack11, [3*R;-3.5*2*pi*R/n]];

function rgb = hex2rgb(hex)
    hex = strrep(hex, '#', '');
    rgb = sscanf(hex, '%2x%2x%2x', [1, 3]) / 255;
end

p1 = fill(gear1(1,:,:), gear1(2,:,:), hex2rgb("#29D69A"));
p2 = fill(gear2(1,:,:), gear2(2,:,:), hex2rgb("#D69A29"));
p3 = fill(rack1(1,:),rack1(2,:), hex2rgb("#9A29D6"));

while true
  for i=1:1000
    alp1 = i/1000*2*pi/n;
    alp2 = -0.5*alp1;

    gear1 = Rot(alp1)*gear11;
    gear2 = Rot(alp2)*gear21;

    gear1 += [R; 0];
    gear2 += [-2*R; 0];
    rack1 = rack11 + [0; alp1*R];

    set(p1, "XData", gear1(1,:,:), "YData", gear1(2,:,:));
    set(p2, "XData", gear2(1,:,:), "YData", gear2(2,:,:));
    set(p3, "XData", rack1(1,:,:), "YData", rack1(2,:,:));

    drawnow;
  endfor
endwhile
