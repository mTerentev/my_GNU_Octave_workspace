R = 1;
ro = 1;
ri = 1;
alp = 20*pi/180;
n = 10;
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

rack = @(t) [reshape(ftooth(abs(t-2-2*floor(t/2-0.5))),1,s(t)); reshape(t,1,s(t))]*pi*R/n;

u_res = 3000;
v_res = 3000;

axis equal;
hold on;
axis([-3*R 3*R -2*R 2*R]);

%_gear = CurvilinearGear(n, R, rack, u_res, v_res);

p1 = fill(gear1(1,:),gear1(2,:),"r");
p2 = fill(gear2(1,:),gear2(2,:),"b");

p3 = plot([0,0],[-2,2],"g", "linewidth", 3);
p4 = plot([0,0],[-2,2],"c", "linewidth", 3);

frames = 10000;

while 1
  alp1 = 0;
  alp2 = pi;
  for i=1:frames
    t = i/frames - 0.5;

    x = pi/5 * 2*t;

    set(p3, 'xdata', [x,x]);

    alp1 = 5*x + pi;
    w1 = 2*pi;
    lin = 2*pi*(R-x);

    alp2 = 2*pi*(R-pi/5 * 2*t)/(R+pi/5 * 2*t)/frames;
    %alp2 = -2 * pi * t + 10 * R * log(5 * R + 2 * pi * t)- 4*pi;

    x2 = alp2/5;

    set(p4, 'xdata', [x2,x2]);

    gear1 = Rot(alp1)*_gear;
    gear2 = Rot(-alp2)*_gear;

    gear1 += [R; 0];
    gear2 += [-R-(x2-x);0];

    set(p1, "xdata", gear1(1,:), "ydata", gear1(2,:));
    set(p2, "xdata", gear2(1,:), "ydata", gear2(2,:));

    drawnow;
  endfor
endwhile
