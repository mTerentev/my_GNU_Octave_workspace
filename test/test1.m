R = 1;
ro = 1.2;
ri = 0.8;
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

tr_x = @(t) ones(size(t))*R + t/5;
tr_y = @(t) -(t.*R + t.*t.*R./10);
tr_rot = @(t) t;

transform = {tr_x, tr_y, tr_rot};

su = linspace(-pi, pi, u_res);
sv = linspace(-0.7*n, 1.3*n, v_res);

_gear = CurvilinearGear(n, R, rack, transform, su, sv);

p1 = fill(_gear(1,:),_gear(2,:),"r");
p2 = fill(_gear(1,:),_gear(2,:),"b");

%p3 = plot([0,0],[-2,2],"g", "linewidth", 3);
%p4 = plot([0,0],[-2,2],"c", "linewidth", 3);

frames = 1000;

while 1
  for i=1:frames
    t = i/frames;
    c = 2*pi;
    w = 2*pi;
    alp1 = -w*t;
    lin1 = -transform{2}(-alp1-pi) + transform{2}(-pi);
    lin2 = c - lin1
    alp2 = fsolve(@(ang) -transform{2}(-ang) + transform{2}(-pi) - lin2, ang=0);

    x = -transform{1}(-alp1-pi) + R;

    %set(p3, 'xdata', [x,x]);

    x2 = transform{1}(-alp2);

    %set(p4, 'xdata', [x2,x2]);

    gear1 = Rot(alp1)*_gear;
    gear2 = Rot(alp2-0.18)*_gear;

    gear1 += [R; 0];
    gear2 += [-(x2-x);0];
    %gear2 += [-R;0];

    set(p1, "xdata", gear1(1,:), "ydata", gear1(2,:));
    set(p2, "xdata", gear2(1,:), "ydata", gear2(2,:));

    drawnow;
  endfor
endwhile

