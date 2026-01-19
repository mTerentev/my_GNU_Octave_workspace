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
rack1 = @(t) rack(t);
rack2 = @(t) [-rack(t)(1,:); rack(t)(2,:)];



u_res = 5000;
v_res = 5000;

axis equal;
hold on;
axis([-3*R 3*R -2*R 2*R]);

tr_x = @(t) ones(size(t))*R + t/5;
tr_y = @(t) -(t.*R + t.*t.*R./10);
tr_rot = @(t) t;

transform = {tr_x, tr_y, tr_rot};

su = linspace(-pi, pi, u_res);
sv = linspace(-0.7*n, 1.22*n, v_res);

%_gear1 = CurvilinearGear(n, R, rack1, transform, su, sv);
%_gear2 = CurvilinearGear(n, R, rack2, transform, su, sv);

p1 = fill(_gear1(1,:),_gear1(2,:),"r");
p2 = fill(_gear2(1,:),_gear2(2,:),"b");

_rack1 = rack1(sv);
_rack2 = rack2(sv);

%p3 = plot(_rack1(1,:),_rack1(2,:),"g", "linewidth", 3);
%p4 = plot(_rack2(1,:),_rack2(2,:),"c", "linewidth", 3);

frames = 1000;

beta = atan(0.2);

rec1 = zeros([frames, size(_gear1)]);
rec2 = zeros([frames, size(_gear2)]);

for i=1:frames-10
  disp(i);
  t = i/frames;
  w = 2*pi;
  alp1 = -w*t;
  lin1 = -transform{2}(-alp1-pi) + transform{2}(-pi) + 0.09;
  alp2 = fsolve(@(ang) +transform{2}(-ang) - transform{2}(pi) - lin1, ang=0);

  x = -transform{1}(-alp1-pi-beta) + R;

  x2 = transform{1}(-alp2);

  gear1 = Rot(alp1)*_gear1;
  gear2 = Rot(alp2)*_gear2;

  gear1 += [R; 0];
  gear2 += [-(x2-x);0];

  rec1(i,:,:,:) = gear1;
  rec2(i,:,:,:) = gear2;

endfor


while 1
  for i=1:frames-10
    set(p1, "xdata", rec1(i,1,:,:)(:), "ydata", rec1(i,2,:,:)(:));
    set(p2, "xdata", rec2(i,1,:,:)(:), "ydata", rec2(i,2,:,:)(:));
    drawnow;
    pause(0.01);
  endfor
endwhile

