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

u_res = 100;
v_res = 1000;

%axis fix;
axis equal;
hold on;
sv = linspace(-n/2, n/2, v_res);

gear = CurvilinearGear(n, R, rack, u_res, v_res);

%fill(gear(1,:),gear(2,:),"r");

