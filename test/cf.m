function f = cf(x,y,t,s,rc)
    f = s*sqrt(rc^2-(t-x).^2) + y;
  end
