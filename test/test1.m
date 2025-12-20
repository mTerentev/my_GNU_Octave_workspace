
r = @(u,v) [u*v; u^2 + v^2; sin(u) + cos(v)];

% Creates function returning 3x1 vector
f = matlabFunction(sym(r));

% Use it:
result = f(2, 3:4)  % Returns 3x1 vector