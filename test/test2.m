global Rot
s = @(x) size(x,2);
Rot = @(alp) reshape([
    reshape(cos(alp),1,1,s(alp)) ,reshape(-sin(alp),1,1,s(alp));
    reshape(sin(alp),1,1,s(alp)), reshape(cos(alp),1,1,s(alp));
  ],2,2,s(alp));

function C = batchMTimesV(A, v)
% A': (rows, cols, batch1) 3D array
% v: (rows, batch2) vector
% C: (rows, batch1, batch2) result

  batch1 = size(A,3);
  batch2 = size(v,2);

  A = reshape(A,[2,2*batch1]);
  C = A'*v;

  % C = permute(reshape(C,[2,batch1, batch2]),[1,3,2]);
  C
  C = reshape(C,[2,batch1, batch2])


end

A = [Rot(0:5)]
v = [[1;2],[3;4],[5;6]]

b = batchMTimesV(A,v);