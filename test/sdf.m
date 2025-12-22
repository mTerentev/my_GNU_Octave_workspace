% Reshape and use arrayfun
M_reshaped = reshape(M, 3, 3, []);
V_reshaped = reshape(M1, 3, 1, []);
R = zeros(3, 1, size(M_reshaped,3));

for k = 1:size(M_reshaped,3)
    R(:,:,k) = M_reshaped(:,:,k) * V_reshaped(:,:,k);
end

R = reshape(R, 3, 1, 101, 101);
