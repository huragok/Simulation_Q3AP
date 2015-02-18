function idxs = get_idxs(q, Q, order)

q_residual = q;
idxs = zeros(6, 1);
for d = 1 : 6
    idxs(order(d)) = mod(q_residual, Q);
    q_residual = floor(q_residual / Q);
end
idxs = idxs + 1;