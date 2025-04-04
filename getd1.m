
function d1 = getd1(S_i, K, r, q, sigma, T)
d1 = (log (S_i / K) + (r - q) * T) / (sigma * sqrt(T)) + (sigma * sqrt(T)) / 2;
end