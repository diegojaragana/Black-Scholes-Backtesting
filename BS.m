
function V0 = BS(S_i, K, r, q, T, d1, d2, epsilon) 
V0 = epsilon * S_i * exp(-q * T) * normcdf(epsilon * d1) - epsilon * K * exp(-r * T) * normcdf(epsilon * d2);
end