function delta_i = getBSdelta(S_i, K, r, q, sigma, T, epsilon)
    %d1 = (log (S_i / K) + (r - q) * T) / (sigma * sqrt(T)) + (sigma * sqrt(T)) / 2;
    d1 = getd1(S_i, K, r, q, sigma, T);
    delta_i = epsilon * exp(-q * T) * normcdf(epsilon * d1);
end