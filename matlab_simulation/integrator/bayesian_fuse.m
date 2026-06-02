function p_fused = bayesian_fuse(p_mi, p_cvsa, t_since_reset, half_life)
% BAYESIAN_FUSE  Hybrid fusion: cosine-annealed LOP (no plateau).
%
%   alpha(t):
%     0 <= t < half_life  →  alpha = 0.5*(1+cos(pi*t/half_life))
%                               alpha(0) = 1  (full CVSA influence)
%                               alpha(half_life/2) = 0.5  (equal weight)
%                               alpha(half_life) = 0  (pure MI)
%     t >= half_life       →  alpha = 0  (pure MI)
%
%   LOP:
%     prior(c) ∝ P_CVSA(c)^alpha
%     lop(c)   ∝ P_MI(c) * prior(c)   (normalised)
%
%   Behaviour:
%     agree     → LOP boosts above both inputs
%     disagree  → products nearly cancel → near-uniform
%     alpha = 0 → prior = uniform → pure MI
%
%   p_mi, p_cvsa  [1 x n]   sLDA probabilities
%   t_since_reset           seconds since the last event 781
%   half_life               total decay duration in seconds
%   p_fused       [1 x n]   posterior (sums to 1)

    if t_since_reset < half_life
        alpha = 0.5 * (1 + cos(pi * t_since_reset / half_life));
    else
        alpha = 0.0;
    end

    prior = p_cvsa .^ alpha;
    prior = prior / sum(prior);
    lop   = p_mi .* prior;
    lop   = lop  / sum(lop);

    p_fused = lop;
end
