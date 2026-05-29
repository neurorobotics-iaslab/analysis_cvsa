function p_fused = bayesian_fuse(p_mi, p_cvsa, t_since_reset, t_hold, half_life)
% BAYESIAN_FUSE  Hybrid fusion: plateau + cosine-annealed LOP.
%
%   alpha(t):
%     t <= t_hold                     → alpha = 1   (full CVSA influence)
%     t_hold < t <= t_hold + half_life → alpha = 0.5*(1+cos(pi*(t-t_hold)/half_life))
%     t > t_hold + half_life           → alpha = 0   (pure MI)
%
%   LOP:
%     prior(c) ∝ P_CVSA(c)^alpha
%     lop(c)   ∝ P_MI(c) * prior(c)   (normalised)
%
%   Behaviour:
%     agree  → LOP boosts above both inputs
%     disagree (symmetric) → products cancel → uniform naturally
%     alpha = 0 → prior = uniform → pure MI
%
%   p_mi, p_cvsa  [1 x n]   sLDA probabilities
%   t_since_reset           seconds since the last event 781
%   t_hold                  seconds CVSA stays at full influence (plateau)
%   half_life               seconds for the cosine decay after the plateau
%   p_fused       [1 x n]   posterior (sums to 1)

    if t_since_reset <= t_hold
        alpha = 1.0;
    elseif t_since_reset <= t_hold + half_life
        t_decay = t_since_reset - t_hold;
        alpha = 0.5 * (1 + cos(pi * t_decay / half_life));
    else
        alpha = 0.0;
    end

    prior = p_cvsa .^ alpha;
    prior = prior / sum(prior);
    lop   = p_mi .* prior;
    lop   = lop  / sum(lop);

    p_fused = lop;
end
