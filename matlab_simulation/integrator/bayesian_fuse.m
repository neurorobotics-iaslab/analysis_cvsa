function p_fused = bayesian_fuse(p_mi, p_cvsa, t_since_reset, half_life)
% BAYESIAN_FUSE  Hybrid fusion identical to test_full_pipeline.m compute_fusion:
%   alpha(t) = 0.5 * (1 + cos(pi * min(t, half_life) / half_life))
%   tempered prior:   P_prior(c) ∝ P_CVSA(c)^alpha
%   fused posterior:  P_out(c)   ∝ P_MI(c) * P_prior(c)
%
%   p_mi, p_cvsa  [1 x 2]   sLDA probabilities (class 1, class 2)
%   t_since_reset seconds since the last event 781
%   half_life     seconds for the CVSA prior to decay to uniform (2.5)
%
%   p_fused       [1 x 2]   posterior over (class 1, class 2)
    alpha = 0.5 * (1 + cos(pi * min(t_since_reset, half_life) / half_life));
    prior = p_cvsa .^ alpha;
    prior = prior / sum(prior);
    fused = p_mi .* prior;
    p_fused = fused / sum(fused);
end
