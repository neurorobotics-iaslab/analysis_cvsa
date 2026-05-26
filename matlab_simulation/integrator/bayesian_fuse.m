function p_fused = bayesian_fuse(p_mi, p_cvsa, t_since_reset, half_life)
% BAYESIAN_FUSE  Hybrid fusion with agreement gate, identical to Integrator.cpp:
%
%   Step 1 — LOP (Logarithmic Opinion Pool):
%     alpha(t) = 0.5 * (1 + cos(pi * min(t, half_life) / half_life))
%     prior(c) ∝ P_CVSA(c)^alpha   (tempered CVSA prior)
%     lop(c)   ∝ P_MI(c) * prior(c)
%
%   Step 2 — Agreement gate:
%     agree_raw = dot(P_CVSA, P_MI)
%     agree_w   = max(0, (agree_raw - 1/n) / (1 - 1/n))   in [0,1]
%     neutral_w = (1 - agree_w) * alpha
%     p_fused   = (1 - neutral_w) * lop + neutral_w * (1/n)
%
%   Behaviour:
%     agree  + any alpha  → neutral_w ≈ 0 → pure LOP  (amplifies agreement)
%     disagree + alpha=1  → neutral_w = 1 → uniform   (blocks buffer at t=0)
%     disagree + alpha=0  → neutral_w = 0 → pure MI   (CVSA faded out)
%
%   p_mi, p_cvsa  [1 x n]   sLDA probabilities
%   t_since_reset           seconds since the last event 781
%   half_life               seconds for CVSA prior to decay to uniform (2.5)
%   p_fused       [1 x n]   posterior (sums to 1)

    alpha = 0.5 * (1 + cos(pi * min(t_since_reset, half_life) / half_life));

    % --- LOP ---
    prior = p_cvsa .^ alpha;
    prior = prior / sum(prior);
    lop   = p_mi .* prior;
    lop   = lop  / sum(lop);

    % --- Agreement gate ---
    n_cls     = numel(p_mi);
    agree_raw = sum(p_cvsa .* p_mi);
    uniform_p = 1.0 / n_cls;
    agree_w   = max(0, (agree_raw - uniform_p) / (1.0 - uniform_p));
    neutral_w = (1 - agree_w) * alpha;
    p_fused   = (1 - neutral_w) * lop + neutral_w * uniform_p;
    % Stays normalised: convex combination of two distributions that sum to 1.
end
