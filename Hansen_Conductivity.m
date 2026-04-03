function [k_ratio, k_SI, k_n_ratio, k_r_ratio] = Hansen_Conductivity( ...
    T, R, M0, Mi, x_all, Sij, Spr, C_over_R_vec, ...
    dlnKp_O2to2O, dlnKp_N2to2N, dlnKp_Oion, dlnKp_Nion)
% Hansen_Conductivity.m  (Eq.79 diagonal-only + Eq.84 double-sum)
% -------------------------------------------------------------------------
% Goal: reproduce the stable "old" Hansen behavior (no 85–87),
% using:
%   - Neutral term: Hansen Eq.(79) with DIAGONAL-only mean-free-path model
%   - Reactive term: Hansen Eq.(84) with Spr = S'_ij/S0 (Table V lumping)
%
% Species order: [O2 N2 O N O+ N+ e]
% -------------------------------------------------------------------------

tiny = 1e-30;

Mi_vec = Mi(:);
x_vec  = x_all(:);
C_vec  = C_over_R_vec(:);

% ---- sanitize mole fractions ----
x_vec(~isfinite(x_vec)) = 0;
x_vec = max(x_vec, 0);
sx = sum(x_vec);
if sx <= 0
    x_vec = [0.2; 0.8; 0; 0; 0; 0; 0];
    sx = 1.0;
end
x_vec = x_vec / sx;

% ---- sanitize matrices ----
Sij(~isfinite(Sij)) = 0;
Spr(~isfinite(Spr)) = 0;

% ==============================================================
% 1) Neutral conductivity k_n/k0  (Hansen Eq.79 + Table V diagonal model)
% ==============================================================
% DIAGONAL-ONLY lambda_i/lambda_0
lambda_ratio = zeros(7,1);
for i = 1:7
    Sii = max(Sij(i,i), tiny);         % diagonal only
    lambda_ratio(i) = sqrt(M0/Mi_vec(i)) * 1.0/sqrt(Sii);
end

% Hansen Eq.80 style sum (what you used before)
k_n_ratio = 0.0;
for i = 1:7
    Ci_factor = (4/19)*C_vec(i) + 9/19;
    k_n_ratio = k_n_ratio + x_vec(i) * lambda_ratio(i) * Ci_factor;
end

% ==============================================================
% 2) Reactive conductivity k_r/k0  (Hansen Eq.84, per reaction)
% ==============================================================
a_O2to2O = [-1;  0;  2;  0; 0; 0; 0];
a_N2to2N = [ 0; -1;  0;  2; 0; 0; 0];
a_Oion   = [ 0;  0; -1;  0; 1; 0; 1];  % O -> O+ + e
a_Nion   = [ 0;  0;  0; -1; 0; 1; 1];  % N -> N+ + e

k_r_ratio = 0.0;
k_r_ratio = k_r_ratio + reaction_term_eq84(a_O2to2O, dlnKp_O2to2O, T, M0, Mi_vec, x_vec, Spr, tiny);
k_r_ratio = k_r_ratio + reaction_term_eq84(a_N2to2N, dlnKp_N2to2N, T, M0, Mi_vec, x_vec, Spr, tiny);
k_r_ratio = k_r_ratio + reaction_term_eq84(a_Oion,   dlnKp_Oion,   T, M0, Mi_vec, x_vec, Spr, tiny);
k_r_ratio = k_r_ratio + reaction_term_eq84(a_Nion,   dlnKp_Nion,   T, M0, Mi_vec, x_vec, Spr, tiny);

% ==============================================================
% 3) Total
% ==============================================================
k_ratio = k_n_ratio + k_r_ratio;
k_SI = [];

if ~isfinite(k_ratio)
    k_ratio = NaN;
end

end


% -------------------------------------------------------------------------
function kr = reaction_term_eq84(a_vec, dlnKp_dT, T, M0, Mi_vec, x_vec, Spr, tiny)
% One-reaction contribution to k_r/k0 using Hansen Eq.(84) double sum.
%
% NOTE: This form is the "old" one that reproduces your stable peak because
% it clamps the denominator away from zero.

if ~isfinite(dlnKp_dT) || abs(dlnKp_dT) < 1e-20
    kr = 0.0;
    return;
end

TD = T * dlnKp_dT;
if ~isfinite(TD)
    kr = 0.0;
    return;
end

denom = 0.0;

for i = 1:7
    ai = a_vec(i);
    if ai == 0, continue; end

    xi = max(x_vec(i), tiny);
    Mi_i = Mi_vec(i);

    for j = 1:7
        aj = a_vec(j);
        if aj == 0 && ai == 0
            continue;
        end

        xj = x_vec(j);
        Mi_j = Mi_vec(j);

        Spr_ij = Spr(i,j);
        if Spr_ij == 0, continue; end

        mass_factor = sqrt( (Mi_i*Mi_j) / (M0*(Mi_i + Mi_j)) );

        denom = denom + mass_factor * (Spr_ij * ai/xi) * (ai*xj - aj*xi);
    end
end

% critical stabilization (this is what prevented the crazy spikes)
denom = max(denom, tiny);

kr = (12*sqrt(2)/95) * (TD^2) / denom;

if ~isfinite(kr)
    kr = 0.0;
end

end
