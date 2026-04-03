function eta_ratio = Hansen_Viscosity(M0, Mi, x_all, Sij)
% Hansen_Viscosity.m  (FINAL, Hansen-consistent)
% --------------------------------------------------------------
% Computes eta_ratio = eta/eta0 using Hansen Table-V lumping style.
%
% Species order: [O2, N2, O, N, O+, N+, e-]
%
% Hansen coupling rules for viscosity:
%   - neutral–neutral        YES
%   - neutral–ion            YES
%   - ion–ion                YES
%   - ion–electron           YES
%   - electron–ion           YES
%   - electron–electron      YES
%   - neutral–electron       NO
%   - electron–neutral       NO
%
% Important:
%   - NaNs in Sij are treated as "not available" → 0 contribution
%   - No artificial 1e-30 inversions
%   - If s_i → 0, that species contributes 0
% --------------------------------------------------------------

% -------------------- checks --------------------
if ~isscalar(M0) || ~isfinite(M0) || M0 <= 0
    error('Hansen_Viscosity:BadM0','M0 must be a positive finite scalar.');
end

Mi_vec = Mi(:);
x_vec  = x_all(:);

if numel(Mi_vec) ~= 7 || numel(x_vec) ~= 7
    error('Hansen_Viscosity:BadSize','Mi and x_all must each have 7 entries.');
end
if ~isequal(size(Sij),[7 7])
    error('Hansen_Viscosity:BadSij','Sij must be 7x7.');
end
if any(~isfinite(Mi_vec)) || any(Mi_vec <= 0)
    error('Hansen_Viscosity:BadMi','Mi must be positive and finite.');
end

% -------------------- sanitize mole fractions --------------------
x_vec(~isfinite(x_vec)) = 0;
x_vec = max(x_vec,0);
sx = sum(x_vec);
if sx <= 0
    x_vec = [0.2; 0.8; 0; 0; 0; 0; 0];
else
    x_vec = x_vec / sx;
end

% -------------------- sanitize Sij --------------------
S = Sij;
S(~isfinite(S)) = 0;   % NaNs → 0 contribution

% -------------------- active species --------------------
x_cut = 1e-14;
active = (x_vec > x_cut);

if ~any(active)
    eta_ratio = NaN;
    return;
end

% -------------------- mass factor A(i,j) --------------------
Mi_over_Mj = Mi_vec ./ (Mi_vec.');
A = sqrt( (1 + Mi_over_Mj) / 2 );

% -------------------- Hansen coupling mask --------------------
% neutrals: 1–4
% ions:     5–6
% electron: 7

mask = true(7,7);

% neutral–electron OFF
mask(1:4,7) = false;

% electron–neutral OFF
mask(7,1:4) = false;

% everything else ON

% -------------------- compute lambda ratios --------------------
lambda_ratio = zeros(7,1);

for i = 1:7
    if ~active(i), continue; end

    s_i = 0.0;
    for j = 1:7
        if ~active(j), continue; end
        if ~mask(i,j), continue; end

        sij = S(i,j);
        if sij == 0, continue; end

        s_i = s_i + x_vec(j) * sij * A(i,j);
    end

    % SAFE handling: do not invert tiny values
    if s_i <= 0 || ~isfinite(s_i)
        lambda_ratio(i) = 0.0;
    else
        lambda_ratio(i) = 1.0 / s_i;
    end
end

% -------------------- mixture sum --------------------
eta_ratio = 0.0;
for i = 1:7
    if ~active(i), continue; end
    if lambda_ratio(i) == 0, continue; end

    eta_ratio = eta_ratio + x_vec(i) * sqrt(Mi_vec(i)/M0) * lambda_ratio(i);
end

if ~isfinite(eta_ratio)
    eta_ratio = NaN;
end

end
