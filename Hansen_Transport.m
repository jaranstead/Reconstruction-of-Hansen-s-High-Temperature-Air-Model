function [Mi, x_all, Sij, Spr, See_over_S0, dbg] = Hansen_Transport( ...
    T, M0, Mmix, rhogas_SI, Rgas_SI, cp_mass_SI, ...
    xO2, xN2, xO, xN, xOplus, xNplus, xe, ...
    dlnKp_O2to2O, dlnKp_N2to2N, dlnKp_Oion, dlnKp_Nion) %#ok<INUSD>
% Hansen_Transport.m (FINAL / NAN-SAFE)
% -------------------------------------------------------------------------
% Hansen-style transport using TABLE V (lumped collision integrals).
%
% Species order:
%   [O2  N2  O  N  O+  N+  e]
%
% Outputs:
%   Mi          [kg/mol]
%   x_all       mole fractions
%   Sij         UNPRIMED  S_ij / S0
%   Spr         PRIMED    S'_ij / S0
%   See_over_S0 electron-electron term = See / S0 (already lnΛ multiplied)
%   dbg         diagnostics
% -------------------------------------------------------------------------

% -------------------- constants (cgs) --------------------
kB_cgs = 1.380649e-16;      % erg/K
e_cgs  = 4.80320425e-10;    % statC
pi_c   = pi;

% -------------------- molar masses [kg/mol] --------------------
Mi = [0.0319988 0.0280134 0.0159994 0.0140067 ...
      0.0159989 0.0140070 5.4858e-7];

% -------------------- SI -> cgs --------------------
rho_cgs   = rhogas_SI .* 1.0e-3;   % kg/m^3 -> g/cm^3
Rgas_cgs  = Rgas_SI   .* 1.0e4;    % J/(kg K) -> erg/(g K)
p_cgs     = rho_cgs .* Rgas_cgs .* T;           % dyn/cm^2
n_tot_cgs = p_cgs ./ (kB_cgs * T);             % cm^-3

% -------------------- mole fractions --------------------
x_all = [xO2, xN2, xO, xN, xOplus, xNplus, xe];
if sum(x_all) <= 0 || any(~isfinite(x_all))
    x_all = [0.2 0.8 0 0 0 0 0];
end
x_all = max(x_all,0);
x_all = x_all ./ sum(x_all);

% -------------------- Table V --------------------
[S0, SN2N, SNN, SNe, See_over_S0_lnL, SprN2N, SprNN] = Hansen_TableV(T); %#ok<ASGLU>

% -------------------- Coulomb log --------------------
xe_use  = max(x_all(7), 0.0);
n_e_cgs = max(xe_use * n_tot_cgs, 1.0e-80);

arg = (3.0./(2.0.*e_cgs.^3)) * sqrt((kB_cgs.^3 .* T.^3) ./ (pi_c .* n_e_cgs));
arg = max(arg, 1.0e-80);
ln_kappa_raw = log(arg);

% Hansen clamp (consistent with what you were doing)
ln_kappa = max(min(ln_kappa_raw, 12.0), 2.0);

% -------------------- diagnostics --------------------
dbg = struct();
dbg.T = T;
dbg.p_cgs = p_cgs;
dbg.n_tot_cgs = n_tot_cgs;
dbg.xe = xe_use;
dbg.n_e_cgs = n_e_cgs;
dbg.ln_kappa_raw = ln_kappa_raw;
dbg.ln_kappa = ln_kappa;

% -------------------- electron-electron term --------------------
% Table V provides See/(S0 lnΛ). Multiply by lnΛ to get See/S0.
See_over_S0 = See_over_S0_lnL .* ln_kappa;

% Safety: if Table-V electron cols are blank (low T), keep them OFF (0), not NaN
if ~isfinite(SNe),            SNe = 0;            end
if ~isfinite(See_over_S0),    See_over_S0 = 0;    end

% -------------------- build matrices --------------------
Sij = build_Sij(SN2N, SNN, SNe, See_over_S0);
Spr = build_Spr(SprN2N, SprNN, SNe, See_over_S0);

end


% ======================================================================
%  Hansen Table V data + interpolation (NaN-aware for cols 5–6)
% ======================================================================
function [S0, SN2N_over_S0, SNN_over_S0, SNe_over_S0, See_over_S0_lnL, ...
          SprN2N_over_S0, SprNN_over_S0] = Hansen_TableV(T)

% Columns:
% 1: T
% 2: S0
% 3: S(N2-N)/S0
% 4: S(N-N)/S0
% 5: S(N-e)/S0
% 6: S(e-e)/(S0 lnΛ)
% 7: S'(N2-N)/S0
% 8: S'(N-N)/S0

TableV = [
   500      38.4  0.946  0.894  NaN    NaN     0.877  0.761
  1000      34.9  0.920  0.838  NaN    NaN     0.843  0.703
  1500      33.7  0.889  0.785  NaN    NaN     0.817  0.652
  2000      33.2  0.886  0.742  NaN    NaN     0.794  0.611
  2500      32.8  0.846  0.705  NaN    NaN     0.775  0.578
  3000      32.6  0.830  0.676  NaN    NaN     0.759  0.551
  3500      32.4  0.815  0.650  NaN    NaN     0.745  0.527
  4000      32.3  0.803  0.628  NaN    NaN     0.733  0.507
  4500      32.2  0.792  0.608  NaN    NaN     0.722  0.489
  5000      32.1  0.782  0.591  NaN    NaN     0.712  0.473
  5500      32.0  0.773  0.575  0.397  89.90   0.703  0.458
  6000      32.0  0.764  0.561  0.380  75.60   0.695  0.445
  6500      31.9  0.757  0.548  0.366  64.50   0.688  0.433
  7000      31.9  0.750  0.536  0.353  55.70   0.681  0.422
  7500      31.9  0.743  0.524  0.342  48.60   0.674  0.412
  8000      31.8  0.737  0.514  0.331  42.80   0.668  0.402
  8500      31.8  0.731  0.504  0.321  37.90   0.662  0.393
  9000      31.8  0.725  0.495  0.330  33.80   0.657  0.385
  9500      31.8  0.720  0.486  0.304  30.40   0.652  0.377
 10000      31.8  0.715  0.478  0.297  27.40   0.647  0.370
 10500      31.7  0.710  0.470  0.290  24.90   0.642  0.363
 11000      31.7  0.706  0.463  0.283  22.70   0.637  0.356
 11500      31.7  0.701  0.456  0.281  20.80   0.633  0.350
 12000      31.7  0.697  0.448  0.270  19.00   0.629  0.342
 12500      31.7  0.693  0.443  0.266  17.60   0.625  0.338
 13000      31.7  0.689  0.437  0.261  16.27   0.621  0.332
 13500      31.7  0.689  0.431  0.256  15.10   0.618  0.327
 14000      31.7  0.689  0.426  0.252  14.00   0.616  0.322
 14500      31.6  0.689  0.420  0.247  13.09   0.613  0.316
 15000      31.6  0.689  0.415  0.243  12.24   0.610  0.312
];

T_tab = TableV(:,1);

S0               = interp1(T_tab, TableV(:,2), T, 'linear', 'extrap');
SN2N_over_S0     = interp1(T_tab, TableV(:,3), T, 'linear', 'extrap');
SNN_over_S0      = interp1(T_tab, TableV(:,4), T, 'linear', 'extrap');

% NaN-aware interpolation for electron columns (5–6)
interp_nan = @(x,y,xq) interp1(x(isfinite(y)), y(isfinite(y)), xq, 'linear', 'extrap');
SNe_over_S0      = interp_nan(T_tab, TableV(:,5), T);
See_over_S0_lnL  = interp_nan(T_tab, TableV(:,6), T);

SprN2N_over_S0   = interp1(T_tab, TableV(:,7), T, 'linear', 'extrap');
SprNN_over_S0    = interp1(T_tab, TableV(:,8), T, 'linear', 'extrap');

% Final safety
if ~isfinite(SN2N_over_S0),    SN2N_over_S0 = 0; end
if ~isfinite(SNN_over_S0),     SNN_over_S0  = 0; end
if ~isfinite(SNe_over_S0),     SNe_over_S0  = 0; end
if ~isfinite(See_over_S0_lnL), See_over_S0_lnL = 0; end
if ~isfinite(SprN2N_over_S0),  SprN2N_over_S0 = 0; end
if ~isfinite(SprNN_over_S0),   SprNN_over_S0  = 0; end

end


% ======================================================================
%  Build Sij and Spr matrices (dimensionless, grouped Hansen style)
% ======================================================================
function Sij = build_Sij(SN2N_over_S0, SNN_over_S0, SNe_over_S0, See_over_S0)
% UNPRIMED collision matrix Sij = S_ij / S0
%
% Species order:
%   1: O2   2: N2   3: O   4: N   5: O+   6: N+   7: e-

Sij = zeros(7,7);

% Molecule–molecule: S0
Sij(1,1) = 1.0;
Sij(2,2) = 1.0;
Sij(1,2) = 1.0;   Sij(2,1) = 1.0;

% Molecule–atom / molecule–ion: SN2N
sN2N = SN2N_over_S0;

Sij(2,3) = sN2N; Sij(3,2) = sN2N;
Sij(2,4) = sN2N; Sij(4,2) = sN2N;
Sij(2,5) = sN2N; Sij(5,2) = sN2N;
Sij(2,6) = sN2N; Sij(6,2) = sN2N;

Sij(1,3) = sN2N; Sij(3,1) = sN2N;
Sij(1,4) = sN2N; Sij(4,1) = sN2N;
Sij(1,5) = sN2N; Sij(5,1) = sN2N;
Sij(1,6) = sN2N; Sij(6,1) = sN2N;

% Atom–atom, atom–ion, ion–ion: SNN
sNN = SNN_over_S0;

Sij(3,3) = sNN;
Sij(4,4) = sNN;
Sij(3,4) = sNN; Sij(4,3) = sNN;

Sij(3,5) = sNN; Sij(5,3) = sNN;
Sij(3,6) = sNN; Sij(6,3) = sNN;
Sij(4,5) = sNN; Sij(5,4) = sNN;
Sij(4,6) = sNN; Sij(6,4) = sNN;

Sij(5,5) = sNN;
Sij(6,6) = sNN;
Sij(5,6) = sNN; Sij(6,5) = sNN;

% Electron collisions: SNe and See
for j = 1:6
    Sij(7,j) = SNe_over_S0;
    Sij(j,7) = SNe_over_S0;
end
Sij(7,7) = See_over_S0;

% Final NaN guard
Sij(~isfinite(Sij)) = 0;

end


function Spr = build_Spr(SprN2N_over_S0, SprNN_over_S0, SNe_over_S0, See_over_S0)
% PRIMED collision matrix Spr = S'_ij / S0
%
% Species order:
%   1: O2   2: N2   3: O   4: N   5: O+   6: N+   7: e-

Spr = zeros(7,7);

% Molecule–molecule (primed): S0
Spr(1,1) = 1.0;
Spr(2,2) = 1.0;
Spr(1,2) = 1.0; Spr(2,1) = 1.0;

% Molecule–atom / molecule–ion: SprN2N
sN2N = SprN2N_over_S0;

Spr(2,3) = sN2N; Spr(3,2) = sN2N;
Spr(2,4) = sN2N; Spr(4,2) = sN2N;
Spr(2,5) = sN2N; Spr(5,2) = sN2N;
Spr(2,6) = sN2N; Spr(6,2) = sN2N;

Spr(1,3) = sN2N; Spr(3,1) = sN2N;
Spr(1,4) = sN2N; Spr(4,1) = sN2N;
Spr(1,5) = sN2N; Spr(5,1) = sN2N;
Spr(1,6) = sN2N; Spr(6,1) = sN2N;

% Atom–atom, atom–ion, ion–ion: SprNN
sNN = SprNN_over_S0;

Spr(3,3) = sNN;
Spr(4,4) = sNN;
Spr(3,4) = sNN; Spr(4,3) = sNN;

Spr(3,5) = sNN; Spr(5,3) = sNN;
Spr(3,6) = sNN; Spr(6,3) = sNN;
Spr(4,5) = sNN; Spr(5,4) = sNN;
Spr(4,6) = sNN; Spr(6,4) = sNN;

Spr(5,5) = sNN;
Spr(6,6) = sNN;
Spr(5,6) = sNN; Spr(6,5) = sNN;

% Electrons: Hansen does NOT define primed electron sections
% Keep unprimed values for completeness only
for j = 1:6
    Spr(7,j) = SNe_over_S0;
    Spr(j,7) = SNe_over_S0;
end
Spr(7,7) = See_over_S0;

% Final NaN guard
Spr(~isfinite(Spr)) = 0;

end
