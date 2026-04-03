% file: Hansen_main.m
clc; clear; close all;
format long

% Symbolic calculations based on HANSEN's paper (APPROXIMATIONS FOR THE
% THERMODYNAMIC AND TRANSPORT PROPERTIES OF HIGH-TEMPERATURE AIR By
% C. FREDERICK HANSEN - NASA TR R-50 - 1959
%
% Temperature is in Kelvin and the pressure is in atmospheres

syms T p R h k Avo rhobar QN2 QO2 QO QN QOplus QNplus Qeminus
syms QpN2 QpO2 QpO QpN QpOplus QpNplus Qpeminus
syms QcN2 QcO2 QcO QcN QcOplus QcNplus Qceminus
syms lnKpOtoOpluseminus lnKpNtoNpluseminus 
syms DOion DNion 
syms ELHS HLHS E H E0
syms MN2 MO2 MO MN MOplus MNplus Meminus 
syms eps1  eps2 eps3  M0 x sumofx


% -------------------- PARTITION FUNCTIONS (Hansen) -----------------------
lnQN2 = 7/2 * log (T) - 0.42 - log( 1 - exp (-3390 / T ) ) - log (p);
lnQO2 = 7/2 * log (T) + 0.11 - log( 1 - exp ( -2270 / T ) ) + ...
        log ( 3 + 2 * exp ( -11390 / T ) + exp ( -18990 / T ) )  - log ( p );
lnQO  = 5/2 * log (T ) + 0.5 + ...
        log ( 5 + 3 * exp ( -228 / T) + exp (-326 / T ) + ...
              5 * exp ( -22800 / T ) + exp ( -48600 / T ) ) - log (p );
lnQN  = 5/2 * log ( T ) + 0.30 + ...
        log ( 4 + 10* exp ( -27700 / T ) + 6 * exp ( -41500 / T ) ) - log ( p );
lnQOplus  = 5/2 * log ( T ) + 0.50 + ...
            log ( 4 + 10*exp( -38600 / T ) + 6 * exp ( -58200 / T ) ) - log ( p ); 
lnQNplus  = 5/2 * log ( T ) + 0.30 + ...
            log ( 1 + 3 * exp ( -70.6 / T ) + 5 * exp ( -188.9 / T ) + ...
                  5 * exp ( -22000 / T ) + exp ( -47000 / T ) + ...
                  5 * exp ( -67900 / T ) ) - log ( p ); 
lnQeminus = 5/2 * log ( T ) - 14.24 - log ( p );

C1 = log(p) - log(R) - log(T);
lnQN2_c = lnQN2 + C1;
lnQO2_c = lnQO2 + C1;
lnQO_c  = lnQO + C1;
lnQN_c  = lnQN + C1;
lnQOplus_c  = lnQOplus + C1;
lnQNplus_c  = lnQNplus + C1;
lnQeminus_c = lnQeminus + C1;

C2 = log(p);
lnQN2_p = lnQN2 + C2;
lnQO2_p = lnQO2 + C2;
lnQO_p  = lnQO + C2;
lnQN_p  = lnQN + C2;
lnQOplus_p  = lnQOplus + C2;
lnQNplus_p  = lnQNplus + C2;
lnQeminus_p = lnQeminus + C2;

% -------------------- ENERGIES, ENTHALPIES, ENTROPIES --------------------
DN2   = 113200 * k;
DO2   = 59000  * k; 
DNplus= 168800 * k;
DOplus= 158000 * k;

E0N2 = 0; E0O2 = 0;
E0N  = DN2/2;  E0O  = DO2/2;
E0Nplus = E0N + DNplus;
E0Oplus = E0O + DOplus;
E0eminus = 0;

EN2     = (R*T*T*diff(lnQN2_c,T)     + E0N2*Avo);
EO2     = (R*T*T*diff(lnQO2_c,T)     + E0O2*Avo);
EO      = (R*T*T*diff(lnQO_c,T)      + E0O*Avo);
EN      = (R*T*T*diff(lnQN_c,T)      + E0N*Avo);
EOplus  = (R*T*T*diff(lnQOplus_c,T)  + E0Oplus*Avo);
ENplus  = (R*T*T*diff(lnQNplus_c,T)  + E0Nplus*Avo);
Eeminus = (R*T*T*diff(lnQeminus_c,T) + E0eminus*Avo);

HN2     = (R*T*T*diff(lnQN2_p,T)     + E0N2*Avo);
HO2     = (R*T*T*diff(lnQO2_p,T)     + E0O2*Avo);
HO      = (R*T*T*diff(lnQO_p,T)      + E0O*Avo);
HN      = (R*T*T*diff(lnQN_p,T)      + E0N*Avo);
HOplus  = (R*T*T*diff(lnQOplus_p,T)  + E0Oplus*Avo);
HNplus  = (R*T*T*diff(lnQNplus_p,T)  + E0Nplus*Avo);
Heminus = (R*T*T*diff(lnQeminus_p,T) + E0eminus*Avo);

SN2     = lnQN2_p     + (HN2    -E0N2*Avo)/(R*T);
SN      = lnQN_p      + (HN     -E0N*Avo )/(R*T);
SO2     = lnQO2_p     + (HO2    -E0O2*Avo)/(R*T);
SO      = lnQO_p      + (HO     -E0O*Avo )/(R*T);
SNplus  = lnQNplus_p  + (HNplus -E0Nplus*Avo)/(R*T);
SOplus  = lnQOplus_p  + (HOplus -E0Oplus*Avo)/(R*T);
Seminus = lnQeminus_p + (Heminus-E0eminus*Avo)/(R*T);

% ---------------------- EQUILIBRIUM CONSTANTS ----------------------------
lnKpO2to2O       = -59000   / T + 2*lnQO_p  - lnQO2_p;
lnKpN2to2N       = -113200  / T + 2*lnQN_p  - lnQN2_p;
lnKpOtoOplusemin = -158000  / T + lnQOplus_p + lnQeminus_p - lnQO_p;
lnKpNtoNplusemin = -168800  / T + lnQNplus_p + lnQeminus_p - lnQN_p;

% d(lnKp)/dT (kept, but not used inside Hansen_Transport)
dlnKp_O2to2O_sym = diff(lnKpO2to2O, T);
dlnKp_N2to2N_sym = diff(lnKpN2to2N, T);
dlnKp_Oion_sym   = diff(lnKpOtoOplusemin, T);
dlnKp_Nion_sym   = diff(lnKpNtoNplusemin, T);

Kp1   = exp(lnKpO2to2O);
Kp2   = exp(lnKpN2to2N);
KpNto = exp(lnKpNtoNplusemin);
KpOto = exp(lnKpOtoOplusemin);
Kp3   = 0.2 * KpOto + 0.8 * KpNto;

CON   = 1/(R*T);
Kc1   = Kp1 * CON;
Kc2   = Kp2 * CON;
KcNto = KpNto * CON;
KcOto = KpOto * CON;
Kc3   = 0.2*KcOto + 0.8*KcNto;

Z = 1 + eps1 + eps2 + 2*eps3;

xO2    = (0.2 - eps1)        / Z; 
xN2    = (0.8 - eps2)        / Z;
xO     = (2*eps1 - 0.4*eps3) / Z; 
xN     = (2*eps2 - 1.6*eps3) / Z;
xOplus = 0.4*eps3                / Z;
xNplus = 1.6*eps3                / Z;
xeminus= 2*eps3              / Z;

Energy   = xN2*EN2 + xO2*EO2 + xO*EO + xN*EN + ...
           xOplus*EOplus + xNplus*ENplus + xeminus*Eeminus;
Enthalpy = xN2*HN2 + xO2*HO2 + xO*HO + xN*HN + ...
           xOplus*HOplus + xNplus*HNplus + xeminus*Heminus;

deps1dT_p = diff(log(Kp1), T) / (2/eps1 - 1/(1+eps1) + 1/(0.2-eps1));
deps2dT_p = diff(log(Kp2), T) / (2/eps2 - 1/(1.2+eps2) + 1/(0.8-eps2));
deps3dT_p = diff(log(Kp3), T) / (2/eps3 - 1/(1+eps3) + 1/(1-eps3));

deps1dT_rho = diff(log(Kc1), T) / (2/eps1 + 1/(0.2-eps1));
deps2dT_rho = diff(log(Kc2), T) / (2/eps2 + 1/(0.8-eps2));
deps3dT_rho = diff(log(Kc3), T) / (2/eps3 + 1/(1-eps3));

dZdT_p   = deps1dT_p   + deps2dT_p   + 2*deps3dT_p;
dZdT_rho = deps1dT_rho + deps2dT_rho + 2*deps3dT_rho;

% Cvmix and Cpmix are Z*Cv/R and Z*Cp/R (Eqs. 51–52)
Cvmix =  Z * (  xN2*diff(EN2,T)/R + xO2*diff(EO2,T)/R + ...
                xO *diff(EO ,T)/R + xN *diff(EN ,T)/R + ...
                xOplus*diff(EOplus ,T)/R + ...
                xNplus*diff(ENplus ,T)/R + ...
                xeminus*diff(Eeminus,T)/R ) ...
      + T * ( ...
          (EN2/(R*T)) * ( dZdT_rho*xN2 + Z*( diff(xN2,eps1)*deps1dT_rho + diff(xN2,eps2)*deps2dT_rho + diff(xN2,eps3)*deps3dT_rho ) ) ...
        + (EO2/(R*T)) * ( dZdT_rho*xO2 + Z*( diff(xO2,eps1)*deps1dT_rho + diff(xO2,eps2)*deps2dT_rho + diff(xO2,eps3)*deps3dT_rho ) ) ...
        + (EO /(R*T)) * ( dZdT_rho*xO  + Z*( diff(xO ,eps1)*deps1dT_rho + diff(xO ,eps2)*deps2dT_rho + diff(xO ,eps3)*deps3dT_rho ) ) ...
        + (EN /(R*T)) * ( dZdT_rho*xN  + Z*( diff(xN ,eps1)*deps1dT_rho + diff(xN ,eps2)*deps2dT_rho + diff(xN ,eps3)*deps3dT_rho ) ) ...
        + (EOplus /(R*T)) * ( dZdT_rho*xOplus + Z*( diff(xOplus,eps1)*deps1dT_rho + diff(xOplus,eps2)*deps2dT_rho + diff(xOplus,eps3)*deps3dT_rho ) ) ...
        + (ENplus /(R*T)) * ( dZdT_rho*xNplus + Z*( diff(xNplus,eps1)*deps1dT_rho + diff(xNplus,eps2)*deps2dT_rho + diff(xNplus,eps3)*deps3dT_rho ) ) ...
        + (Eeminus/(R*T)) * ( dZdT_rho*xeminus+ Z*( diff(xeminus,eps1)*deps1dT_rho + diff(xeminus,eps2)*deps2dT_rho + diff(xeminus,eps3)*deps3dT_rho ) ) ...
      );

Cpmix =  Z * ( xN2*(diff(EN2,T)/R + 1) + ...
               xO2*(diff(EO2,T)/R + 1) + ...
               xO *(diff(EO ,T)/R + 1) + ...
               xN *(diff(EN ,T)/R + 1) + ...
               xOplus *(diff(EOplus ,T)/R + 1) + ...
               xNplus *(diff(ENplus ,T)/R + 1) + ...
               xeminus*(diff(Eeminus,T)/R + 1) ) ...
      + T * ( ...
        (EN2/(R*T)+1)*( dZdT_p*xN2 + Z*( diff(xN2,eps1)*deps1dT_p + diff(xN2,eps2)*deps2dT_p + diff(xN2,eps3)*deps3dT_p ) ) ...
      + (EO2/(R*T)+1)*( dZdT_p*xO2 + Z*( diff(xO2,eps1)*deps1dT_p + diff(xO2,eps2)*deps2dT_p + diff(xO2,eps3)*deps3dT_p ) ) ...
      + (EO /(R*T)+1)*( dZdT_p*xO  + Z*( diff(xO ,eps1)*deps1dT_p + diff(xO ,eps2)*deps2dT_p + diff(xO ,eps3)*deps3dT_p ) ) ...
      + (EN /(R*T)+1)*( dZdT_p*xN  + Z*( diff(xN ,eps1)*deps1dT_p + diff(xN ,eps2)*deps2dT_p + diff(xN ,eps3)*deps3dT_p ) ) ...
      + (EOplus /(R*T)+1)*( dZdT_p*xOplus + Z*( diff(xOplus,eps1)*deps1dT_p + diff(xOplus,eps2)*deps2dT_p + diff(xOplus,eps3)*deps3dT_p ) ) ...
      + (ENplus /(R*T)+1)*( dZdT_p*xNplus + Z*( diff(xNplus,eps1)*deps1dT_p + diff(xNplus,eps2)*deps2dT_p + diff(xNplus,eps3)*deps3dT_p ) ) ...
      + (Eeminus/(R*T)+1)*( dZdT_p*xeminus+ Z*( diff(xeminus,eps1)*deps1dT_p + diff(xeminus,eps2)*deps2dT_p + diff(xeminus,eps3)*deps3dT_p ) ) ...
      );

gamma = Cpmix / Cvmix;
Mmix  = xO2*MO2 + xN2*MN2 + xO*MO + xN*MN + ...
        xOplus*MOplus + xNplus*MNplus + xeminus*Meminus;
Rgas  = R / Mmix;
rhogas= p / (Rgas*T);

S = Z * ( xO2*SO2 + xN2*SN2 + xO*SO + xN*SN + ...
          xOplus*SOplus + xNplus*SNplus + xeminus*Seminus ...
        - ( xO2*log(xO2) + xN2*log(xN2) + xO*log(xO) + xN*log(xN) + ...
            xOplus*log(xOplus) + xNplus*log(xNplus) + xeminus*log(xeminus) ) ...
        - log(p) );

multiplier   = (1 + (T/Z)*dZdT_rho) / (1 + (T/Z)*dZdT_p);
a2_rho_o_p   = gamma * multiplier;
speedofsound = sqrt(a2_rho_o_p * (Z*R*T)/M0);




    CvO2_over_R_func    = diff(EO2    , T) / R ;
    CvN2_over_R_func    = diff(EN2    , T) / R ;
    CvO_over_R_func     = diff(EO     , T) / R ;
    CvN_over_R_func     = diff(EN     , T) / R ;
    CvOplus_over_R_func = diff(EOplus , T) / R ;
    CvNplus_over_R_func = diff(ENplus , T) / R ;
    Cve_over_R_func     = diff(Eeminus, T) / R ;





% ------------------- NUMERICAL LOOP --------------------------------------
ttt = 500 : 50 : 15000;
for i=1:length(ttt)
    clc
    fprintf("Loading: %f %%",i/length(ttt)*100)
    
    R = 8.314;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p = 1;                  % pressure [atm]%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    T = ttt(i);
    C1 = log(p) - log(R) - log(T);
    C2 = log(p);
    k  = 1.380649e-23;

    MO2     = 0.0319988;
    MN2     = 0.0280134;  
    MO      = 0.0159994; 
    MN      = 0.0140067; 
    MOplus  = 0.0159989; 
    MNplus  = 0.014007; 
    Meminus = 5.4858e-07;
    M0      = 0.029;
    Avo     = 6.02214076e23;
    h       = 6.62607015e-34;

    % eps1, eps2, eps3 from Hansen algebraic relations
    eps1_sym = ( -0.8 + sqrt(0.64 + 0.8*(1 + 4*p/Kp1)) ) / (2*(1 + 4*p/Kp1));
    eps2_sym = ( -0.4 + sqrt(0.16 + 3.84*(1 + 4*p/Kp2)) ) / (2*(1 + 4*p/Kp2));
    eps3_sym = ( 1 + p/Kp3 )^(-1/2);
    
    eps1 = double(subs(eps1_sym));
    eps2 = double(subs(eps2_sym));
    eps3 = double(subs(eps3_sym));
    
    % --- Mixture properties for transport
    Mmix_num = double(subs(Mmix));
    Rgas_num = double(subs(Rgas));

    h_mass_SI_num = double(subs(Enthalpy)) / Mmix_num;   % (J/mol) / (kg/mol) = J/kg
    h_mass_SI(i)  = h_mass_SI_num;

    p_SI      = p * 101325;
    rhogas_SI = p_SI / (Rgas_num * T);

    cp_mass_SI_num = double(subs( (Cpmix / Z) * (R / Mmix) ));  % J/(kg·K)

    xO2_num    = double(subs(xO2));
    xN2_num    = double(subs(xN2));
    xO_num     = double(subs(xO ));
    xN_num     = double(subs(xN ));
    xOplus_num = double(subs(xOplus));
    xNplus_num = double(subs(xNplus));
    xe_num     = double(subs(xeminus));
    
    xO2plot(i)     = xO2_num;
    xN2plot(i)     = xN2_num;
    xOplot(i)      = xO_num;
    xNplot(i)      = xN_num;
    xOplusplot(i)  = xOplus_num;
    xNplusplot(i)  = xNplus_num;
    xeminusplot(i) = xe_num;
    
    dlnKp_O2to2O_num = double(subs(dlnKp_O2to2O_sym));
    dlnKp_N2to2N_num = double(subs(dlnKp_N2to2N_sym));
    dlnKp_Oion_num   = double(subs(dlnKp_Oion_sym));
    dlnKp_Nion_num   = double(subs(dlnKp_Nion_sym));

    % --- Transport coefficients (Hansen Table V)
[Mi, x_all_out, Sij, Spr, See_over_S0] = Hansen_Transport( ...
    T, M0, Mmix_num, rhogas_SI, Rgas_num, cp_mass_SI_num, ...
    xO2_num, xN2_num, xO_num, xN_num, xOplus_num, xNplus_num, xe_num, ...
    dlnKp_O2to2O_num, dlnKp_N2to2N_num, dlnKp_Oion_num, dlnKp_Nion_num);

C_over_R_vec = [
    double(subs( CvO2_over_R_func ));
    double(subs( CvN2_over_R_func ));
    double(subs( CvO_over_R_func ));
    double(subs( CvN_over_R_func ));
    double(subs( CvOplus_over_R_func ));
    double(subs( CvNplus_over_R_func ));
    double(subs( Cve_over_R_func ))
];


% ================= HANSEN VISCOSITY ===========================
% η/η0  (dimensionless)
eta_raw = Hansen_Viscosity(M0, Mi, x_all_out, Sij);     % η/η0
eta_ratio(i) = eta_raw;                                 % store ratio

% Hansen reference viscosity η0(T)   [poise]
eta0_cgs = 1.462e-5 * sqrt(T) / (1 + 112/T);            % Hansen Eq. (68)
eta0_cgs_array(i) = eta0_cgs;

% Absolute viscosity μ(T) [poise, SI]
mu_dim_cgs(i) = eta_raw * eta0_cgs;                     % g/(cm·s)
mu_dim_SI(i)  = 0.1 * mu_dim_cgs(i);                    % Pa·s


% ================= HANSEN CONDUCTIVITY =========================

[k_ratio_raw, ~, k_n_ratio(i), k_r_ratio(i)] = Hansen_Conductivity( ...
    T, R, M0, Mi, x_all_out, Sij, Spr, C_over_R_vec, ...
    dlnKp_O2to2O_num, dlnKp_N2to2N_num, dlnKp_Oion_num, dlnKp_Nion_num);

k_ratio(i) = k_ratio_raw;                               % k/k0 for plotting

% Hansen reference conductivity k0 (Eq. 73)
eta0_SI = 0.1 * eta0_cgs;                               % Pa·s
k0_SI   = (19/4) * (R/M0) * eta0_SI;                    % W/(m·K)

% Dimensional conductivity
k_dim(i) = k_ratio_raw * k0_SI;                         % W/(m·K)


% ================= PRANDTL NUMBER ==============================
Pr(i) = mu_dim_SI(i) * cp_mass_SI_num / k_dim(i);



    % --- Store thermodynamic quantities
    temp(i)             = T;
    sound(i)            = double(vpa(subs(a2_rho_o_p), 6));
    ratio(i)            = double(vpa(subs(gamma),       6));
    cpplot(i)           = double(vpa(subs(Cpmix),       6));
    cvplot(i)           = double(vpa(subs(Cvmix),       6));
    Mmixplot(i)         = double(vpa(subs(Mmix),        6));
    Zplot(i)            = double(vpa(subs(Z),           6));
    Energyplot(i)       = double(vpa(subs(Energy*Z/R/T),6));
    eps1plot(i)         = eps1;
    eps2plot(i)         = eps2;
    eps3plot(i)         = eps3;
    Splot(i)            = double(vpa(subs(S),           6));
    multiplierplot(i)   = double(vpa(subs(multiplier),  6));
    Kc1plot(i)          = double(vpa(subs(Kc1),         6));
    Kc2plot(i)          = double(vpa(subs(Kc2),         6));
    Kc3plot(i)          = double(vpa(subs(Kc3),         6));
    speedofsoundplot(i) = double(vpa(subs(speedofsound),6));



cp_mass_SI_num      = double(subs( (Cpmix / Z) * (R / Mmix) ));
cp_mass_SI(i)       = cp_mass_SI_num;          % <-- new line
Pr(i)               = mu_dim_SI(i) * cp_mass_SI_num / k_dim(i);


end



k0_SI_vec    = (19/4) * (R/M0) .* (0.1 .* eta0_cgs_array);
k_from_Pr    = mu_dim_SI .* cp_mass_SI ./ Pr;
k_ratio_check = k_from_Pr ./ k0_SI_vec;

rel_err = abs((k_ratio_check - k_ratio) ./ k_ratio);
fprintf('Max relative error in conductivity consistency = %.3e\n', max(rel_err));



% =====================================================================
% EXPORT VARIABLES FOR AUTOMATED COMPARISON (used by run_all.m)
% =====================================================================

% Temperature array
T_vals = temp(:);

% Cp_vals: use cp_mass_SI array (actual mass-based Cp)
Cp_vals = cp_mass_SI(:);


% Compressibility factor
Z_vals = Zplot(:);

% Dimensional viscosity [Pa·s]
mu_vals = mu_dim_SI(:);

% Dimensional conductivity [W/(m·K)]
k_vals = k_dim(:);

% Prandtl number
Pr_vals = Pr(:);

% =====================================================================
fprintf('\nHansen_main: Export variables prepared for run_all.m\n');





% =============================
% High-Readability Journal Style
% =============================
set(0,'defaultFigureColor','white');

% Fonts
set(0,'defaultAxesFontName','Times');
set(0,'defaultAxesFontSize',26);      % axis numbers
set(0,'defaultTextFontSize',26);      % axis labels, titles

% Lines
set(0,'defaultLineLineWidth',2.0);    % primary line width
set(0,'defaultAxesLineWidth',2.0);    % thickness of axes

% Legend appearance
set(0,'defaultLegendFontSize',20);
set(0,'defaultLegendLocation','best');

% Figure size (good for 2-column scaling)
set(0,'defaultFigureUnits','inches');
set(0,'defaultFigurePosition',[1 1 6 4.5]);  % width x height


% -------------------- Plots: titles + axis labels + legends (B/W only) --------------

% --- Figure 1: a^2 rho/p and gamma ---
figure(1); hold on; grid on;
plot(temp, sound, '-k',  'LineWidth',1.3);    % solid
plot(temp, ratio, '--k', 'LineWidth',1.3);    % dashed
xlabel('Temperature (K)');
ylabel('a^2 \rho/p   and   \gamma');
title('a^2 \rho/p (solid) and \gamma (dashed)');
legend('a^2 \rho/p','\gamma','Location','Best');
ylim([1 1.7]);

% --- Figure 2: Cv and Cp ---
figure(2); hold on; grid on;
plot(temp, cvplot, '-k',  'LineWidth',1.3);   % solid
plot(temp, cpplot, '--k', 'LineWidth',1.3);   % dashed
xlabel('Temperature (K)');
ylabel('Z C_v / R   and   Z C_p / R');
title('Specific Heats (Hansen Eqs. 51–52)');
legend('Z C_v / R','Z C_p / R','Location','Best');

% --- Figure 3: Mmix ---
figure(3); hold on; grid on;
plot(temp, Mmixplot, '-k','LineWidth',1.3);
xlabel('Temperature (K)');
ylabel('Molecular Weight (kg/mol)');
title('Mixture Molecular Weight');

% --- Figure 5: Z ---
figure(5); hold on; grid on;
plot(temp, Zplot, '-k','LineWidth',1.3);
xlabel('Temperature (K)');
ylabel('Z');
title('Compressibility Factor (Z)');

% --- Figure 6: Energy ---
figure(6); hold on; grid on;
plot(temp, Energyplot, '-k','LineWidth',1.3);
xlabel('Temperature (K)');
ylabel('E Z / (R T)');
title('Non-dimensional Internal Energy');

% --- Figure 7: eps1, eps2, eps3 ---
figure(7); hold on; grid on;
plot(temp, eps1plot, '-k',  'LineWidth',1.3);   % solid
plot(temp, eps2plot, '--k', 'LineWidth',1.3);   % dashed
plot(temp, eps3plot, '-.k', 'LineWidth',1.3);   % dash-dot
xlabel('Temperature (K)');
ylabel('\epsilon parameters');
legend('\epsilon_1','\epsilon_2','\epsilon_3','Location','Best');
title('Dissociation / Ionization Parameters');

% --- Figure 8: Entropy ---
figure(8); hold on; grid on;
plot(temp, Splot, '-k','LineWidth',1.3);
xlabel('Temperature (K)');
ylabel('Z S / R');
title('Mixture Entropy (Hansen Eq. 50)');

% --- Figure 9: multiplier ---
figure(9); hold on; grid on;
plot(temp, multiplierplot, '-k','LineWidth',1.3);
xlabel('Temperature (K)');
ylabel('Multiplier');
title('Thermochemical Multiplier to \gamma');

% --- Figure 10: Speed of sound ---
figure(10); hold on; grid on;
plot(temp, speedofsoundplot, '-k','LineWidth',1.3);
xlabel('Temperature (K)');
ylabel('Speed of sound (m/s)');
title('Equilibrium Speed of Sound');


% --- Figure 11: Species mole fractions ---
figure(11); hold on; grid on; box on;

% neutral species
plot(temp, xO2plot,   '-k',  'LineWidth',2);    % O2
plot(temp, xN2plot,   '--k', 'LineWidth',2);    % N2
plot(temp, xOplot,    '-.k', 'LineWidth',2);    % O
plot(temp, xNplot,    ':k',  'LineWidth',2);    % N

% ions / electrons (thicker)
plot(temp, xOplusplot,  '-k',  'LineWidth',3.5);  % O+
plot(temp, xNplusplot,  '--k', 'LineWidth',3.5);  % N+
plot(temp, xeminusplot, '-.k', 'LineWidth',3.5);  % e-

xlabel('Temperature (K)');
ylabel('Mole Fraction');
title('Equilibrium Mole Fractions — Hansen Model');

ylim([0 1]);      % 0–1% mole fraction
box on;              % explicit (journal-safe)

legend('O_2','N_2','O','N','O^+','N^+','e^-','Location','Best');

% ============================================================
%   Figure 12 — Viscosity Ratio μ/μ0
% ============================================================
figure(12); clf; hold on; grid on;
plot(temp, eta_ratio, '-k', 'LineWidth', 1.6);
xlabel('Temperature (K)');
ylabel('\mu / \mu_0');
title('Hansen Viscosity Ratio');
set(gca,'FontSize',14);


% ============================================================
%   Figure 13 — Conductivity Ratio k/k0
% ============================================================
figure(13); clf; hold on; grid on;
plot(temp, k_ratio, '-k', 'LineWidth', 1.6);
xlabel('Temperature (K)');
ylabel('k / k_0');
title('Hansen Conductivity Ratio');
set(gca,'FontSize',14);


% ============================================================
%   Figure 14 — Prandtl Number
% ============================================================
figure(14); clf; hold on; grid on;
plot(temp, Pr, '-k', 'LineWidth', 1.8);
xlabel('Temperature (K)');
ylabel('Pr');
title('Hansen Prandtl Number');
ylim([0 2]);
set(gca,'FontSize',14);


% ============================================================
%   Figure 15 — Dimensional Viscosities μ(T) and μ0(T)
% ============================================================
figure(15); clf; hold on; grid on;

% μ dimensional (poise or convert to SI if preferred)
plot(temp, mu_dim_cgs, '-k', 'LineWidth', 1.8);

% μ0 reference viscosity
plot(temp, eta0_cgs_array, '--k', 'LineWidth', 1.8);

xlabel('Temperature (K)');
ylabel('Viscosity (poise)');
title('\mu(T) and \mu_0(T)');
legend('\mu(T)', '\mu_0(T)', 'Location', 'Best');
set(gca,'FontSize',14);


% ============================================================
%   Journal-Style Export of Selected Figures (Vector PDF)
% ============================================================
baseDir = fileparts(mfilename('fullpath'));
if isempty(baseDir), baseDir = pwd; end
outDir  = fullfile(baseDir,'fig_export');
if ~exist(outDir,'dir'), mkdir(outDir); end

figList = [1 2 3 5 6 7 8 9 10 11 12 13 14 15];

FigureNames = {
    'Hansen_Fig01_a2rho_over_p_and_gamma.pdf'
    'Hansen_Fig02_Cp_Cv_normalized.pdf'
    'Hansen_Fig03_Mixture_Molecular_Weight.pdf'
    'Hansen_Fig05_Compressibility_Z.pdf'
    'Hansen_Fig06_Internal_Energy_ND.pdf'
    'Hansen_Fig07_Epsilon_Parameters.pdf'
    'Hansen_Fig08_Entropy_ZS_over_R.pdf'
    'Hansen_Fig09_Thermochemical_Multiplier.pdf'
    'Hansen_Fig10_Speed_of_Sound.pdf'
    'Hansen_Fig11_Equilibrium_Mole_Fractions.pdf'
    'Hansen_Fig12_Viscosity_Ratio.pdf'
    'Hansen_Fig13_Conductivity_Ratio.pdf'
    'Hansen_Fig14_Prandtl_Number.pdf'
    'Hansen_Fig15_Dimensional_Viscosities.pdf'
};

for ii = 1:numel(figList)
    fnum = figList(ii);

    if ~ishandle(fnum), continue; end
    figure(fnum);

    title(''); % if you want titles stripped
    fname = fullfile(outDir, FigureNames{ii});

    exportgraphics(gcf, fname, ...
        'ContentType','vector', ...
        'BackgroundColor','white');

    fprintf('Saved FIGURE %d  →  %s\n', fnum, fname);
end


%% ===============================================================
%   SAVE HANSEN RESULTS TO MAT FILE (FINAL, DATA-ONLY)
% ===============================================================

Hansen = struct();

Hansen.T  = temp(:);
Hansen.h  = h_mass_SI(:);
Hansen.Cp = cp_mass_SI(:);
Hansen.Z  = Zplot(:);
Hansen.mu = mu_dim_SI(:);
Hansen.k  = k_dim(:);
Hansen.Pr = Pr(:);

% -----------------------------------------------------------
% Additional Hansen quantities needed for comparison
% -----------------------------------------------------------

Hansen.gamma = ratio(:);
Hansen.a2rho_over_p = sound(:);

Hansen.CvZ_over_R = cvplot(:);
Hansen.CpZ_over_R = cpplot(:);

Hansen.Mmix = Mmixplot(:);
Hansen.EnergyND = Energyplot(:);
Hansen.speed = speedofsoundplot(:);

% --------------------------
% Species mole fractions
% --------------------------
Hansen.X = struct();
Hansen.X.O2    = xO2plot(:);
Hansen.X.N2    = xN2plot(:);
Hansen.X.O     = xOplot(:);
Hansen.X.N     = xNplot(:);
Hansen.X.Oplus = xOplusplot(:);
Hansen.X.Nplus = xNplusplot(:);
Hansen.X.e     = xeminusplot(:);

outfile = fullfile(baseDir, sprintf('Hansen_p_%g.mat', p));
save(outfile,'Hansen');
fprintf('Saved Hansen results → %s\n', outfile);



%% =====================================================================
%  DIMENSIONAL SUMMARY PLOTS (END-OF-FILE ONLY, NO MIXING OF CALCS)
%  - Build dimensional E, H, S, Cv, Cp, R, mu, k, a on the existing T grid
%  - Plot as compact 2x3 "one-page" figure
% =====================================================================

% ---- Basic arrays (force column) ----
T_dim   = temp(:);                 % K
p_atm   = p;                       % (you used p=1 in loop)
p_SI    = p_atm * 101325;          % Pa

h_dim   = h_mass_SI(:);            % J/kg  (already computed)
Cp_dim  = cp_mass_SI(:);           % J/(kg K) (already computed)
mu_dim  = mu_dim_SI(:);            % Pa*s (already computed)
k_dim2  = k_dim(:);                % W/(m K) (already computed)
gam_dim = ratio(:);                % gamma (already computed)
a_dim   = speedofsoundplot(:);     % m/s (already computed)

Mmix_dim = Mmixplot(:);            % kg/mol (already computed)
Z_dim    = Zplot(:);               % compressibility function (dimensionless)

% ---- Reconstruct mixture gas constant Rmix(T) [J/(kg K)] ----
% Use Cp = gamma/(gamma-1) * R  => R = Cp*(gamma-1)/gamma
Rmix_dim = Cp_dim .* (gam_dim - 1) ./ gam_dim;

% ---- Cv(T) consistent with Cp and gamma ----
Cv_dim = Cp_dim ./ gam_dim;        % since gamma = Cp/Cv

% ---- Density and p/rho = R*T (for ideal gas mixture form) ----
rho_dim = p_SI ./ (Rmix_dim .* T_dim);   % kg/m^3
p_over_rho = p_SI ./ rho_dim;            % J/kg  (= R*T)

% ---- Internal energy e(T) from h = e + p/rho ----
E_dim = h_dim - p_over_rho;        % J/kg

% ---- Dimensional entropy s(T) from Hansen S (dimensionless):
% Hansen uses S such that dimensional s = (R/Mmix)*S, but your loop uses:
%   R is universal gas constant and Mmix is mixture molar mass.
% Therefore Rmix = R/Mmix. We'll use reconstructed Rmix_dim above.
if exist('Splot','var') && ~isempty(Splot)
    S_Hansen_dimless = Splot(:);          % Hansen "S" (dimensionless)
    s_dim = Rmix_dim .* S_Hansen_dimless; % J/(kg K)
else
    s_dim = NaN(size(T_dim));
end

% ---- Quick sanity checks (optional prints) ----
fprintf('\n[Dimensional summary] Rmix range = [%.2f, %.2f] J/(kg K)\n', ...
    min(Rmix_dim), max(Rmix_dim));
fprintf('[Dimensional summary] rho range  = [%.3e, %.3e] kg/m^3 (at p=%.2f atm)\n', ...
    min(rho_dim), max(rho_dim), p_atm);

%% =====================================================================
%  ONE-PAGE FIGURE (2x3) — dimensional quantities
% =====================================================================

% Smaller figure for "one page" insert; adjust if you want tighter
figDim = figure(200); clf;
set(figDim,'Units','inches','Position',[1 1 7.5 5.5],'Color','w');

tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

% (1) Cp and Cv
nexttile; hold on; grid on; box on;
plot(T_dim, Cp_dim, '-k','LineWidth',1.6);
plot(T_dim, Cv_dim, '--k','LineWidth',1.6);
xlabel('T (K)'); ylabel('J/(kg·K)');
legend('C_p','C_v','Location','best');
title(''); % no title

% (2) Rmix
nexttile; hold on; grid on; box on;
plot(T_dim, Rmix_dim, '-k','LineWidth',1.6);
xlabel('T (K)'); ylabel('R_{mix} (J/(kg·K))');
title('');

% (3) h and e
nexttile; hold on; grid on; box on;
plot(T_dim, h_dim, '-k','LineWidth',1.6);
plot(T_dim, E_dim, '--k','LineWidth',1.6);
xlabel('T (K)'); ylabel('J/kg');
legend('h','e','Location','best');
title('');

% (4) s (if available)
nexttile; hold on; grid on; box on;
plot(T_dim, s_dim, '-k','LineWidth',1.6);
xlabel('T (K)'); ylabel('s (J/(kg·K))');
title('');

% (5) mu and k (dual axis)
nexttile; grid on; box on;
yyaxis left
plot(T_dim, mu_dim, '-k','LineWidth',1.6);
ylabel('\mu (Pa·s)');
yyaxis right
plot(T_dim, k_dim2, '--k','LineWidth',1.6);
ylabel('k (W/(m·K))');
xlabel('T (K)');
title('');
legend('\mu','k','Location','best');

% (6) gamma and speed of sound
nexttile; grid on; box on;
yyaxis left
plot(T_dim, gam_dim, '-k','LineWidth',1.6);
ylabel('\gamma');
yyaxis right
plot(T_dim, a_dim, '--k','LineWidth',1.6);
ylabel('a (m/s)');
xlabel('T (K)');
title('');
legend('\gamma','a','Location','best');

%% =====================================================================
% OPTIONAL: Export this one-page summary as a PDF into your figs folder
% (Uncomment and pick a filename that mainHansen.tex will include if desired)
% =====================================================================
% outDirSummary = fullfile(baseDir,'figs');   % or 'fig_export'
% if ~exist(outDirSummary,'dir'), mkdir(outDirSummary); end
% exportgraphics(figDim, fullfile(outDirSummary,'figXX_dimensional_summary.pdf'), ...
%     'ContentType','vector','BackgroundColor','white');
% fprintf('Saved dimensional summary page -> %s\n', fullfile(outDirSummary,'figXX_dimensional_summary.pdf'));
