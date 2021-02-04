function [x1 T] = NRTL(y1)

%% This function returns y1, y2, and T for a given x using NRTL

%% NRTL Thermodynamic model by Jeon San

%% Key parameters
P = 100; % Pressure in column in kPa (1bar)
R = 8.314; % gas constant [J/mol]


y1 = y1; % Range of x1 values (mol% of etoh) to solve for other vars
y2 = 1-y1 ; % mol percent h2o

syms T x1 x2 

A_12 = -633; % J/mol
A_21 = 5823.1; % J/mol

tau_12 = A_12/(R*T); % Tau parameter
tau_21 = A_21/(R*T); % Tau parameter

G_12 = exp(-0.3*tau_12); % G12 parameter
G_21 = exp(-0.3*tau_21); % G21 parameter

% Expression for activity coefficients gamma1 and gamma2 in parts

%%% ETOH
A1 = tau_21*(G_21/(x1+x2*G_21))^2; 
B1 = (G_12*tau_12)/((x2+x1*G_12)^2);
gamma_1 = exp((x2^2)*(A1+B1));
%%% H2O
A2 = tau_12*(G_12/(x2+x1*G_12))^2;
B2 = (G_21*tau_21)/((x1+x2*G_21)^2);
gamma_2 = exp((x1^2)*(A2+B2));

Psat_etoh = exp(16.8958-3795.17/(T-42.232)); % vapour pressure of EtOH (kPa)
Psat_h2o = exp(16.3872-3885.7/(T-42.980)); % vapour pressure of EtOH (kPa)

% Equations to be solved
eqn1 = P*y1 == gamma_1*x1*Psat_etoh;
eqn2 = P*y2 == gamma_2*x2*Psat_h2o;
eqn3 = 1 == x1 + x2;

%% Use Raoults law to create initial guess for higher convergence

r1 = P*y1 == x1*Psat_etoh; % Raoult's law for EtOH
r2 = P*y2 == x2*Psat_h2o; % Raoult's law for H2O
r3 = 1 == x1+ x2;

guess_sol = vpasolve([r1 r2 r3],[T x1 x2]);

T_ini = double(guess_sol.T);
x1_ini = double(guess_sol.x1);
x2_ini = double(guess_sol.x2);

eqns = [eqn1 eqn2 eqn3];
sol = vpasolve(eqns,[T x1 x2],[T_ini,x1_ini,x2_ini]);

% Save solved T,y1, and y2 in arrays
T = double(sol.T);
x1 = double(sol.x1);
x2 = double(sol.x2);

end

