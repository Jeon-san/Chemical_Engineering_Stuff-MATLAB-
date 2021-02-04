%% Unit Ops Assignment 2 - Distillation column
%% Name: Jeon San 
%% Mat no.: U1820704G

%% Binary Distillation solver 
clc
clear

zf = input('Please enter feed Ethanol concentration [g EtOH/g mixture]clc:*try around 0.35');%0.35; % Feed EtOH mass percentage
feed_mol = input('Please enter feed flowrate [kmol/hr] *try something around 100-1000:'); %500; % Molar flowrate in [kmol/hr]

mm_etoh = 46.07; % Molar mass of EtOH [g/mol]
mm_h2o = 18.015; % Molar mass of EtOH [g/mol]

%% Solve for Feed Mass Flow Rate,FMFR (Used in heat balance later)

%%%SYSTEM 1
syms etoh_mol h2o_mol

eqn1 = feed_mol == etoh_mol + h2o_mol; % Overall MB
eqn2 = zf == etoh_mol*mm_etoh/(etoh_mol*mm_etoh + h2o_mol*mm_h2o); % EtOH MB

eqns1 = [eqn1,eqn2];
sol = solve(eqns1,[h2o_mol etoh_mol],'ReturnConditions',true); % Solve eqns

etoh_mol = double(sol.etoh_mol); % Mol of etoh in feed
h2o_mol = double(sol.h2o_mol); % Mol of h2o in feed

FMFR = (etoh_mol*mm_etoh + h2o_mol*mm_h2o)/1000; % Feed mass flow rate
%% Target concentration calculation
% 95% Recovery of EtOH and H2O

etoh_D = 0.95*etoh_mol;
etoh_B = 0.05*etoh_mol;
h2o_B = 0.95*h2o_mol;
h2o_D = 0.05*h2o_mol;

xd = etoh_D/(etoh_D + h2o_D); % EtOH Mol percentage of distillate 
xb = etoh_B/(etoh_B + h2o_B); % EtOH Mol percentage of bottom product

%% Weight Percentage (Used in heat balance later)
xd_weight = xd*(mm_etoh)/(xd*(mm_etoh)+(1-xd)*(mm_h2o)); % EtOH Weight percentage of distillate
xb_weight = xb*(mm_etoh)/(xb*(mm_etoh)+(1-xb)*(mm_h2o)); % EtOH Weight percentage of bottom product

%% Solve for D and B mass flow rate (Used in heat balance later)

%%% SYSTEM 2
syms D B

eqn3 = FMFR == D+B;
eqn4 = FMFR*zf == D*xd_weight + B*xb_weight;

eqns2 = [eqn3,eqn4];
sol2 = solve(eqns2,[D B],'ReturnConditions',true); % Solve eqns

D = double(sol2.D); % Mass flow rate of distillate
B = double(sol2.B); % Mass flow rate of bottom product

%% Solve for D and B molar flow rate 
%%% SYSTEM 3

syms D_mol B_mol

eqn5 = feed_mol == D_mol + B_mol;
eqn6 = etoh_mol == D_mol*xd + B_mol*xb;

eqns3 = [eqn5,eqn6];
sol3 = solve(eqns3,[D_mol B_mol],'ReturnConditions',true); % Solve eqns

D_mol = double(sol3.D_mol); % Mass flow rate of distillate
B_mol = double(sol3.B_mol); % Mass flow rate of bottom product


%% Calculate q using EB and feed EtOH mol fraction

H_feed = 40; % Enthalpy of feed
q = (510-H_feed)/(510-66); % Enthalphy values read from H_xy diagram
zf_mol = etoh_mol/feed_mol; % EtOH molar fraction in feed

%% Guess R (reflux) 

R = 1.4 %% Initial guess of R <<<<<<<<<<<<<<<<<<<<<<<<<
imp=0;
h = waitbar(imp, 'Distillation column started...');  
error = 1; % set error to initiate iteration
counter = 0;
while abs(error)>0.001
   
    %%%%%%% Solve for rectification and feed op line intersection
    %% Solve for Boil-up ratio 

    %%% SYSTEM 5
    syms V V_bar vb lo

    eqn7 = R == lo/D_mol;
    eqn8 = V == lo + D_mol;
    eqn9 = vb == V_bar/B_mol;
    eqn10 = V == V_bar - feed_mol*(q-1);

    eqns4 = [eqn7 eqn8 eqn9 eqn10];
    sol4 = vpasolve(eqns4,[V V_bar vb lo]);

    V = double(sol4.V);
    V_bar = double(sol4.V_bar);
    vb = double(sol4.vb);
    lo = double(sol4.lo);
    syms y_f x_f

    %% Find intersection
    % Rectifying op line
    rec = y_f == (R/(R+1))*x_f + (1/(R+1))*xd;
    % Feed op line
    feed = y_f == (q/(q-1))*x_f - zf_mol/(q-1);
    
    sol_f = vpasolve([rec feed],[y_f x_f]);
    y_f = double(sol_f.y_f);
    x_f = double(sol_f.x_f); % IMPORTANT: change from rectification to stripping here

   

    % vapour from top tray = distillate composition

    y = xd; % Vapour at tray 1 is has mole frac equal to distillate due to total condenser

    %% Going down the column
    
    tray_no = 17; % Number of trays
    x_track = []; 
    T_track = [];
    y_track = [];
    for j = 1:(tray_no-1) % 16 Trays, therefore 16-1 = 15
        
        % Clean track variables to prevent accumulation
        
        
        [x T] = NRTL(y); % x1 and T using NRTL model
        x_track(j) = x; % Record down x
        T_track(j) = T; % Record down T
        if x<x_f
            y = ((vb+1)/vb)*x-(1/vb)*xb; % y1 of next tray (Stripping op line)
            y_track(j) = y; % Record down y
            opline_track(j) = 1;
        else 
            y = (R/(R+1))*x + (1/(R+1))*xd; % y1 of next tray (Rectifying op line)
            y_track(j) = y; % Record down y
            opline_track(j) = 0;
        end
    [x_final T_final] = NRTL(y);   %16th Tray 
    x_track = [x_track x_final]; 
    T_track = [T_track T_final];
    y_track = [xd y_track]; 
    end


b_recovery = ((1-x_final)*B_mol)/h2o_mol;

counter = counter+1; % number of times passed through loop
error = (0.95-b_recovery)/0.95;
error_array(counter) = error;
error_ini = error_array(1);
error_fin = error_array(end);
imp = abs((error_ini-error_fin)/error_ini);
waitbar(imp,h,'Hi Prof. Poernomo! Please wait......');

con_sped = 2;
R = R + R*error *con_sped
end
close(h)

%% With minimum reflux ratio, find boil-up

R_true =R; % Assign reflux to new 'True Reflux ratio'


%% Solve for Boil-up ratio 

%%% SYSTEM 6
syms V V_bar vb lo

eqn11 = R_true == lo/D_mol;
eqn12 = V == lo + D_mol;
eqn13 = vb == V_bar/B_mol;
eqn14 = V == V_bar - feed_mol*(q-1);

eqns6 = [eqn7 eqn8 eqn9 eqn10];
sol6 = vpasolve(eqns6,[V V_bar vb lo]);

V_true = double(sol4.V);
V_bar_true = double(sol4.V_bar);
vb_true = double(sol4.vb);
lo_true = double(sol4.lo);

fprintf('--------------------------\n\n')
fprintf('The reflux and boil-up ratios are %.3f and %.3f',R_true,vb_true);
fprintf('\n\n--------------------------')

%% Find ideal feed entry location

%%%%%%% Solve for rectification and feed op line intersection

syms y_f_true x_f_true

% Rectifying op line
rec = y_f_true == (R_true/(R_true+1))*x_f_true + (1/(R_true+1))*xd;
% Feed op line
feed = y_f_true == (q/(q-1))*x_f_true - zf_mol/(q-1);

sol_f = vpasolve([rec feed],[y_f_true x_f_true]);
y_f_true = double(sol_f.y_f_true);
x_f_true = double(sol_f.x_f_true); % IMPORTANT: change from rectification to stripping here


hold on

figure(1)
plot(1:tray_no,x_track,'-s');
yline(x_f_true);

tray = ['Tray ';'Tray ';'Tray ';'Tray ';'Tray ';'Tray ';'Tray ';'Tray ';'Tray ';'Tray ';'Tray ';'Tray ';'Tray ';'Tray ';'Tray ';'Tray ';'Cond.';];
a = [1:17]'; b = strcat(tray,num2str(a)); c = cellstr(b);
dx = 0.03; dy = 0.03; % displacement so the text does not overlay the data points
text([1:tray_no]+dx, x_track+dy, c,'FontSize',7);
ylabel('EtOH molar fraction');
xlabel('Tray #  (left to right is condenser to boiler)');

hold off


fprintf('\n\n')
fprintf('The ideal Tray location is tray 15 as its the first point after the intersection of operating lines');
fprintf('\n\n--------------------------')


%% Overall Heat balance to get Cond. and Reboiler duty

%% Assume total condensation


H_d_vap = 280; % read off from graph [kcal/kg]
H_d_liq = 50; % read off from graph [kcal/kg]
H_evap = H_d_vap-H_d_liq; % Enthapy of evaporation

V_mass = ((V_true*xd*mm_etoh)+(V_true*(1-xd)*mm_h2o))/1000; %kg/hr
Qc = V_mass*H_evap; % Condenser duty [kcal/hr]

% Overall energy balance to get boiler duty

H_d = H_d_liq; % Enthalpy of sat. liq. distillate
H_b = 95; % Enthalpy of bottom prod from graph


syms Qb

heat = FMFR*H_feed -Qc + Qb == D*H_d + B*H_b;

Qb = vpasolve(heat,Qb);

fprintf('\n\n')
fprintf('The condenser duty is %.4f kcal/hr',Qc);
fprintf('\n\n')
fprintf('The boiler duty is %.4f kcal/hr',Qb);
fprintf('\n\n--------------------------')


hold on

figure('Renderer', 'painters', 'Position', [10 10 900 600])
plot(1:tray_no,T_track,'-s');

tray = ['Tray ';'Tray ';'Tray ';'Tray ';'Tray ';'Tray ';'Tray ';'Tray ';'Tray ';'Tray ';'Tray ';'Tray ';'Tray ';'Tray ';'Tray ';'Tray ';'Cond.';];
a = [1:17]'; b = strcat(tray,num2str(a)); c = cellstr(b);
dx = 0.03; dy = 1; % displacement so the text does not overlay the data points
text([1:tray_no]+dx, T_track+dy, c,'FontSize',7);
pbaspect([2.5 1 1])
title('Temperature profile of distillation column')
ylabel('Temperature(k)');
xlabel('Tray #  (left to right is condenser to boiler)');


hold off