function MDFN_3c_Charge_edited
%% Code by Eloise Tredenick 2023, University of Oxford
%% The following paper must be cited when using or modifying this code - Tredenick et al 2024 J. Electrochem. Soc. https://doi.org/10.1149/1945-7111/ad5767

%% 2026 update, Eloise Tredenick along with University of Canberra Capstone IT unit, 
% improved run time (without plotting) from 65 sec to 36 sec. 
% Group members: Karma Tobgyel, Kinley Wangchuk, Sonam Tobgay, Naveed Khan, Sangay Choden
%and Tenzin Lhaden. New code produces identical solution, improved by
%changing for loop in RHS to vectorised version, putting indecies outside
%RHS and changing tols.

%Half cell MDFN model. p1=NMC, p2=LFP. Same code used for discharge and
%other C rates D_e and Kappa_e functions, D_s constant

% Fast version
N_r = 55;
N_x = 35;

V_0_cs0 = 3; %Initial voltage for charge
I_sign = 1;  %Sign of current. -1 when discharging, +1 when charging

%% If discharging instead
%V_0_cs0 =  4.193; % Initial voltage for discharge
%I_sign = -1;%in discharge
%v_min = 2.5;
%plus change end event function

%% Parameters for NMC - p1
eps_e_p1 = 0.31; %Porosity
eps_CBD_p1 = 0.11; %CBD volume fraction
brugg_p1 = 1.6 ;%Tortuosity
Rs_p1 = 4.942e-6; %Radius of solid particles in m
sigma_p1 = 5;   %Electronic conductivity in solid in S/m
Ds_p1_bulk = 3.9952e-14; %Diffusivity of Li ions in solid in m^2/s
k_p1_ref = 1e-10;%Reaction rate constant in m^{2.5}s^{-1}mol^{-0.5}

%% Parameters for LFP - p2
eps_e_p2 = 0.2625; %Porosity
eps_CBD_p2 = 0.11; %CBD volume fraction
brugg_p2 = 2.1;%Tortuosity
Rs_p2 = 4.283e-7; %Radius of solid particles in m
sigma_p2 =  5;  %Electronic conductivity in solid in S/m
Ds_p2_bulk = 3e-16;%Diffusivity of Li ions in solid in m^2/s
k_p2_ref = 8e-13;%Reaction rate constant in m^{2.5}s^{-1}mol^{-0.5}

%% Parameters for Separator
brugg_s = 1.5; %Tortuosity
eps_e_s = 0.45;%porosity

Rc = 9.5;  %Contact resistance [ohm.m^2]
capacity005C = 165; %Cell capacity in mAh/g at 0.05C
cell_weightNMCLFP = 0.035061; %in grams
mAh = 5.761875; %mAh of cell to convert to applied current
Crate = 3; %3C charge
current_A = mAh * Crate *1e-3; %Applied current in Amps
Area_electrode	= 1.54e-4; %Surface area of electrode in m^2
i_app =  current_A/Area_electrode;
v_max = 4.2;% max voltage cut-off
N = N_x;
Ns = N_x - 2; %cs doesn't need to be evaluated at boundaries of x
N2 = Ns * N_r; %Number of nodes for particle radius is the radius nodes * x nodes
t0 = 0;
tf = 12*60*60;  %t fianl in seconds - keep high for 0.1C
tspan1 = [t0 tf];

%% Typical values - used for non-dimensionalisation to scale order of magnitudes
D_char_e        = 5.34e-10;     % Typical De in m^2/s
Kappa_char      = 1;            %Typical kappa in S/m
c_char          = 1e3;          %Typical value of lithium-ion concentration in electrolyte in mol/m3
ce0             = c_char;

%% NMC OCP
load('NMC_OCP1_fit.mat','U_Dischargep1'); %fit based on exp data

%% LFP OCP
U_Dischargep2 = @(y) 3.413025424326315  + 0.001./y + 0.001./(y-1) ;     %is mean of LFP data in middle. y=cs/csmax
%NMC - now find what cs/csmax should be based on initial voltage
high_NMC = 0.923;
xxx=linspace(0,high_NMC,100000);
U1_check = U_Dischargep1(xxx); %Now rearrange to find y in U_Dischargep1
[U_I,~] = find(U1_check<=V_0_cs0,1,'first');
cs0_on_csmax_p1 = xxx(U_I);
%LFP
lhs = 1000 * (V_0_cs0 - 3.413025424326315 ); %Now rearrange to find y in U_Dischargep2
cs0_on_csmax_p2 =   (- sqrt (lhs ^2 +4) + lhs + 2 )/(2*lhs);

%% CONSTANTS AND TEMPERATURE
R = 8.31446261815324;  %[Pa m^3/ (K mol)]
Fara = 96485.33212;  %Faraday constant [C.mol-1]
Temp_c = 20; %Room temperature
T_ref = Temp_c + 273.15;           % Reference temperature [K]

tplus = 0.37; % Lithium transference number.

Ls = 16e-6;       % Thickness of separator, m
Lp1 = 44e-6 ;         % Thickness of cathode, m
Lp2 = 44e-6 ;         % Thickness of cathode, m

%% INITIL CONDITIONS & CONCENTRATIONS
cs_p1_max =  48700;
cs_p10 = cs_p1_max*cs0_on_csmax_p1;

cs_p2_max =  22806;
cs_p20 = cs_p2_max * cs0_on_csmax_p2;

a_p1 = 3* (1 -  eps_e_p1- eps_CBD_p1) / Rs_p1;
a_p2 = 3* (1 -  eps_e_p2- eps_CBD_p2) / Rs_p2;

L =  Ls + Lp1 + Lp2; 	% Thickness of the cell,m
L1 = Ls;
L2 = Lp1 + Ls;

%% De and kappa e functions
p1 =   -7.55e-21 ;
p2 =   8.231e-17 ;
p3 =  -3.401e-13 ;
p4 =   5.464e-10 ;
De_bulk_fun = @(cebar) p1.*(cebar.*ce0).^3 + p2.*(cebar.*ce0).^2 + p3.*(cebar.*ce0) + p4 ;

q1 =   5.746e-18  ;
q2 =  -1.025e-13  ;
q3 =   7.082e-10 ;
q4 =  -2.236e-06  ;
q5 =     0.00273 ;
q6 =   -0.003002  ;
Ke_bulk_fun =  @(cebar) q1.*(cebar .* ce0).^5 + q2.*(cebar .* ce0).^4 + q3.*(cebar .* ce0).^3 + q4.*(cebar .* ce0).^2 + q5.*(cebar .* ce0) + q6;
%% Mesh - ions
nodes_s = linspace(0, L1/L, N_x);
dx_s = nodes_s(2) - nodes_s(1);
faces_s = ([nodes_s(1) nodes_s] + [nodes_s nodes_s(end)]) ./ 2;
Dx_s = diff(faces_s(:));

nodes_p1 = linspace(  L1/L  ,  L2/L , N_x);
dx_p1 = nodes_p1(2) - nodes_p1(1);
faces_p1 = ([nodes_p1(1) nodes_p1] + [nodes_p1 nodes_p1(end)]) ./ 2;
Dx_p1 = diff(faces_p1(:));

nodes_p2 = linspace(L2/L, L/L, N_x);
dx_p2 = nodes_p2(2) - nodes_p2(1);
faces_p2 = ([nodes_p2(1) nodes_p2] + [nodes_p2 nodes_p2(end)]) ./ 2;
Dx_p2 = diff(faces_p2(:));

%% Mesh - r
nodesRp1 = linspace(0 , 1 , N_r);
drp1 = nodesRp1(2) - nodesRp1(1);
nodesRp2 = linspace(0 , 1 , N_r);
drp2 = nodesRp2(2) - nodesRp2(1);

%DIM VERSION
nodes_sepdim = nodes_s .* L;
nodes_p1dim = nodes_p1 .* L;
nodes_p2dim = nodes_p2 .* L;

nodes_dim = [nodes_sepdim nodes_p1dim nodes_p2dim];

%% Dimensionless Groupings and constants
taud = L^2 / D_char_e;
tspan_nodim1 = tspan1 ./ taud;
Mm1 = (1 - tplus) * L * i_app / (D_char_e * Fara * eps_e_p1 * c_char) ;
Mm2 = (1 - tplus) * L * i_app / (D_char_e * Fara * eps_e_p2 * c_char) ;
FonRT = Fara / (R * T_ref);
U1 = 2* L * a_p1 * k_p1_ref * c_char^0.5 * cs_p1_max *Fara/ i_app;
U2 = 2* L * a_p2 * k_p2_ref * c_char^0.5 * cs_p2_max *Fara/ i_app;
%cs paras
B1 = taud * Ds_p1_bulk / (Rs_p1 ^2);
B2 = taud * Ds_p2_bulk / (Rs_p2 ^2);
Ff1 = i_app * Rs_p1 / (L * Fara * a_p1 * Ds_p1_bulk * cs_p1_max);
Ff2 = i_app * Rs_p2 / (L * Fara * a_p2 * Ds_p2_bulk * cs_p2_max);
%phis paras
H1 = i_app * Fara * L/ (sigma_p1 * R * T_ref );
H2 = i_app * Fara * L/ (sigma_p2 * R * T_ref );
%phie paras
P = 2 * (1 - tplus);
%multilayer
sigp1 = sigma_p2 / sigma_p1;
%new ke de function values
kap_p1 = eps_e_p1 ^ brugg_p1 /Kappa_char ;
De_paras_p1 = eps_e_p1 ^ brugg_p1/D_char_e ;
Des_paras = eps_e_s ^ brugg_s /D_char_e ;
kap_s = eps_e_s ^ brugg_s /Kappa_char;
kap_p2 = eps_e_p2 ^ brugg_p2 /Kappa_char ;
De_paras_p2 = eps_e_p2 ^ brugg_p2/D_char_e ;

Ss = i_app * L * Fara / (R * T_ref * Kappa_char);
S_LP = (1 - tplus) * i_app * L / ( Fara * D_char_e * c_char );

%% Dimensionless initial conditions
c_bar_s_p10 =  cs_p10 ./ cs_p1_max .* ones(1, N2);
c_bar_s_p_surf100 = c_bar_s_p10(end);
c_bar_s_p20 =  cs_p20 ./ cs_p2_max .* ones(1, N2);
c_bar_s_p_surf200 = c_bar_s_p20(end);
c_bar_e_p10 = ce0./c_char.* ones(1, N_x);
c_bar_e_p20 = ce0./c_char.* ones(1, N_x);
c_bar_e_s0 = ce0./c_char.* ones(1, N_x);
phi_bar_ep10 = zeros(1, N_x);
phi_bar_ep20 = zeros(1, N_x);
phi_bar_es0 = zeros(1, N_x);

U_p1_bar0 =  FonRT .* U_Dischargep1(c_bar_s_p_surf100);
U_p2_bar0 =  FonRT .* U_Dischargep2(c_bar_s_p_surf200);

phi_bar_sp10 = U_p1_bar0 .* ones(1, N_x);
phi_bar_sp20 = U_p2_bar0 .* ones(1, N_x);
Jbar_p10 =  U1 .* sqrt(c_bar_e_p10(2:N-1) .* c_bar_s_p_surf100 .* (1 - c_bar_s_p_surf100)) ...
    .* sinh( 0.5 .* (phi_bar_sp10(2:N-1) - phi_bar_ep10(2:N-1)  - U_p1_bar0));
Jbar_p20 =  U2 .* sqrt(c_bar_e_p20(2:N-1) .* c_bar_s_p_surf200 .* (1 - c_bar_s_p_surf200)) ...
    .* sinh( 0.5 .* (phi_bar_sp20(2:N-1) - phi_bar_ep20(2:N-1)  - U_p2_bar0));

%% ICs
ICs_0_orig = [c_bar_s_p10 c_bar_s_p20 phi_bar_sp10 phi_bar_sp20 phi_bar_ep10 phi_bar_ep20 phi_bar_es0 c_bar_e_p10 c_bar_e_p20 c_bar_e_s0  Jbar_p10 Jbar_p20];

%Make M the mass matrix for ODE15s
%for cs
D_pat_xrd = eye(N2, N2);
for ii = 1:Ns
    D_pat_xrd(ii*N_r - (N_r-1),ii*N_r - (N_r-1)) = 0;
    D_pat_xrd(N_r*ii,N_r*ii) = 0;
end
%for ce
D_pat_x1            = eye(N, N);
D_pat_x1(1,1)       = 0;
D_pat_x1(end,end)   = 0;

icL = length(ICs_0_orig);
M   = sparse(icL,icL);

M(1:N2,1:N2) = D_pat_xrd;
M(N2+1:2*N2,N2+1:2*N2) = D_pat_xrd;
M( 2 * N2 + 5 * N  + 1 : 2 * N2 + 6 * N,  2 * N2 + 5 * N  + 1 : 2 * N2 + 6 * N) = D_pat_x1;
M(2 * N2 + 6 * N  + 1 : 2 * N2 + 7 * N, 2 * N2 + 6 * N  + 1 : 2 * N2 + 7 * N ) = D_pat_x1;
M(2 * N2 + 7 * N  + 1 : 2 * N2 + 8 * N ,  2 * N2 + 7 * N  + 1 : 2 * N2 + 8 * N  ) = D_pat_x1;

%% Precompute indices for state vector
idx_csp1   = 1:N2;
idx_csp2   = N2+1:2*N2;
idx_phisp1 = 2*N2 + 1 : 2*N2 + N;
idx_phisp2 = 2*N2 + N + 1 : 2*N2 + 2*N;
idx_phiep1 = 2*N2 + 2*N + 1 : 2*N2 + 3*N;
idx_phiep2 = 2*N2 + 3*N + 1 : 2*N2 + 4*N;
idx_phies  = 2*N2 + 4*N + 1 : 2*N2 + 5*N;
idx_cep1   = 2*N2 + 5*N + 1 : 2*N2 + 6*N;
idx_cep2   = 2*N2 + 6*N + 1 : 2*N2 + 7*N;
idx_ces    = 2*N2 + 7*N + 1 : 2*N2 + 8*N;
idx_Jp1    = 2*N2 + 8*N + 1 : 2*N2 + 8*N + Ns;
idx_Jp2    = 2*N2 + 8*N + Ns + 1 : 2*N2 + 8*N + 2*Ns;

%% Precompute radial terms
rp1_sq_outer = (nodesRp1(2:N_r)').^2;
rp1_sq_inner = (nodesRp1(2:N_r-1)').^2;
rp2_sq_outer = (nodesRp2(2:N_r)').^2;
rp2_sq_inner = (nodesRp2(2:N_r-1)').^2;
inv_eps_e_s = 1 / eps_e_s;
inv_eps_e_p1 = 1 / eps_e_p1;
inv_eps_e_p2 = 1 / eps_e_p2;
inv_dx_p1 = 1 / dx_p1;
inv_dx_p2 = 1 / dx_p2;
inv_dx_s = 1 / dx_s;
inv_Dx_p1_mid = 1 ./ Dx_p1(2:N-1);
inv_Dx_p2_mid = 1 ./ Dx_p2(2:N-1);
inv_Dx_s_mid = 1 ./ Dx_s(2:N-1);

%% Solve with ODE15s, for rhs at bottom
reltol = 1e-3;
abstol = 1e-5;
opts = odeset('RelTol', reltol, ...
    'AbsTol', abstol, ...
    'Mass', M, ...
    'Events', @(t,phi_sol) DFN_event(t,phi_sol));
[t_sol,phi_solution] = ode15s(@rhs,tspan_nodim1, ICs_0_orig, opts);

%% Dimensionalise the parameters
time_steps = length(t_sol(1:end,1));
num_plots = 15;
every = round(time_steps / num_plots);

fontsizeL = 12;
fontsize = 14;
linewidth = 3;
markersize = 10;

%Extract the solution
csp1_nondim      = phi_solution(:, idx_csp1);
csp2_nondim      = phi_solution(:, idx_csp2);
% phiphisp1_nondim = phi_solution(:, idx_phisp1);
phiphisp2_nondim = phi_solution(:, idx_phisp2);
phiep1_nondim    = phi_solution(:, idx_phiep1);
phiep2_nondim    = phi_solution(:, idx_phiep2);
phies_nondim     = phi_solution(:, idx_phies);
cep1_nondim      = phi_solution(:, idx_cep1);
cep2_nondim      = phi_solution(:, idx_cep2);
ces_nondim       = phi_solution(:, idx_ces);
Jbarp1           = phi_solution(:, idx_Jp1);
Jbarp2           = phi_solution(:, idx_Jp2);

%Dimensional versions
Csp_dim1 = cs_p1_max .* csp1_nondim;
phiep_dim1 =    phiep1_nondim./ FonRT;
phies_dim =   phies_nondim./ FonRT;
J_cathode_plot1 = Jbarp1 .* (i_app / (L * a_p1));
cep_dim1 = c_char  .*  cep1_nondim;
ces_dim = c_char  .*  ces_nondim;
% phiphisp_dim1 =  (R * T_ref / Fara).* phiphisp1_nondim;
Csp_dim2 = cs_p2_max .* csp2_nondim;
phiep_dim2 =    phiep2_nondim./ FonRT;
J_cathode_plot2 = Jbarp2 .* (i_app / (L * a_p2));
cep_dim2 = c_char  .*  cep2_nondim;
phiphisp_dim2 =  (R * T_ref / Fara).* phiphisp2_nondim;

phie_0_real = [phi_bar_es0 phi_bar_ep10 phi_bar_ep20 ]./ FonRT;

%Calc cs
surf_idx   = N_r * (1:Ns);
% centre_idx = N_r * (1:Ns) - (N_r - 1);

csp_surf1   = Csp_dim1(:, surf_idx);
% csp_centre1 = Csp_dim1(:, centre_idx);
csp_surf2   = Csp_dim2(:, surf_idx);
% csp_centre2 = Csp_dim2(:, centre_idx);

t_sec = taud .* t_sol;
V_cell = phiphisp_dim2(:,end) + Rc*current_A* I_sign;

cap_mAhg = current_A .* t_sec * 1e3./ ( 3600 * cell_weightNMCLFP );
cap_norm = cap_mAhg ./capacity005C;

%% Plotting
colmat = parula(6);
figure;
subplot(1,2,1);
set(gcf,'DefaultAxesColorOrder',colmat)
plot(cap_norm, V_cell, 'LineWidth', linewidth,'color', colmat(1,:));hold on;
set(xlabel('Capacity (-)'), 'FontSize', fontsize);
set(ylabel('Voltage (V)'), 'FontSize', fontsize);
set(gca,'FontSize',fontsize,'FontName','Times New Roman');

subplot(1,2,2);
set(gcf,'DefaultAxesColorOrder',colmat)
plot(cap_mAhg, V_cell, 'LineWidth', linewidth,'color', colmat(1,:));hold on;
set(xlabel('Capacity (mAh/g)'), 'FontSize', fontsize);
set(ylabel('Voltage (V)'), 'FontSize', fontsize);
set(gca,'FontSize',fontsize,'FontName','Times New Roman');

% ce
colmat = parula(num_plots+7);
figure;
subplot(2,1,1);
set(gcf,'DefaultAxesColorOrder',colmat)
plot113 = plot(nodes_dim,[ ces_dim(1,1:end)'; cep_dim1(1,1:end)'; cep_dim2(1,1:end)']  , '-kx', 'MarkerSize', markersize,'LineWidth', linewidth);hold on;
plot(nodes_dim(1:every:end,:), [ces_dim(1:every:end,:)'; cep_dim1(1:every:end,:)'; cep_dim2(1:every:end,:)'], 'LineWidth', linewidth);
plot443 = plot(nodes_dim,[ ces_dim(end,1:end)'; cep_dim1(end,1:end)'; cep_dim2(end,1:end)'], 'LineWidth', linewidth);
xline(0,'--k',{'LMP'},'LabelHorizontalAlignment','left');xline(0,'--k',{'s'});xline(L1,'--k',{'p1'});xline(L2,'--k',{'p2'});xline(L,'--k',{'CC'});
set(xlabel('x - m (s| p1| p2)'), 'FontSize', fontsize);
set(ylabel('c_e - mol/m^3'), 'FontSize', 6);
set(legend([ plot113 plot443],{'IC','End time soln'}, 'Location','best'), 'FontSize',6);
set(gca, 'FontSize', fontsize);

subplot(2,3,4);
set(gcf,'DefaultAxesColorOrder',colmat)
plot(t_sec, ces_dim(:,1)','LineWidth', linewidth);hold on;
plot(t_sec, ces_dim(:,end)', '-b', 'MarkerSize', markersize,'LineWidth', linewidth);
set(xlabel('time - sec'), 'FontSize', fontsize);
set(ylabel('c_{e,s} - mol/m^3'), 'FontSize', 6);
set(legend('L','R', 'Location','best'), 'FontSize',8);
set(gca, 'FontSize', fontsize);

subplot(2,3,5);
set(gcf,'DefaultAxesColorOrder',colmat)
plot(t_sec, cep_dim1(:,1)','LineWidth', linewidth);hold on;
plot(t_sec, cep_dim1(:,end)', '-b', 'MarkerSize', markersize,'LineWidth', linewidth);
set(xlabel('time - sec'), 'FontSize', fontsize);
set(ylabel('c_{e,p1} - mol/m^3'), 'FontSize', 6);
set(legend('Sep','p1|p2', 'Location','best'), 'FontSize',8);
set(gca, 'FontSize', fontsize);

subplot(2,3,6);
set(gcf,'DefaultAxesColorOrder',colmat)
plot(t_sec, cep_dim2(:,1)','LineWidth', linewidth);hold on;
plot(t_sec, cep_dim2(:,end)', '-b', 'MarkerSize', markersize,'LineWidth', linewidth);
set(xlabel('time - sec'), 'FontSize', fontsize);
set(ylabel('c_{e,p2} - mol/m^3'), 'FontSize', 6);
set(legend('p1|p2','CC', 'Location','best'), 'FontSize',8);
set(gca, 'FontSize', fontsize);

% phie
colmat = parula(num_plots+7);
figure;
subplot(2,1,1);
set(gcf,'DefaultAxesColorOrder',colmat)
plot1133 = plot(nodes_dim,phie_0_real', '-kx', 'MarkerSize', markersize,'LineWidth', linewidth);hold on;
plot(nodes_dim(1:every:end,:), [ phies_dim(1:every:end,:)'; phiep_dim1(1:every:end,:)'; phiep_dim2(1:every:end,:)'], 'LineWidth', linewidth);
plot4435 = plot(nodes_dim,[ phies_dim(end,1:end)'; phiep_dim1(end,1:end)'; phiep_dim2(end,1:end)'], 'LineWidth', linewidth);
xline(0,'--k',{'LMP'},'LabelHorizontalAlignment','left');xline(0,'--k',{'s'});xline(L1,'--k',{'p1'});xline(L2,'--k',{'p2'});xline(L,'--k',{'CC'});
set(xlabel('x - m ( s| p1| p2)'), 'FontSize', fontsize);
set(ylabel('\phi_e - V'), 'FontSize', 6);
set(legend([ plot1133 plot4435],{'IC','End time soln'}, 'Location','best'), 'FontSize',6);
set(gca, 'FontSize', fontsize);

subplot(2,3,4);
set(gcf,'DefaultAxesColorOrder',colmat)
plot(t_sec, phies_dim(:,1)','LineWidth', linewidth);hold on;
plot(t_sec, phies_dim(:,end)', '-b', 'MarkerSize', markersize,'LineWidth', linewidth);
set(xlabel('time - sec'), 'FontSize', fontsize);
set(ylabel('\phi_{e,s} - V'), 'FontSize', 6);
set(legend('L','R', 'Location','best'), 'FontSize',8);
set(gca, 'FontSize', fontsize);

subplot(2,3,5);
set(gcf,'DefaultAxesColorOrder',colmat)
plot(t_sec, phiep_dim1(:,1)', 'LineWidth', linewidth);hold on;
plot(t_sec, phiep_dim1(:,end)', '-b', 'MarkerSize', markersize,'LineWidth', linewidth);
set(xlabel('time - sec'), 'FontSize', fontsize);
set(ylabel('\phi_{e,p1} - V'), 'FontSize', 6);
set(legend('Sep','p1|p2', 'Location','best'), 'FontSize',8);
set(gca, 'FontSize', fontsize);

subplot(2,3,6);
set(gcf,'DefaultAxesColorOrder',colmat)
plot(t_sec, phiep_dim2(:,1)','LineWidth', linewidth);hold on;
plot(t_sec, phiep_dim2(:,end)', '-b', 'MarkerSize', markersize,'LineWidth', linewidth);
set(xlabel('time - sec'), 'FontSize', fontsize);
set(ylabel('\phi_{e,p2} - V'), 'FontSize', 6);
set(legend('p1|p2','CC', 'Location','best'), 'FontSize',8);
set(gca, 'FontSize', fontsize);

%% J and cs
figure;
%J
colmat = parula(num_plots+7);
subplot(2,2,1);
set(gcf,'DefaultAxesColorOrder',colmat)
plot(nodes_p1dim(2:N-1)./1e-6, J_cathode_plot1(1,:)', '-kx', 'MarkerSize', markersize,'LineWidth', linewidth);hold on;
plot(nodes_p1dim(2:N-1)./1e-6, J_cathode_plot1(1:every:end,:)', 'LineWidth', linewidth);
plot(nodes_p1dim(2:N-1)./1e-6, J_cathode_plot1(end,:)', 'LineWidth', linewidth);
set(xlabel('x_{p1} - \mum'), 'FontSize', fontsize,'FontName','Times New Roman');
set(ylabel('J_{p1} - A/m^2'), 'FontSize', fontsize,'FontName','Times New Roman');
set(gca,'FontSize',fontsize,'FontName','Times New Roman');

subplot(2,2,2);
set(gcf,'DefaultAxesColorOrder',colmat)
plot(nodes_p2dim(2:N-1)./1e-6, J_cathode_plot2(1,:)', '-kx', 'MarkerSize', markersize,'LineWidth', linewidth);hold on;
plot(nodes_p2dim(2:N-1)./1e-6, J_cathode_plot2(1:every:end,:)', 'LineWidth', linewidth);
plot(nodes_p2dim(2:N-1)./1e-6, J_cathode_plot2(end,:)', 'LineWidth', linewidth);
set(xlabel('x_{p2} - \mum'), 'FontSize', fontsize,'FontName','Times New Roman');
set(ylabel('J_{p2} - A/m^2'), 'FontSize', fontsize,'FontName','Times New Roman');
set(gca,'FontSize',fontsize,'FontName','Times New Roman');

%cs
colmat = parula(num_plots+7);
subplot(2,2,3);
set(gcf,'DefaultAxesColorOrder',colmat)
plot(nodes_p1dim(2:N-1)./1e-6, csp_surf1(1,:)./cs_p1_max', '-kx', 'MarkerSize', markersize,'LineWidth', linewidth);hold on;
plot(nodes_p1dim(2:N-1)./1e-6, csp_surf1(1:every:end,:)./cs_p1_max', 'LineWidth', linewidth);
plot(nodes_p1dim(2:N-1)./1e-6, csp_surf1(end,:)./cs_p1_max', 'LineWidth', linewidth);
set(xlabel('x_{p1} - \mum'), 'FontSize', fontsize,'FontName','Times New Roman');
set(ylabel('c_{s,p1,surf}/c_{s,p1}^{max}'), 'FontSize', fontsize,'FontName','Times New Roman');
set(gca,'FontSize',fontsize,'FontName','Times New Roman');

subplot(2,2,4);
set(gcf,'DefaultAxesColorOrder',colmat)
plot113=plot(nodes_p2dim(2:N-1)./1e-6, csp_surf2(1,:)./cs_p2_max', '-kx', 'MarkerSize', markersize,'LineWidth', linewidth);hold on;
plot(nodes_p2dim(2:N-1)./1e-6, csp_surf2(1:every:end,:)./cs_p2_max', 'LineWidth', linewidth);
plot443=plot(nodes_p2dim(2:N-1)./1e-6, csp_surf2(end,:)./cs_p2_max', 'LineWidth', linewidth);
set(xlabel('x_{p2} - \mum'), 'FontSize', fontsize,'FontName','Times New Roman');
set(ylabel('c_{s,p2,surf}/c^{max}_{s,p2}'), 'FontSize', fontsize,'FontName','Times New Roman');
set(legend([ plot113 plot443],{'Initial condition','Final solution'}, 'Location','best'), 'FontSize',fontsizeL,'FontName','Times New Roman');
set(gca,'FontSize',fontsize,'FontName','Times New Roman');




%% Plot cs in r
nodesdimRp1 = nodesRp1 .* Rs_p1;
%csp
colmat = parula(num_plots+5);
figure;
subplot(2,2,1);
set(gcf,'DefaultAxesColorOrder',colmat)
plot(nodesdimRp1./1e-6, Csp_dim1(1,N2-N_r+1:N2)'./cs_p1_max, '-kx', 'MarkerSize', markersize,'LineWidth', linewidth);hold on;
plot(nodesdimRp1./1e-6, Csp_dim1(1:every:end,N2-N_r+1:N2)./cs_p1_max, 'LineWidth', linewidth);
plot(nodesdimRp1./1e-6, Csp_dim1(end,N2-N_r+1:N2)./cs_p1_max , 'LineWidth', linewidth);
set(xlabel('\it r\rm_{p1} (m)'), 'FontSize', fontsize);
set(ylabel('\it c\rm_{s,p1}/ \it c^{max}_{s,p1} (mol/m^3)'), 'FontSize', 6);
set(title('Separator'), 'FontSize', 6);
set(gca, 'FontSize', fontsize);

% cs at x=0+1
subplot(2,2,2);
set(gcf,'DefaultAxesColorOrder',colmat)
plot111 = plot(nodesdimRp1./1e-6, Csp_dim1(1,1:N_r)'./cs_p1_max, '-kx', 'MarkerSize', markersize,'LineWidth', linewidth);hold on;
plot(nodesdimRp1./1e-6, Csp_dim1(1:every:end,1:N_r)./cs_p1_max, 'LineWidth', linewidth);
plot441 = plot(nodesdimRp1./1e-6, Csp_dim1(end,1:N_r)./cs_p1_max , 'LineWidth', linewidth);
set(xlabel('\it r\rm_{p1} (m)'), 'FontSize', fontsize);
set(ylabel('\it c\rm_{s,p1}/ \it c^{max}_{s,p1} (mol/m^3)'), 'FontSize', 6);
set(title('p1|p2 adjacent'), 'FontSize', 6);
set(legend([plot111 plot441],{'Initial Condition','Final Solution'}, 'Location','best'), 'FontSize',8);
set(gca, 'FontSize', fontsize);

nodesdimRp2 = nodesRp2 .* Rs_p2;
% Plot cs in r - done
%csp
colmat = parula(num_plots+5);
%figure;
subplot(2,2,3);
set(gcf,'DefaultAxesColorOrder',colmat)
plot(nodesdimRp2./1e-6, Csp_dim2(1,N2-N_r+1:N2)'./cs_p2_max, '-kx', 'MarkerSize', markersize,'LineWidth', linewidth);hold on;
plot(nodesdimRp2./1e-6, Csp_dim2(1:every:end,N2-N_r+1:N2)./cs_p2_max, 'LineWidth', linewidth);
plot(nodesdimRp2./1e-6, Csp_dim2(end,N2-N_r+1:N2) ./cs_p2_max, 'LineWidth', linewidth);
set(xlabel('\it r\rm_{p2} (m)'), 'FontSize', fontsize);
set(ylabel('\it c\rm_{s,p2}/ \it c^{max}_{s,p2} (mol/m^3)'), 'FontSize', 6);
set(title('p1|p2 adjacent'), 'FontSize', 6);
set(gca, 'FontSize', fontsize);

% cs at x=0+1
subplot(2,2,4);
set(gcf,'DefaultAxesColorOrder',colmat)
plot(nodesdimRp2./1e-6, Csp_dim2(1,1:N_r)'./cs_p2_max, '-kx', 'MarkerSize', markersize,'LineWidth', linewidth);hold on;
plot(nodesdimRp2./1e-6, Csp_dim2(1:every:end,1:N_r)./cs_p2_max, 'LineWidth', linewidth);
plot(nodesdimRp2./1e-6, Csp_dim2(end,1:N_r) ./cs_p2_max, 'LineWidth', linewidth);
set(xlabel('\it r\rm_{p2} (m)'), 'FontSize', fontsize);
set(ylabel('\it c\rm_{s,p2}/ \it c^{max}_{s,p2} (mol/m^3)'), 'FontSize', 6);
set(title('CC'), 'FontSize', 6);
set(gca, 'FontSize', fontsize);

%% Add in the SOC
SOC_p1 = zeros(time_steps,Ns);
SOC_p2 = zeros(time_steps,Ns);
SOCinXRp1 = zeros(time_steps,1);
SOCinXRp2 = zeros(time_steps,1);
for tt = 1:time_steps
    for ii = 1:Ns
        cp1 =  nodesdimRp1.^2 .*(1- csp1_nondim(tt,(  ii*N_r - (N_r-1) : N_r*ii ) ) );
        SOC_p1(tt,ii)  = 3 ./  (Rs_p1 ^3)  .* trapz(nodesdimRp1,cp1);
        cp2 =  nodesdimRp2.^2 .*(1- csp2_nondim(tt,(  ii*N_r - (N_r-1) : N_r*ii ) ) );
        SOC_p2(tt,ii)  = 3 ./  (Rs_p2 ^3)  .* trapz(nodesdimRp2,cp2);
    end
end
nodesones = linspace(0,1, length(nodes_p1dim(2:N-1)));
for i = 1:time_steps
    SOCinXRp1(i,1)  =  trapz(nodesones,SOC_p1(i,:)) ;
    SOCinXRp2(i,1)  =  trapz(nodesones,SOC_p2(i,:)) ;
end
%% SOC combined
low_NMC = 0.2721;
ratio_scaled_SOCNMCmax =( max(SOCinXRp1) - (1-high_NMC) )/ (high_NMC - low_NMC);
ratio_scaled_SOCNMCmin = ( min(SOCinXRp1) - (1-high_NMC) )/ (high_NMC - low_NMC);
SOCp1 =   rescale(SOCinXRp1, ratio_scaled_SOCNMCmin,ratio_scaled_SOCNMCmax) ;
SOCinboth = 0.5.*SOCp1 + 0.5.*SOCinXRp2;

colmat = parula(5);
figure;
set(gcf,'DefaultAxesColorOrder',colmat)
plot1112=plot(t_sec,SOCp1,'LineWidth', linewidth);  hold on;
plot1113=plot(t_sec,SOCinXRp2,'LineWidth', linewidth);
plot1111=plot(t_sec,SOCinboth,'LineWidth', linewidth);
set(ylabel('SOC (-)'), 'FontSize', fontsize);
set(xlabel('Time (s)'), 'FontSize', fontsize);
set(legend([plot1111 plot1112 plot1113 ],{'NMC','LFP', 'Combined'}, 'Location','best'), 'FontSize',8);
set(gca, 'FontSize', fontsize);

    function out = rhs(~, X)
        %solves PDE equations using finite central differences and control volumes
        %% Initialise time derivatives and functions
        Fphiep1= zeros(N, 1);
        Fphiep2= zeros(N, 1);
        Fphies = zeros(N, 1);
        cep1_dt= zeros(N, 1);
        cep2_dt= zeros(N, 1);
        ces_dt= zeros(N, 1);
        Fphisp1= zeros(N, 1);
        Fphisp2= zeros(N, 1);

        %% Vectors of PDE quantities
        c_sp1 = X(idx_csp1,1);
        c_sp2 = X(idx_csp2,1);
        phisp1 = X(idx_phisp1,1);
        phisp2 = X(idx_phisp2,1);
        phiep1 = X(idx_phiep1,1);
        phiep2  = X(idx_phiep2,1);
        phies = X(idx_phies,1);
        cep1 =  X(idx_cep1,1);
        cep2 =  X(idx_cep2,1);
        ces =   X(idx_ces,1);
        Jbar_p1 =  X(idx_Jp1,1);
        Jbar_p2  =  X(idx_Jp2,1);

        %% Fluxes and at faces
        grad_phisp1 = diff(phisp1) .* inv_dx_p1;
        grad_phisp2= diff(phisp2) .* inv_dx_p2;

        Kappa_e_p1_bar = kap_p1 .* Ke_bulk_fun(cep1);
        Kappa_p1_faces =(Kappa_e_p1_bar(2:N) + Kappa_e_p1_bar(1:N-1)) ./ 2;
        grad_phiep1 = diff(phiep1) .* inv_dx_p1;
        flux_phiep1 = - Kappa_p1_faces .* grad_phiep1;
        gradlogcep1 = diff(cep1) .* (inv_dx_p1 ./ cep1(2:N));
        fluxcep1log = -  Kappa_p1_faces .* P .* gradlogcep1;
        D_e_p1_bar = De_paras_p1 .* De_bulk_fun(cep1);
        De_p1_faces =(D_e_p1_bar(2:N) + D_e_p1_bar(1:N-1)) ./ 2;
        grad_cep1 = diff(cep1) .* inv_dx_p1;
        flux_cep1 = - De_p1_faces .* grad_cep1;

        Kappa_e_p2_bar = kap_p2 .* Ke_bulk_fun(cep2);
        Kappa_p2_faces =(Kappa_e_p2_bar(2:N) + Kappa_e_p2_bar(1:N-1)) ./ 2;
        grad_phiep2 = diff(phiep2) .* inv_dx_p2;
        flux_phiep2 = - Kappa_p2_faces .* grad_phiep2;
        gradlogcep2 = diff(cep2) .* (inv_dx_p2 ./ cep2(2:N));
        fluxcep2log = -  Kappa_p2_faces .* P .* gradlogcep2;
        D_e_p2_bar = De_paras_p2 .* De_bulk_fun(cep2);
        De_p2_faces =(D_e_p2_bar(2:N) + D_e_p2_bar(1:N-1)) ./ 2;
        grad_cep2 = diff(cep2) .* inv_dx_p2;
        flux_cep2 = - De_p2_faces .* grad_cep2;

        Kappa_e_s_bar = kap_s .* Ke_bulk_fun(ces) ;
        Kappa_s_faces =(Kappa_e_s_bar(2:N) + Kappa_e_s_bar(1:N-1)) ./ 2;
        grad_phies = diff(phies) .* inv_dx_s;
        flux_phies = - Kappa_s_faces .* grad_phies ;
        gradlogces = diff(ces) .* (inv_dx_s ./ ces(2:N));
        fluxceslog = -  Kappa_s_faces .* P .* gradlogces;
        D_e_s_bar = Des_paras .* De_bulk_fun(ces) ;
        De_s_faces =(D_e_s_bar(2:N) + D_e_s_bar(1:N-1)) ./ 2;
        grad_ces = diff(ces) .* inv_dx_s;
        flux_ces = - De_s_faces .* grad_ces ;

        %% cs - at each x and r point ie P2D
        c_sp1_mat = reshape(c_sp1, N_r, Ns);
        c_sp2_mat = reshape(c_sp2, N_r, Ns);
        grad_csp1 = diff(c_sp1_mat, 1, 1) ./ drp1;
        grad_csp2 = diff(c_sp2_mat, 1, 1) ./ drp2;
        csp1_dt_mat = zeros(N_r, Ns);
        csp2_dt_mat = zeros(N_r, Ns);
        csp1_dt_mat(2:N_r-1,:) = -B1 .* diff(-(rp1_sq_outer .* grad_csp1), 1, 1) ./ (drp1 .* rp1_sq_inner);
        csp2_dt_mat(2:N_r-1,:) = -B2 .* diff(-(rp2_sq_outer .* grad_csp2), 1, 1) ./ (drp2 .* rp2_sq_inner);
        csp1_dt_mat(1,:) = grad_csp1(1,:);
        csp2_dt_mat(1,:) = grad_csp2(1,:);
        csp1_dt_mat(N_r,:) = grad_csp1(N_r-1,:) + Ff1 .* Jbar_p1';
        csp2_dt_mat(N_r,:) = grad_csp2(N_r-1,:) + Ff2 .* Jbar_p2';
        cs_p1_surf = c_sp1_mat(end,:)';
        cs_p2_surf = c_sp2_mat(end,:)';
        csp1_dt = csp1_dt_mat(:);
        csp2_dt = csp2_dt_mat(:);

        %% phisp1
        Fphisp1(2:N-1) = diff(grad_phisp1) .* inv_dx_p1 - H1.*Jbar_p1;
        Fphisp1(1) = grad_phisp1(1);
        Fphisp1(N) = grad_phisp1(N-1) - sigp1 * grad_phisp2(1);

        %% phisp2
        Fphisp2(2:N-1) = diff(grad_phisp2) .* inv_dx_p2 - H2.*Jbar_p2;
        Fphisp2(1) = phisp2(1) - phisp1(N) ;
        Fphisp2(N) = grad_phisp2(N-1) - H2 * I_sign ;

        %% ces
        ces_dt(2:N-1) = -inv_eps_e_s .* diff(flux_ces) .* inv_Dx_s_mid;
        ces_dt(1) = flux_ces(1) + S_LP * I_sign ;
        ces_dt(N) = flux_ces(N-1) -  flux_cep1(1)  ;

        %% cep1
        cep1_dt(2:N-1) = -inv_eps_e_p1 .* diff(flux_cep1) .* inv_Dx_p1_mid + Mm1 .* Jbar_p1;
        cep1_dt(1) = cep1(1) - ces(N);
        cep1_dt(N) = flux_cep1(N-1) -  flux_cep2(1)   ;

        %% cep2
        cep2_dt(2:N-1) = -inv_eps_e_p2 .* diff(flux_cep2) .* inv_Dx_p2_mid + Mm2 .* Jbar_p2;
        cep2_dt(1)  = cep2(1) -  cep1(N) ;
        cep2_dt(N)  = grad_cep2(N-1) ;

        %% phies
        Fphies(2:N-1) = diff(-flux_phies + fluxceslog) .* inv_Dx_s_mid;
        Fphies(1) = phies(1)  ;
        Fphies(N) = - flux_phies(N-1) + fluxceslog(N-1) - Ss * I_sign;

        %% phiep1
        Fphiep1(2:N-1) = diff(- flux_phiep1 + fluxcep1log) .* inv_Dx_p1_mid + Ss .* Jbar_p1 ;
        Fphiep1(1) = phiep1(1) - phies(N);
        Fphiep1(N) = - flux_phiep1(N-1) + fluxcep1log(N-1) + flux_phiep2(1) - fluxcep2log(1);

        %% phiep2
        Fphiep2(2:N-1) = diff(- flux_phiep2 + fluxcep2log) .* inv_Dx_p2_mid + Ss .* Jbar_p2 ;
        Fphiep2(1) = phiep2(1) - phiep1(N);
        Fphiep2(N) =  grad_phiep2(N-1) ;

        %% U and J functions
        Ubarp1  =  FonRT .* U_Dischargep1(cs_p1_surf);
        Ubarp2  =  FonRT .* U_Dischargep2(cs_p2_surf);
        FJp1  = Jbar_p1  - U1 .* sqrt(cep1(2:N-1)  .* cs_p1_surf  .* (1 - cs_p1_surf )) .* sinh( 0.5 .* (phisp1(2:N-1) - phiep1(2:N-1) - Ubarp1));
        FJp2 = Jbar_p2  - U2 .* sqrt(cep2(2:N-1)  .* cs_p2_surf  .* (1 - cs_p2_surf )) .* sinh( 0.5 .* (phisp2(2:N-1) - phiep2(2:N-1)  - Ubarp2));

        %% Collect derivatives and Functions
        out = [csp1_dt; csp2_dt; Fphisp1; Fphisp2;  Fphiep1; Fphiep2; Fphies; cep1_dt; cep2_dt; ces_dt;  FJp1 ;FJp2];
    end
    function [value,isterminal,direction] = DFN_event(~,phi_sol)
        %V
        Vevent = (R * T_ref / Fara) * real(phi_sol(idx_phisp2(end))) + Rc*current_A* I_sign;
        value       = v_max - Vevent; %charge
        % value        =  Vevent  - v_min; %discharge
        isterminal   =   1;
        direction    =  -1;
    end
end