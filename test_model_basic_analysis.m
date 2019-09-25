%% 24 September 2019 Miroslav Gasparek
%%% Modeling of the interaction of Luteinizing Hormonone Releasing
%%% Hormone (LHRH), Luteinizing Hormone (LH) and Testosterone (T)
%
% The work is based on the following paper:
%
% (1) Smith, W. R. (1980). HYPOTHALAMIC REGULATION OF PITUITARY SECRETION OF LUTEINIZING HORMONE-
% II FEEDBACK CONTROL OF GONADOTROPIN SECRETION*. 
% Bulletin of Mathematical Biology (Vol. 42). 
% Retrieved from https://link.springer.com/content/pdf/10.1007%2FBF02462366.pdf
%
% Model equations have the following form
%
%   dR/dt = c - h * T * H(1 - (T - c/h)) - b1 * R
%   dL/dt = g1 * R - b2 * L
%   dT/dt = g2 * L - b3 * T
% 
% Where H(x) is a Heaviside step function:
% H(x <= 0) = 0
% H(x >  0) = 1
% 
% c, h, b1, b2, b3, g1, g2 are constants
%
% In this script, we conduct the basic analysis of the testosterone
% secretion system, considering the normal physiological behavior,
% effect of the castration, and effect of the ext. testosterone influx.

clear;clc;close all;

addpath('subroutines');
fprintf('Subroutines added to path.\n====================\n\n');

%% Basic testosterone secretion model with pulsatile secretion
%%% Store the default values of parameters %%%
pars_def = test_model_parameters();
pars = pars_def;

%%% Define the duration of simulation %%%
tstart = 0; % Hours
tfinal = 48; % Duration of simulation in hours

%%% Initial concentrations of hormones %%%
LHRH_init = 1; % ng/ml
LH_init = 25; % ng/ml
T_init = 5; % ng/ml

y0 = [LHRH_init; 
      LH_init; 
      T_init];

%-------------------%

% Run the simulation
[tout, yout, teout, yeout, ieout] = test_solve_ode([tstart, tfinal], y0, pars);

% Display results
fprintf('LHRH-LH-T secretion model\n=========================================== \n\n')
fprintf('The mean testosterone value: %.2f ng/ml \n', mean(yout(:,3)));


% Plotting
% Set up the colors for plotting of each hormone
c_LHRH = [0, 0, 0];
c_LH = [0.7, 0, 0];
c_T = [0, 0.5, 0.9];

figure(1)
hold on
plot(tout,yout(:,1),'Color', c_LHRH,'LineWidth',2)
plot(tout, yout(:,2),'Color', c_LH,'LineWidth',2)
plot(tout, yout(:,3),'Color', c_T,'LineWidth',2)
hold off

fig = gcf;
fig.Position = [291   280   709   518];
ax = gca;
ax.FontSize=15;

xlabel('Time (h)','fontsize',20,'interpreter','latex');
ylabel('Concentrations','fontsize',20,'interpreter','latex');
title('Time evolution of the Testosterone secretion system','fontsize',20,'interpreter','latex');
legend('LHRH (ng/ml)', 'LH (ng/ml)', 'T (ng/ml)','fontsize',15,'interpreter','latex');

% Plot the "phase plane" of LH and T
figure(2)
plot(yout(:,2), yout(:,3),'Color','black','LineWidth', 2) 

fig = gcf;
fig.Position = [291   280   709   518];
ax = gca;
ax.FontSize=15;

xlabel('LH (ng/ml)','fontsize',20,'interpreter','latex');
ylabel('T (ng/ml)','fontsize',20,'interpreter','latex');
title('Testosterone vs. Luteinizing Hormone concentrations','fontsize',20,'interpreter','latex');



%% Testosterone secretion model in the case of castration
% In the case of castration, the production rate of testosterone in 
% the gonadal glands is zero
pars.g2 = 0; % h^-1

% Run the simulation
[tout, yout, teout, yeout, ieout] = test_solve_ode([tstart, tfinal], y0, pars);

% Display results
fprintf('LHRH-LH-T secretion model (castration) \n')
fprintf('The mean testosterone value: %.2f ng/ml \n', mean(yout(:,3)));

% Plotting
figure(3)
hold on
plot(tout,yout(:,1),'Color', c_LHRH,'LineWidth',2)
plot(tout, yout(:,2),'Color', c_LH,'LineWidth',2)
plot(tout, yout(:,3),'Color', c_T,'LineWidth',2)
hold off

fig = gcf;
fig.Position = [291   280   709   518];
ax = gca;
ax.FontSize=15;

xlabel('Time (h)','fontsize',20,'interpreter','latex');
ylabel('Concentrations','fontsize',20,'interpreter','latex');
title('Time evolution of the Testosterone secretion system (castration)','fontsize',20,'interpreter','latex');
legend('LHRH (ng/ml)', 'LH (ng/ml)', 'T (ng/ml)','fontsize',15,'interpreter','latex');


%% Testosterone secretion model in the case of external testosterone input
% In the case of external input, wT is non-zero
pars = pars_def;
pars.wT = 20; % ng/(ml h)
[tout, yout, teout, yeout, ieout] = test_solve_ode([tstart, tfinal], y0, pars);

% Check if the Testosterone output oscillates and compute its 
% steady state mean
frac_var = 0.8;
frac_mean = 0.5;
var_thres = 1.0;
[iout, T_mean] = check_steady_state(yout(:,3), frac_var, frac_mean, var_thres);

% Display results
fprintf('LHRH-LH-T secretion model (ext. input), wT = %.2f ng/(ml h) \n', pars.wT);
fprintf('The mean testosterone value: %.2f ng/ml \n', T_mean);

% Print out the stability of the system
if iout == 0
    fprintf('Solution oscillates (limit cycle).\n')
elseif iout == 1
    fprintf('Solution approaches steady-state (fixed point).\n')
end

% Plotting
figure(4)
hold on
plot(tout,yout(:,1),'Color', c_LHRH,'LineWidth',2)
plot(tout, yout(:,2),'Color', c_LH,'LineWidth',2)
plot(tout, yout(:,3),'Color', c_T,'LineWidth',2)
hold off

fig = gcf;
fig.Position = [291   280   709   518];
ax = gca;
ax.FontSize=15;

xlabel('Time (h)','fontsize',20,'interpreter','latex');
ylabel('Concentrations','fontsize',20,'interpreter','latex');
title({['Time evolution of the Testosterone secretion system,'],...
    ['ext. testosterone influx: $w_{T} = ',num2str(pars.wT),' \ ng \ ml^{-1} \ h^{-1}$']},'fontsize',20,'interpreter','latex');
legend('LHRH (ng/ml)', 'LH (ng/ml)', 'T (ng/ml)','fontsize',15,'interpreter','latex');


%% Dependence of the mean testosterone level on the external test. input
pars = pars_def;
% Define the range of external inputs
wT_range = 0:0.1:20; % ng/(ml h)

% Get the empty vector for the mean hormone concentrations
LHRH_mean_vec = zeros(size(wT_range));
LH_mean_vec = zeros(size(wT_range));
T_mean_vec = zeros(size(wT_range));

% Set the longer timespan 
tspan = [0, 72]; % hours

% Calculate the value of ext. test. at which the LHRH production is shut
wT_crit = (pars.c * pars.b3)/pars.h;

% Run the simulations to obtain the mean value of the testosterone
for i = 1:length(wT_range)
    pars.wT = wT_range(i);
    [tout, yout, teout, yeout, ieout] = test_solve_ode(tspan, y0, pars);
    LHRH_mean_vec(i) = mean(yout(:,1));
    LH_mean_vec(i) = mean(yout(:,2));
    T_mean_vec(i) = mean(yout(:,3));
end

% Plotting
figure(5)
hold on
plot(wT_range, LHRH_mean_vec,'Color', c_LHRH,'LineWidth',2)
plot(wT_range, LH_mean_vec,'Color', c_LH,'LineWidth',2)
plot(wT_range, T_mean_vec,'Color', c_T,'LineWidth',2)

% Plot the vertical line at the critical value of ext. test. input
YLim = get(gca, 'YLim');
plot([wT_crit, wT_crit], YLim,'k--','LineWidth',3)
hold off

fig = gcf;
fig.Position = [291   280   709   518];
ax = gca;
ax.FontSize=15;

xlabel('External Testosterone influx rate, $w_{T} \ (ng \ ml^{-1} \ h^{-1})$','fontsize',20,'interpreter','latex')
ylabel('Mean horomone concentration','fontsize',20,'interpreter','latex')
title('Mean hormone concentrations for varying external Testosterone input ','fontsize',20,'interpreter','latex')
legend('$\langle LHRH \rangle$ (ng/ml)', '$\langle LH \rangle$ (ng/ml)', '$\langle T \rangle$ (ng/ml)','fontsize',15,'interpreter','latex');