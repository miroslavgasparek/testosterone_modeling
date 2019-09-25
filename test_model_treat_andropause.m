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
% In this script, we look at the manipulation of the levels of 
% hormones through the external (constant) dosing of these hormones.
% We look at the external influx of Testosterone (T)
% and Luteinizing Hormone (LH)

clear;clc;close all;

addpath('subroutines');
fprintf('Subroutines added to path.\n====================\n\n');

%% The example simulation of the ext. Testosterone influx
% First initialize the standard set of parameters
pars_def = test_model_parameters();
pars = pars_def;

% Pick the parameters values for which the levels of Testosterone are below
% the clinically plausible threshold
  
pars.g1 = 1.0278;
pars.g2 = 0.0719;
pars.wT = 15; % ng/(ml h)

% Set the timespan 
tstart = 0; % hours
tfinal = 48; % hours 

%%% Initial concentrations of hormones %%%
LHRH_init = 1; % ng/ml
LH_init = 25; % ng/ml
T_init = 5; % ng/ml

y0 = [LHRH_init; 
      LH_init; 
      T_init];

% Integrate the ODE system
[tout, yout, teout, yeout, ieout] = test_solve_ode([tstart, tfinal], y0, pars);

% Check if the Testosterone level oscillates and calculate its mean value
% at the selected interval
frac_var = 0.9;
frac_mean = 0.7;
var_thres = 1.0;
[T_iout, T_mean] = check_steady_state(yout(:,3), frac_var, frac_mean, var_thres);

% Display the simulation results
fprintf('LHRH-LH-T secretion model (ext. Test. influx)\n=========================================== \n\n')
fprintf('Mean Testosterone value: %.2f ng/ml \n', T_mean);

% Calculate the threshold for the expected steady-state of Testosterone 
wT_min = pars.c*pars.b3/pars.h;

% Display the threshold
fprintf('Min.ext. Testosterone influx to drive LHRH & LH to zero: %.2f ng/(ml h) \n\n', wT_min);
% Display the results predicted by the model
fprintf('Expected mean Testosterone value: %.2f ng/ml \n', pars.wT/pars.b3);

% Plotting
% Set up the colors for plotting of each hormone
c_LHRH = [0, 0, 0];
c_LH = [0.7, 0, 0];
c_T = [0, 0.5, 0.9];

% Plot the levels of the hormones and the predicted steady-state of
% Testosterone
figure(8)
hold on
plot(tout,yout(:,1),'Color', c_LHRH,'LineWidth',2);
plot(tout, yout(:,2),'Color', c_LH,'LineWidth',2);
plot(tout, yout(:,3),'Color', c_T,'LineWidth',2);
plot(tout,pars.wT/pars.b3*ones(1,length(tout)),'k--','LineWidth',3);
xlabel('Time (h)','fontsize',20,'interpreter','latex');
ylabel('Concentrations','fontsize',20,'interpreter','latex');
title({['Time evolution of the Testosterone secretion system,'],...
    ['$g_{1} = ',num2str(pars.g1),'\ h^{-1}$, $g_{2} = ',...
    num2str(pars.g2),'\ h^{-1}$, $w_{T}$ = ',num2str(pars.wT),...
    ' $ng \ ml^{-1} \ h^{-1}$']},'fontsize',20,'interpreter','latex');
hold off
legend('LHRH (ng/ml)', 'LH (ng/ml)', 'T (ng/ml)','fontsize',15,'interpreter','latex');

fig = gcf;
fig.Position = [291   280   709   518];

%% Effect of the varying external Testosterone influx
% Set the parameters to their default values 
pars = pars_def;

% Set the timespan 
tstart = 0; % hours
tfinal = 48; % hours 

%%% Initial concentrations of hormones %%%
LHRH_init = 1; % ng/ml
LH_init = 25; % ng/ml
T_init = 5; % ng/ml

y0 = [LHRH_init; 
      LH_init; 
      T_init];
  
  
% Select the vector of values of the ext. Testosterone inputs
wT_range = [9, 12, 15, 20];

% Calculate the threshold for the expected steady-state of Testosterone 
wT_min = pars.c*pars.b3/pars.h;

% Integrate the ODE for the given value of ext. influx wT
% and simultaneously plot the figure

figure(9)
hold on
for i = 1:length(wT_range)
    
    pars.wT = wT_range(i);
    % Integrate the ODE system
    [tout, yout, teout, yeout, ieout] = test_solve_ode([tstart, tfinal], y0, pars);
    
    % Plot the figure
    txt = ['$w_{T} = ',num2str(pars.wT),'\ ng \ ml^{-1} \ h^{-1}$'];
    p1 = plot(tout, yout(:,3),'LineWidth',2,'DisplayName',txt);
    c = get(p1, 'Color');
    p2 = plot(tout,pars.wT/pars.b3*ones(1,length(tout)),'--','Color',c, 'LineWidth',3,'HandleVisibility','off');
    
end

plot(tout, wT_min/pars.b3*ones(1,length(tout)),'k-.','LineWidth',3,'DisplayName','Steady-state thres.');
hold off

% Set the legend
l = legend;
set(l,'fontsize',15, 'interpreter','latex')


fig = gcf;
fig.Position = [291   280   709   518];
ax = gca;
ax.XLim = [0, tout(end)];
ax.FontSize=15;

% Set the axis labels and the title
xlabel('Time (h)','fontsize',20,'interpreter','latex');
ylabel('Testosterone concentration (ng/ml)','fontsize',20,'interpreter','latex');
title('Testosterone concentration for varying external Testosterone influx rates','fontsize',20,'interpreter','latex')


%% Effect of the different values of LH and T secretion rates with ext. influx
% Set the parameters to their default values 
pars = pars_def;
pars.wT = 15; % ng / (ml h)

% Set the timespan 
tstart = 0; % hours
tfinal = 48; % hours 

%%% Initial concentrations of hormones %%%
LHRH_init = 1; % ng/ml
LH_init = 25; % ng/ml
T_init = 5; % ng/ml

y0 = [LHRH_init; 
      LH_init; 
      T_init];

% Define the range of the pitituary sensitivities relative to 
% physiological pitituary sensitivity (1/50 up to 1)
g1_range_extT = [1, 10]; % 1/h

% Define the range of the gonadal sensitivities relative to 
% physiological gonadal sensitivity (1/50 up to 1)
g2_range_extT = [0.07, 0.7]; % 1/h

% Run the simulation and plot the testosterone levels
figure(10)
hold on
for i = 1:length(g1_range_extT)
    
    pars.g1 = g1_range_extT(i);
    for j = 1:length(g2_range_extT)
    
        
        pars.g2 = g2_range_extT(j);
        % Integrate the ODE system
        [tout, yout, teout, yeout, ieout] = test_solve_ode([tstart, tfinal], y0, pars);

        % Plot the figure
        txt = ['$g_{1} = ',num2str(pars.g1),'\ h^{-1},$',' $g_{2} = ',num2str(pars.g2),'\ h^{-1}$'];
        plot(tout, yout(:,3),'LineWidth',2,'DisplayName',txt);
    end
end
plot(tout,pars.wT/pars.b3*ones(1,length(tout)),'--','Color','Black', 'LineWidth',3);
hold off

% Set the legend
l = legend;
set(l,'fontsize',15, 'interpreter','latex')

fig = gcf;
fig.Position = [291   280   709   518];
ax = gca;
ax.XLim = [0, tout(end)];
ax.FontSize=15;
ax.Legend.String{5} = 'Expected steady-state T level';

% Set the axis labels and title
xlabel('Time (h)','fontsize',20,'interpreter','latex');
ylabel('Testosterone concentration (ng/ml)','fontsize',20,'interpreter','latex');
title({['Testosterone levels for varying pitituary gonadal secretion rates,'],...
    ['external Testosterone input: $w_{T} = ',...
    num2str(pars.wT),'\ ng \ ml^{-1} \ h^{-1}$']},'fontsize',20,'interpreter','latex')

%% Effect of the constant external influx of all the hormones
% Set the values of the pitituary/gonadal secretion rates
pars = pars_def;

% Set the timespan 
tstart = 0; % hours
tfinal = 48; % hours 

%%% Initial concentrations of hormones %%%
LHRH_init = 1; % ng/ml
LH_init = 25; % ng/ml
T_init = 5; % ng/ml

y0 = [LHRH_init; 
      LH_init; 
      T_init];
  
% Reduce the pitituary/gondal secretion rates
f_andro_LH = 10;
f_andro_T = 10;
pars.g1 = pars.g1/f_andro_LH;
pars.g2 = pars.g2/f_andro_T;

% In this case, we aim to achieve the following concentrations
% of the hormones
LHRH_target = 0; % ng/ml
LH_target = 25; % ng/ml
T_target = 15; % ng/ml

% Calculate the appropriate external fluxes
% wR_target = LHRH_target * pars.b1;
% wL_target = LH_target * pars.b2 - pars.b2 * pars.g1 * LHRH_target;
% wT_target = pars.b3 * T_target - (pars.g1 * pars.g2 * LHRH_target)/pars.b2 + pars.g1 * LHRH_target;

wR_target = LHRH_target * pars.b1;
wL_target = LH_target * pars.b2 - pars.g1 * LHRH_target;
wT_target = pars.b3 * T_target - pars.g2 * LH_target;

% Display the results
fprintf('LHRH-LH-T secretion model (ext. influx of all hormones)\n=========================================== \n\n')

% Print the reduction of the secretion rates of LH and T
fprintf('The LH secretion rate is %.2f times the original secretion rate\n', 1/f_andro_LH)
fprintf('The T secretion rate is %.2f times the original secretion rate\n', 1/f_andro_T)

% Print the desired steady-state values of hormones
fprintf('\nDesired steady-state value of LHRH: %.2f ng/ml \n', LHRH_target);
fprintf('Desired steady-state value of LH: %.2f ng/ml \n', LH_target);
fprintf('Desired steady-state value of T: %.2f ng/ml \n', T_target);

% Print the required values of the external fluxes of the hormones
fprintf('\nExt. LHRH influx needed to get steady-state hormone conc.: %.2f ng/(ml h)\n', wR_target);
fprintf('Ext. LH influx needed to get steady-state hormone conc.: %.2f ng/(ml h)\n', wL_target);
fprintf('Ext. T influx needed to get steady-state hormone conc.: %.2f ng/(ml h)\n', wT_target);

% Check the conditions for the steady-state validity
lhs_condition = pars.g1 * pars.g2 * wR_target +...
    pars.b1 * pars.g2 * wL_target + pars.b1 * pars.b2 * wT_target;

rhs_condition = (pars.c * pars.b1 * pars.b2 * pars.b3)/pars.h;

% Check if the steady state solution is feasible
% If it is feasible, set the external influxes of the hormones
% to their calculated values
if lhs_condition > rhs_condition
    fprintf('\nThe steady-state solution is feasible.\n')
    pars.wR = wR_target;
    pars.wL = wL_target;
    pars.wT = wT_target;
    
else
    fprintf('\nThe steady-state solution is not feasible.\n')
end

% Integrate the ODEs
[tout, yout, teout, yeout, ieout] = test_solve_ode([tstart, tfinal], y0, pars);

% Check if the Testosterone output oscillates and compute its 
% steady state mean
frac_var = 0.8;
frac_mean = 0.5;
var_thres = 1.0;
[iout, T_mean] = check_steady_state(yout(:,3), frac_var, frac_mean, var_thres);

% Plotting
figure(11)
hold on
plot(tout,yout(:,1),'Color', c_LHRH,'LineWidth',2)
plot(tout, yout(:,2),'Color', c_LH,'LineWidth',2)
plot(tout, yout(:,3),'Color', c_T,'LineWidth',2)
plot(tout, LHRH_target*ones(1,length(tout)),'k--','LineWidth',2.5);
plot(tout, LH_target*ones(1,length(tout)),'k--','LineWidth',2.5);
plot(tout, T_target*ones(1,length(tout)),'k--','LineWidth',2.5);
hold off

fig = gcf;
fig.Position = [291   280   709   518];
ax = gca;
ax.FontSize=15;

xlabel('Time (h)','fontsize',20,'interpreter','latex');
ylabel('Concentrations','fontsize',20,'interpreter','latex');
title({['Time evolution of the Testosterone secretion system.'],...
    ['Pitituary secretion rate: ',num2str(1/f_andro_LH),'$g_{1}^{healthy}$',', gonadal secretion rate: ',num2str(1/f_andro_T),'$g_{2}^{healthy}$'],...
    ['External hormone influx rates: $w_{R} = ',num2str(pars.wR),...
    ' \ ng \ ml^{-1} \ h^{-1}, w_{L} = ', num2str(pars.wL),...
    ' \ ng \ ml^{-1} \ h^{-1}, w_{T} = ', num2str(pars.wT), ' \ ng \ ml^{-1} \ h^{-1}$']},...
    'fontsize',15,'interpreter','latex');
legend('LHRH (ng/ml)', 'LH (ng/ml)', 'T (ng/ml)','fontsize',15,'interpreter','latex');


