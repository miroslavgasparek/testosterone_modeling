%% 25 September 2019 Miroslav Gasparek
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
% In this script, we analyze the effect of onset of the andropause, 
% which is modeled through the decrease in the pitituary secretion rate of LH and gonadal
% secretion rate of T.

clear;clc; close all;

addpath('subroutines');
fprintf('Required subroutines added to path.\n====================\n\n');

%% Modeling the effect of the onset of andropause 
%%% Initial set up %%%
% Set parameters
pars_def = test_model_parameters();
pars = pars_def;

%%{

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

%%% Running the loop %%% 
% Define the number of points
num_points = 50;

% Define the range of the pitituary sensitivities relative to 
% physiological pitituary sensitivity (1/50 up to 1)
g1_range = linspace(pars.g1/num_points, pars.g1, num_points); % 1/h

% Define the range of the gonadal sensitivities relative to 
% physiological gonadal sensitivity (1/50 up to 1)
g2_range = linspace(pars.g2/num_points, pars.g2, num_points); % 1/h

% Get the empty matrix for the mean hormone concentrations
LHRH_mean_mat = zeros(size(g1_range,2), size(g2_range,2));
LH_mean_mat = zeros(size(g1_range,2), size(g2_range,2));
T_mean_mat = zeros(size(g1_range,2), size(g2_range,2));

% Get the matrix to check whether the Test. concentration oscillates
T_stab_mat = zeros(size(g1_range,2), size(g2_range,2));

% Define the values for the checking of occurence of the hormonal conc.
% oscillations and for the computation of the mean value 
% of the hormones
frac_var = 0.9;
frac_mean = 0.5;
var_thres = 2.0;

% Run the simulations
for i = 1:length(g1_range)
    pars.g1 = g1_range(i);
    for j = 1:length(g2_range)
        pars.g2 = g2_range(j);
        [tout, yout, teout, yeout, ieout] = test_solve_ode([tstart, tfinal], y0, pars);
        
        % Check if the output oscillates and compute the mean steady-state
        % Testosterone levels
        [~, LHRH_mean] = check_steady_state(yout(:,1), frac_var, frac_mean, var_thres);
        [~, LH_mean] = check_steady_state(yout(:,2), frac_var, frac_mean, var_thres);
        [T_iout, T_mean] = check_steady_state(yout(:,3), frac_var, frac_mean, var_thres);
    
        % Store the computed results into the matrix
        LHRH_mean_mat(i,j) = LHRH_mean;
        LH_mean_mat(i,j) = LH_mean;
        T_mean_mat(i,j) = T_mean;
        T_stab_mat(i,j) = T_iout;
    end
    i
end

% Plotting
figure(6)
imagesc(T_mean_mat);
% Select the colormap
colormap bone;
cbar = colorbar;
set(cbar, 'YDir','reverse');

titleString = 'Mean Test. levels (ng/ml)';
ylabel(cbar, titleString,'Fontsize',15,'interpreter','latex');

% Set the ticks of the axes
xticklabels = round(g1_range(2:4:end),2);
xticks = linspace(2,size(T_mean_mat,2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);

yticklabels = g2_range(2:4:end);
yticks = linspace(2,size(T_mean_mat,2), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);

fig = gcf;
fig.Position = [291   280   709   518];
ax = gca;
ax.XAxis.TickLabelFormat = '%.2f';
ax.YAxis.TickLabelFormat = '%.2f';
ax.FontSize=15;

xlabel('Relative LH Secretion Rate $(h^{-1})$','fontsize',20,'interpreter','latex')
ylabel('Relative Testosterone Secretion Rate $(h^{-1})$','fontsize',20,'interpreter','latex')
title('Mean steady-state Testosterone level for different LHRH and LH secretion rates','Fontsize',18,'interpreter','latex')

% Plotting of the oscillation colormap
figure(7)
imagesc(T_stab_mat);
% Select the colormap
cmap = gray(2);
colormap(cmap);
cbar1 = colorbar;
cbar1.Ticks = [0.25, 0.75];
cbar1.TickLabels = {'Oscillations','Steady State'};
cbar1.TickLabelInterpreter = 'latex';

% Set the ticks of the axes
% Set the ticks of the axes
xticklabels = round(g1_range(2:4:end),2);
xticks = linspace(2,size(T_mean_mat,2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);

yticklabels = g2_range(2:4:end);
yticks = linspace(2,size(T_mean_mat,2), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);

fig = gcf;
fig.Position = [291   280   709   518];
ax = gca;
ax.XAxis.TickLabelFormat = '%.2f';
ax.YAxis.TickLabelFormat = '%.2f';
ax.FontSize=15;

xlabel('Relative LH Secretion Rate $(h^{-1})$','fontsize',20,'interpreter','latex')
ylabel('Relative Testosterone Secretion Rate $(h^{-1})$','fontsize',20,'interpreter','latex')

%}
%% Some weird results for certain parameters values
% Try to simulate this, you will see some interesting behavior! 
% g1_weird = 6.3056;
% g2_weird = 0.1828;
% pars.g1 = g1_weird;
% pars.g2 = g2_weird;
