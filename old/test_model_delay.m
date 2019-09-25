%%% 06 September 2019 Miroslav Gasparek
%%% Simple Delay Differential Equation (DDE) model 
% of the testosterone secretion

% Based on Ruan, S., & Wei, J. (2001). 
% On the zeros of a third degree exponential polynomial 
% with applications to a delayed model for the control of testosterone secretion. 
% IMA Journal of Mathemathics Applied in Medicine and Biology.

clear; clc; close all;

%% Basic initial model of test. production
% Define the parameters of the model
pars.c = 100; % pg/ml
pars.g1 = 10; % h^-1
pars.g2 = 0.7; % h^-1
pars.b1 = 1.29; % h^-1
pars.b2 = 0.97; % h^-1
pars.b3 = 1.39; % h^-1
pars.h = 1; % h^-1

% Define the delay in the effect of LH on the Test secretion
tau0 = 0.2905; % h

% Define length of simulation in days
days = 5;

% Define timespan of the solution in hours
tspan = [0, days*24];

% Define the initial values of the hormone concentrations:
LHRH_init = 12; % pg/ml
LH_init = 100; % pg/ml
T_init = 70; % pg/ml

% Store the initial values in a vector
init_vals = [LHRH_init;
             LH_init;
             T_init];

% Solve the DDE
sol1 = dde23(@(t,y,Z) test_dde(t, y, Z, pars), tau0, @(t) history(t,init_vals), tspan);

% Unpack the 'sol' structure
t_sol = sol1.x;
LHRH_t = sol1.y(1,:);
LH_t = sol1.y(2,:);
T_t = sol1.y(3,:);

disp(mean(T_t))
%%% Plotting
figure(1)
hold on
sgtitle({['Time evolution of the Testosterone secretion system, ','\tau_{0} = ', num2str(tau0),' h']},'Fontsize',15)
subplot(3,1,1)
plot(t_sol, LHRH_t,'LineWidth',2,'Color','k')
xlabel('Time (hours)','fontsize',15)
ylabel('LHRH (pg/ml)','fontsize',15)
subplot(3,1,2)
plot(t_sol, LH_t,'LineWidth',2,'Color','b')
xlabel('Time (hours)','fontsize',15)
ylabel('LH (pg/ml)','fontsize',15)
subplot(3,1,3)
plot(t_sol, T_t,'LineWidth',2,'Color','r')
xlabel('Time (hours)','fontsize',15)
ylabel('Testosterone (pg/ml)','fontsize',15)

fig = gcf;
fig.Position = [440   292   681   506];

%% Modeling the effect of castration
% The effect of castration is such that g2 = 0, i. e. there is no
% testosterone production

% Set the production rate of testosterone to zero
pars.g2 = 0; % h^-1

% Solve the DDE
sol2 = dde23(@(t,y,Z) test_dde(t, y, Z, pars), tau0, @(t) history(t,init_vals), tspan);

% Unpack the 'sol' structure
t_sol = sol2.x;
LHRH_t = sol2.y(1,:);
LH_t = sol2.y(2,:);
T_t = sol2.y(3,:);

%%% Plotting
figure(2)
hold on
sgtitle({['Time evolution of the Testosterone secretion system (castration), ','\tau_{0} = ', num2str(tau0),' h']},'Fontsize',15)
subplot(3,1,1)
plot(t_sol, LHRH_t,'LineWidth',2,'Color','k')
xlabel('Time (hours)','fontsize',15)
ylabel('LHRH (pg/ml)','fontsize',15)
subplot(3,1,2)
plot(t_sol, LH_t,'LineWidth',2,'Color','b')
xlabel('Time (hours)','fontsize',15)
ylabel('LH (pg/ml)','fontsize',15)
subplot(3,1,3)
plot(t_sol, T_t,'LineWidth',2,'Color','r')
xlabel('Time (hours)','fontsize',15)
ylabel('Testosterone (pg/ml)','fontsize',15)

fig = gcf;
fig.Position = [440   292   681   506];

%% Modelling the onset of the puberty


%% Auxillary functions
function dydt = test_dde(t, y, Z, pars)
%%% Function to describe the DDEs 

    % Define the delay approximation
    ylag1 = Z(:,1);
    
    dydt = zeros(3,1);
    
    % Define the system of the DDEs
    dydt = [pars.c - pars.h * y(3) - pars.b1 * y(1);
            pars.g1 * y(1) - pars.b2 * y(2);
            pars.g2 * ylag1(2) - pars.b3 * y(3)];
        
end

function s = history(t, init_vals)
    %%% Function to define the solution history
    
    % Essentially, defines the initial conditions in this case
    s = [init_vals(1); 
         init_vals(2); 
         init_vals(3)];
end