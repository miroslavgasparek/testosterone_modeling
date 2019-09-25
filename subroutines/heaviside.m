%% 18 September 2019 Miroslav Gasparek
% Event function for Heaviside step function

function [value, isterminal, direction] = heaviside(t, y, pars)
    % Function the defines the conditions for the occurence of the
    % switching event according to the ODE model.
    % This is needed for the stopping of the integration at the
    % discontinuity that Heaviside step function represents.
    
    % Essentially, defines the initial conditions in this case
    value = y(3) - pars.c/pars.h;
    isterminal = 1;
    direction = 0;
    
end