%% 18 September 2019 Miroslav Gasparek
% ODE Model of the testosterone secretion system

function dydt = test_ode_pulse(t, y, pars)
% Function describes the ODEs model of the pulsatile testosterone
% secretion with the (potential) constant external input wT
% The model is based on:
% (1) Smith, W. R. (1980). HYPOTHALAMIC REGULATION OF PITUITARY SECRETION OF LUTEINIZING HORMONE-
% II FEEDBACK CONTROL OF GONADOTROPIN SECRETION*. 
% Bulletin oj Mathematical Biology (Vol. 42). 
% Retrieved from https://link.springer.com/content/pdf/10.1007%2FBF02462366.pdf

    % Initialize the function to zeros
    dydt = zeros(3,1);
    
    % Define the system of the ODEs
    dydt = [(pars.wR + pars.c - pars.h * y(3)) * (1 - pars.flag) - pars.b1 * y(1);
            pars.wL + pars.g1 * y(1) - pars.b2 * y(2);
            pars.wT + pars.g2 * y(2) - pars.b3 * y(3)];
        
end
