%% 19 September 2019 Miroslav Gasparek
% The default parameters of the model of the interaction 
% of the Luteinizing Hormonone Releasing Hormone (LHRH), 
% Luteinizing Hormone (LH) and Testosterone (T)
% 
% 
% The parameters are taken from: 
% (1) Smith, W. R. (1980). HYPOTHALAMIC REGULATION OF PITUITARY SECRETION OF LUTEINIZING HORMONE-
% II FEEDBACK CONTROL OF GONADOTROPIN SECRETION*. 
% Bulletin oj Mathematical Biology (Vol. 42). 
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

function pars = test_model_parameters()
    pars.c = 100; % ng/(ml h), maximum LHRH secretion rate
    pars.g1 = 10; % h^-1,  LH secretion rate
    pars.g2 = 0.7; % h^-1, T secretion rate
    pars.b1 = 1.29; % h^-1, LHRH decay rate
    pars.b2 = 0.97; % h^-1, LH decay rate
    pars.b3 = 1.39; % h^-1, T decay rate
    pars.h = 12; % h^-1, (roughly) proportionality constant for the
                 % strength of the testestosterone inhibition
                 
    pars.wR = 0; % ng/(ml h), influx rate of the exogenous LHRH
    pars.wL = 0; % ng/(ml h), influx rate of the exogenous LH
    pars.wT = 0; % ng/(ml h), influx rate of the exogenous T

end

