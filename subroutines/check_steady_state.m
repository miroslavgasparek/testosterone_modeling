%% 19 September 2019 Miroslav Gasparek
% Function checks if the output signal oscillates 

function [iout, y_mean_ss] = check_steady_state(y_in, frac_var, frac_mean, var_thres)
% Check (roughly) if the input vector oscillates
% by checking the variance of the fraction 'frac_var' of the simulation time
% 
% Also computes the mean value of input vector y_in from the selected sample of
% data at the interval [frac_mean*t_final, t_final]
% 
%%% Inputs %%% 
% y_in: Input vector
% frac_var: Relative length of the simulation considered in the simulation,
% number between 0 and 1
% frac_mean: The relative length of simulation considered for the
% calculation of the y_in mean value.
% var_thres: Maximum admissible variance for the steady state

%%% Outputs %%%
% iout: Index describing whether the system oscillates
% iout = 0 (solution oscillates)
% iout = 1 (solution roughly goes to steady state)
%
% y_mean_ss: The mean value of y
    
    % Select the portion of the y_in vector from which this is computed
    y_for_var = y_in(floor(frac_mean * length(y_in)));
    
    % Evaluate the variance and compare it to threshold
    if var( y_in(y_in >= y_for_var) ) >= var_thres
        iout = 0; % Solution oscillates
    else
        iout = 1; % Solution does not oscillate
    end
    
    %  Compute mean value of y_in at the selected interval [t > frac_mean*tfinal]
    % Select the portion of the y_in vector from which this is computed
    y_for_mean = y_in(floor(frac_mean * length(y_in)));
    y_mean_ss = mean( y_in(y_in >= y_for_mean) );
end