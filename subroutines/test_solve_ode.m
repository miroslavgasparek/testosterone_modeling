%% 18 September 2019 Miroslav Gasparek
% Solving the ODEs function

function [tout, yout, teout, yeout, ieout] = test_solve_ode(tspan, y0, pars)
% Function solves the system of the ODEs describing the pulsatile secretion
% of LHRH, LH, and T.

%%% Inputs %%%
% tspan: [tstart, tfinal] the timespan of the simulation (i. e. starting
% time and the ending time
% y0: Initial concentrations of the hormones
% pars: The set of parameters as described in (1)
    % Get the starting point and endpoint of the simulation
%%% Outputs %%%
% tout: The time variable of the solution
% yout: The levels of the hormones (state variables)
% teout: The times of the switching events
% yeout: The values of the state variables at the switching times
% ieout: The indexes of the switching events

    tstart = tspan(1);
    tfinal = tspan(2);
    
    % Set the ODE solver
    refine = 4;
    options = odeset('Events',@(t, y)heaviside(t, y, pars),...
       'Refine',refine);

    % Set the first input of the time and state vector
    tout = tstart;
    yout = y0.';

    % Define empty vectors for events info
    teout = [];
    yeout = [];
    ieout = [];

    % Set up the switching variable
    pars.flag = 1;
    
    % Set the initial value of the 'flag' variable
    if (y0(3) - pars.c/pars.h) <= 0
        pars.flag = 0;
    end

    % Loop over the events
    while tout(end) < tfinal
       % Solve until the first terminal event.
       [t,y,te,ye,ie] = ode45(@(t,y)test_ode_pulse(t, y, pars),[tstart, tfinal], y0, options);  

       % Switch the value of the Heaviside step function upon crossing the zero
       % before the integration
       if pars.flag == 0
           pars.flag = 1;
       elseif pars.flag == 1
           pars.flag = 0;
       end

       % Accumulate output.  This could be passed out as output arguments.
       nt = length(t);
       tout = [tout; t(2:nt)];
       yout = [yout; y(2:nt,:)];
       teout = [teout; te];          % Events at tstart are never reported.
       yeout = [yeout; ye];
       ieout = [ieout; ie];

       % Set the new initial conditions, starting from the end of the previous
       % integration
       y0(1) = y(nt,1);
       y0(2) = y(nt,2);
       y0(3) = y(nt,3);

       % A good guess of a valid first timestep is the length of the last valid
       % timestep, so use it for faster computation.  'refine' is 4 by default.
       options = odeset(options,'InitialStep',t(nt)-t(nt-refine),...
          'MaxStep',t(nt)-t(1));

       tstart = t(nt);
    end

end