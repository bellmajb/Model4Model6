%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  The following program solves JacobsCoupled  ODEs and plots Mp  %%% 
%%%                  Jacob Bellman, 1/22/2015                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

ep=0;     % Coupling factor of two oscillators
MDT=22;  % Mass doubling time of cell cycle

% Initialize.m specifies integration specifications
% (IntTime, Initial Values, options)
Initialize

% Solve the ODE system
[T,y] = ode45(@(t,y)JacobsCoupled(t,y,MDT,ep),[0 IntTime],...
              [W0 Fm0 Fp0 WFp0 Mp0 Ma0],options);

% Plot Mp
figure;
plot(T,y(:,5))
xlabel('t', 'FontSize', 14)
ylabel('M_p', 'FontSize', 14)
title({['M_p(t) Solution, ' ...
      'MDT = ', num2str(MDT), ...
      ', \epsilon = ', num2str(ep)]}, 'FontSize', 18)
