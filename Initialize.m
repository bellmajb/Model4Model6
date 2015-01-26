%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   The following program defines values necessary for the    %%%        
%%%           integration of JacobsCoupled Model                %%%
%%%                Jacob Bellman, 1/25/2015                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Options for solving the ode. The Relative Tolerance is required to
% be below 1e-4 for accurate results for the Arnold Tongue
options = odeset('RelTol',1e-5,'AbsTol',[1e-5 1e-5 1e-5 1e-5 1e-5 1e-5]);

% Initial Conditions for solving the ode that
% begin at Fm Max and Mp Max without coupling
W0=0.5296776; Fm0=0.958658; Fp0=0.7892854; WFp0=0.3232268;
Mp0=199.3611; Ma0=143.3817;

IntTime=2000;  % Integration Time