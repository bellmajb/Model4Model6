%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The following set of ODEs simulates 4 molecular components %%%
%%%    of the Circadian clock in Neurospora crassa, and 2      %%%
%%%       molecular components of the division cycle.          %%%
%%%                                                            %%%
%%%          ODE system for the circadian clock:               %%%
%%%         Fm / dt = k3W^n / ( K^n + W^n ) - k4Fm             %%%
%%%         Fp / dt = k5Fm - k6Fp + k7WFp - k8WFp              %%%
%%%           WFp / dt = k8WFp - k7WFp - k9WFp                 %%%
%%%          W / dt = k1 - k2W + k7WFp - k8WFp                 %%%
%%%                                                            %%%
%%%          ODE system for the circadian clock:               %%%
%%%  Mp / dt = ( P / MDT )[ vp - c1Mp - c2Mp - ...             %%%
%%%             c3MpMa^s / ( Km^s + Ma^s ) +c5Ma( L + epsW ) ] %%%
%%%  Ma / dt = ( P / MDT )[ c2Mp - c4Ma + ...                  %%%
%%%            c3MpMa^s / ( Km^s + Ma^s ) - c5Ma( L + epsW ) ] %%%                                                   %%%
%%%                Jacob Bellman, 1/22/2015                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dy = JacobsCoupled(t,y,MDT,ep)

% The following parameter set causes the clock to oscillate
% robustly with a period of about P=22 hours
k1 = 0.1619;  % Constant synthesis of WC-1
k2 = 0.0121;  % Degradation of WC-1
k3 = 5.8968;  % Maximum transcription rate of frq
k4 = 0.1226;  % Degradation of frq mRNA
k5 = 0.3206;  % Translation of FRQ
k6 = 0.0402;  % Degradation of FRQ
k7=0.0175;    % Dissociation of WC-1:FRQ complex
k8 = 0.4740;  % Association of WC-1:FRQ complex
k9 = 0.5256;  % Degradation of WC-1:FRQ complex
m = 7;        % Hill coefficient of transcription
Kw = 0.924;   % Threshold for transcription of frq 
P = 22;       % Period of the circadian clock

% The following parameter set causes the cell cycle to oscillate
% with a mass doubling time of MDT
vp = 48.06;   % Constant synthesis of Mp
c1 = 0.04005; % Degradation of Mp
c2 = 0.22695; % Rate of Mp becoming active and turning into Ma
c3 = 1.068;   % Max rate of autocatalysis of Ma
c4 = 0.2403;  % Degradation of Ma
c5 = 0.1068;  % Rate of Ma becoming inactive and turning into Mp
n = 9;        % Hill coefficient of autocatalysis of Ma
Km = 200;     % Threshold for autocatalysis of Ma
L = 1;        % Linear rate of coupling

dy = zeros(6,1);    % initialize ODE RHS

% Circadian Clock Equations
dy(1) = k1 - k2*y(1) + k7*y(4) - k8*y(1)*y(3);       %WC-1
dy(2) = k3*y(1)^m / ( Kw^m + y(1)^m ) - k4*y(2);     %frq mRNA
dy(3) = k5*y(2) - k6*y(3) + k7*y(4) - k8*y(1)*y(3);  %FRQ protein
dy(4) = k8*y(1)*y(3) - k7*y(4) - k9*y(4);            %WC-1:FRQ complex

% Cell Cycle Equations
dy(5) = ( P / MDT )*( vp - c1*y(5) - c2*y(5) - ( c3*y(5)*y(6)^n ) / ... %Mp
        ( Km^n + y(6)^n ) + c5*y(6)*( L + ep*y(1) ) );  
dy(6) = ( P / MDT )*( c2*y(5) - c4*y(6) + ( c3*y(5)*y(6)^n ) / ...      %Ma
        ( Km^n + y(6)^n ) - c5*y(6)*( L + ep*y(1) ) );