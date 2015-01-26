%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The following program solves JacobsCoupled  ODEs, finds %%% 
%%%  peaks of Mp, and plots the resulting histogram of the  %%%
%%%      distribution of cell cycle lengths for a given     %%%
%%%   coupling strength (ep) and mass doubling time (MDT)   %%%        
%%%                  Jacob Bellman, 1/22/2015               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

clear all

global MDT   % cell cycle uncoupled mass doubling time
global ep    % coupling factor

MDT=24; ep=15;

% options for solving the ode each time
options = odeset('RelTol',1e-5,'AbsTol',[1e-5 1e-5 1e-5 1e-5 1e-5 1e-5]);

% Initial Conditions for solving the ode each time that
% begin at Fm Max and Mp Max without coupling
W0=0.5296776; Fm0=0.958658; Fp0=0.7892854; WFp0=0.3232268;
Mp0=199.3611; Ma0=143.3817;

IntTime=10000;  % Integration Time

% Solve the ODE
[T,y] = ode45(@JacobsCoupled,[0 IntTime],[W0 Fm0 Fp0 WFp0 Mp0 Ma0],options);

% Find all peaks, locations of peaks, and number of peaks
[pks,locs]=findpeaks(y(:,5)); 
numpeaks = length(peaks);

% Find time distances between peaks
peakdiffs = zeros(numpeaks-1,1);
for i=1:numpeaks-1
    
    peakdiffs(i) = T(locs(i+1))-T(locs(i));
    
end

% Determine the width of the plotting frame
mindiffs = min(peakdiffs);
maxdiffs = max(peakdiffs);
q=(maxdiffs-mindiffs)/10;
histlow=mindiffs-q; histhigh=maxdiffs+q; % Smallest and largest values for hist
histsteps=20;   % Number of bars
histedges=linspace(histlow,histhigh,histsteps);  % Edges for the histogram bars

% Integration Time
IntTime=10000;

% Information for Histogram
[counts]=histc(peakdiffs,histedges);

% Plotting information
figure
h=bar(histedges,counts,'histc');
xlabel('Cell Cycle Length (h)','FontSize',14)
ylabel('Count','FontSize',14)
str1=sprintf('Histogram of Cell Cycle Lengths, \n MDT=%d, ep=%0.1f', MDT, ep);
title(str1,'FontSize',18,'FontWeight','bold')
str2=sprintf('examplefig%dby%d.fig',MDT,ep);
saveas(h,str2)

