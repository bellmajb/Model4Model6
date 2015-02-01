%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   The following program solves JacobsCoupled ODEs   %%% 
%%%       for a specified combination of coupling       %%%
%%%  strength (ep) and mass doubling time (MDT), finds  %%%   
%%%  peaks of Mp, and plots the resulting histogram of  %%%
%%%      the distribution of cell cycle lengths.        %%%        
%%%               Jacob Bellman, 1/22/2015              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

clear all   

MDT=24; % cell cycle uncoupled mass doubling time
ep=0.5;  % coupling factor

% Initialize.m specifies integration specifications
% (IntTime, Initial Values, options)
Initialize  

% Solve the ODE
[T,y] = ode45(@(t,y)JacobsCoupled(t,y,MDT,ep),[0 IntTime4Hist],...
              [W0 Fm0 Fp0 WFp0 Mp0 Ma0],options);

% Find the location of all peaks and their values
[pks,locs]=findpeaks(y(:,5)); 
numpeaks = length(peaks);

% Find time distances between peaks
peakdiffs = zeros(numpeaks-1,1);
ii=1:numpeaks-1;
peakdiffs(ii) = T(locs(ii+1))-T(locs(ii));

% Determine the width of the plotting frame
mindiffs = min(peakdiffs);
maxdiffs = max(peakdiffs);

% Take 10% of the median of all cell cycle lengths (q)
% and lengthen the plotting q above and below the max
% and mins, respectively.
q=(maxdiffs-mindiffs)/10;
histlow=mindiffs-q; histhigh=maxdiffs+q;

% Number of bars for the histogram
histsteps=20;                                    

% Edges for the histogram bars
histedges=linspace(histlow,histhigh,histsteps);  

% Calculate the histogram counts for bar heights
[counts]=histc(peakdiffs,histedges);

% Plotting information
figure
h=bar(histedges,counts,'histc');
xlabel('Cell Cycle Length (h)','FontSize',14)
ylabel('Count','FontSize',14)
str1=char({'Histogram of Cell Cycle Lengths', ...
          ['MDT = ' num2str(MDT) ', ep = ' num2str(ep)]});
title(str1,'FontSize',18,'FontWeight','bold')
str2=sprintf(['CellCycleHistMDT' num2str(MDT) 'ep' num2str(ep) '.fig']);
saveas(h,str2)
