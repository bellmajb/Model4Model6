%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  The following program solves JacobsCoupled ODEs  %%% 
%%%   and calculates and displays the coefficient of  %%%
%%%   var for the last NP peaks (both time and value  %%%
%%%   distances) comparing consecutive peaks as well  %%%        
%%%                 as every other peak.              %%%
%%%               Jacob Bellman, 1/22/2015            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NP=20;  % Number of peaks considered when determining entrainment

% Evaluate the last NP peaks of Mp and time distances between them
[pks,locs]=findpeaks(y(:,5));  % Determine local maximums of Mp
TP=length(pks);                % total number of peaks found       

% Vectors for peak diffs
k=1:NP;       
Tdiff1 = T(locs(TP-NP+k))- ... % vector to hold time between peaks
         T(locs(TP-NP+k-1));  
Lastpeaks1 = y((locs(TP-NP+k)),5); % vector to hold peak values of Mp

% calc the coef of var for the last NP peaks and times between peaks
Tcv1=std(Tdiff1)/mean(Tdiff1);
Pcv1=std(Lastpeaks1)/mean(Lastpeaks1);
disp(['CV for last ', num2str(NP), ' time distances between peaks = ', ...
     num2str(abs(Tcv1))])
disp(['CV for last ', num2str(NP), ' peak vals of Mp = ' ...
     , num2str(abs(Pcv1))])

% Evaluate the last NP peaks of Mp (every other one)
% and time distances between them
Tdiff2=zeros(floor(NP/2),1);  % vector to hold time between every other peak
Lastpeaks2=zeros(floor(NP/2),1); % vector to hold every other peak value of Mp
k=1;
while 2*k<=NP     
    Tdiff2(k)=T(locs(TP-NP+2*k))-T(locs(TP-NP+2*k-2));
    Lastpeaks2(k)=y((locs(TP-NP+2*k)),5);
    k=k+1;
end

% calc the coef of var for every other one of the last NP peak values
% as well as the times between them
Tcv2=std(Tdiff2)/mean(Tdiff2);
Pcv2=std(Lastpeaks2)/mean(Lastpeaks2);
disp(['CV for time distances between every other peak = ', ...
     num2str(abs(Tcv2))])
disp(['CV for every other peak val of Mp = ', ...
     num2str(abs(Pcv2))])