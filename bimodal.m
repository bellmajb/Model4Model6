%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Bimodal determines based on peaks and locations  %%%
%%%    of those peaks if a solution is periodic        %%% 
%%%                   and bimodal.                     %%%
%%%              Jacob Bellman, 2/1/2015               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mode,period] = bimodal(T,y,locs,TP,NP,Tthresh,Pthresh)
    
    Tdiff2=zeros(floor(NP/2),1);  % vector to hold time between every other peak
    Lastpeaks2=zeros(floor(NP/2),1); % vector to hold every other peak value of Mp
    k=1;
    while 2*k<=NP     
        Tdiff2(k)=T(locs(TP-NP+2*k))-T(locs(TP-NP+2*k-2));
        Lastpeaks2(k)=y((locs(TP-NP+2*k)),5);
        k=k+1;
    end
    % calc the coef of var for every other peak and peak time difference
    % of the final peaks
    Tcv2=std(Tdiff2)/mean(Tdiff2);
    Pcv2=std(Lastpeaks2)/mean(Lastpeaks2);
    if (abs(Tcv2)<Tthresh && abs(Pcv2)<Pthresh)
        mode=2;    % solution is periodic, bimodal
        % period is set as final peak diff (last and 2nd to last)
        period=T(locs(TP))-T(locs(TP-2));
    else
        mode=0;  % failed to be periodic with two peaks per period
        period=0;
    end
end

