%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Unimodal determines based on peaks and locations  %%%
%%%    of those peaks if a solution is periodic         %%% 
%%%                   and unimodal.                     %%%
%%%              Jacob Bellman, 2/1/2015                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mode,period] = unimodal(T,y,locs,TP,NP,Tthresh,Pthresh)
    
    Tdiff1=zeros(NP,1);  % vector to hold time between peaks
    Lastpeaks1=zeros(NP,1); % vector to hold peak values of Mp
    for k=1:NP     
        Tdiff1(k)=T(locs(TP-NP+k))-T(locs(TP-NP+k-1));
        Lastpeaks1(k)=y((locs(TP-NP+k)),5);        
    end
    % calc the coef of var for the last NP peaks and times between peaks
    Tcv=std(Tdiff1)/mean(Tdiff1);
    Pcv=std(Lastpeaks1)/mean(Lastpeaks1);
    if (abs(Tcv)<Tthresh && abs(Pcv)<Pthresh)
        mode=1;    % solution is periodic, unimodal
        period=T(locs(TP))-T(locs(TP-1));  % period is set as final peak diff
    else
        mode=0;  % failed to be periodic with one peak per period
        period=0;
    end
end

