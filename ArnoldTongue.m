%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  The following program solves JacobsCoupled ODEs for   %%%
%%%    various values of the coupling strength (ep) and    %%%
%%%  mass doubling time (MDT). It then plots a heat map    %%%
%%%  of the resulting Arnold Tongue for the coupled model. %%%
%%%     The period of the cell cycle is determined by      %%%
%%%    analyzing peaks of the inactive form of MPF (Mp).   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all   

MDT=24;  % cell cycle uncoupled mass doubling time
ep=0.5;  % coupling factor

Initialize

% Determines the grid of the Arnold Tongue plot
E1=0; E2=16;  % Smallest and largest values of ep for grid
M1=1; M2=50;  % Smallest and largest values of MDT for grid
Esteps=2; Msteps=2;          % Number of steps of ep and MDT for grid
Tthresh=0.01; Pthresh=0.01;  % Thresholds for periods and peaks to converge

% vectors containing all ep and MDT values
epvec=linspace(E1,E2,Esteps);
MDTvec=linspace(M1,M2,Msteps);

NP=20;   % Number of peaks considered

% Matrix for holding modal numbers (SynchMat(:,:,1))
% and periods (SynchMat(:,:,2)) associated with ep and MDT
SynchMat=zeros(Esteps,Msteps,2);

% Initialize progress bar
perc=0;
h=waitbar(perc/100, sprintf('%d%% done...',perc));

for i=1:Esteps
    ep=epvec(i);
    for j=1:Msteps
        MDT=MDTvec(j);
        [T,y] = ode45(@(t,y)JacobsCoupled(t,y,MDT,ep),[0 IntTime4Hist],...
                      [W0 Fm0 Fp0 WFp0 Mp0 Ma0],options);
        [pks,locs]=findpeaks(y(:,5));  % Determine local maximums of Mp
        TP=length(pks);                % Total number of peaks found
        
        % Make sure there are enough peaks to consider
        if (NP>TP)
            SynchMat(i,j,1)=0;    %solution is not periodic
            SynchMat(i,j,2)=0;
        else
            % See if the solution is periodic with one peak per period
            [mode1,period1]=unimodal(T,y,locs,TP,NP,Tthresh,Pthresh);
            if mode1>0.5
                SynchMat(i,j,1)=mode1;
                SynchMat(i,j,2)=period1;
            else
                % Solution may be bimodal; check:
                [mode2,period2]=bimodal(T,y,locs,TP,NP,Tthresh,Pthresh);
                if mode2>1.5
                    SynchMat(i,j,1)=mode2;
                    SynchMat(i,j,2)=period2;
                else
                    SynchMat(i,j,1)=0;  % solution is not periodic
                    SynchMat(i,j,2)=0;
                end
            end
        end
        
        % update progress bar
        pdone = (((i-1)*Esteps+j)/(Esteps*Msteps))*100; % percent done
        waitbar(pdone/100,h,sprintf('%0.1f%% done...',pdone))
    end
end

close(h)

% Plotting Information
figure
A=SynchMat(:,:,2);
h=image(MDTvec,epvec,A);
set(gca,'YDir','normal');
colormap('hot');
colorbar;
xlabel('Mass Doubling Time (MDT)')
ylabel('Coupling Strength (\epsilon)')
title('Arnold Tongue','FontSize',18,'FontWeight','bold')
t=colorbar('peer',gca);
set(get(t,'title'),'String','Entrained MDT')
str=sprintf(['ArnoldTongue' num2str(Esteps) 'by' num2str(Msteps) '.fig']);
saveas(h,str)