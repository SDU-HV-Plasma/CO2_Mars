close all
clear all

%% INPUT PARAMETERS SECTION

P=4; %Torr
Tg=0.020;%eV

me=9.1e-31;%kg,电子质量

e=1.6e-19;%eV,kB

S_abs=3e7;%W/m3,初始能量沉积

%Parallel plate reactor
r=0.01;
area=pi*r^2;%m2
gap=0.002;%m

gas='He'; %Valid Options: He, He_N2, Ar, CR_Ar, CO2.
ne=1e17; %m-3
dt=1e-9; %sec
%% INITIALIZATION

%Geometry
volume=area*gap;%m3
area_volume=(2*area+2*pi*r*gap)/volume;%m-1


%Gas type
if (strcmp(lower(gas),lower('He_N2')))
    rates=@rates_he_n2;
    [nsp,species]=rates(0,P,[],area_volume,Tg,1);
    %Set initial relative densities
    n=10*[ 0.0100;
   0.0366;
   0.0050;
   0.0100;
   3.1040;
   0.0050;
   0.0047;
   0.0012;
   0.0001];
elseif (strcmp(lower(gas),lower('He')))
    rates=@rates_he3;
    mg=6.6969e-27;%kg，CO2分子质量；
    [nsp,species]=rates(0,P,[],area_volume,Tg,1,me,mg,e,S_abs,dt);
    %Set initial relative densities
    n=[
    1;
    10;
    0;
    0;
    10;
    1;
    1e7];
elseif (strcmp(lower(gas),lower('Ar')))
    rates=@rates_argon;
    [nsp,species]=rates(0,P,[],area_volume,Tg,1);
    %Set initial relative densities
    n=ones(nsp,1);
    
elseif (strcmp(lower(gas),lower('CO2')))
    rates=@rates_CO2;
    [nsp,species]=rates(0,P,[],area_volume,Tg,1);
    %Set initial relative densities
    n=ones(nsp,1);
    
elseif (strcmp(lower(gas),lower('CO2eg')))
    rates=@rates_CO2eg;
    [nsp,species]=rates(0,P,[],area_volume,Tg,1);
    %Set initial relative densities
    n=ones(nsp,1);
    
elseif (strcmp(lower(gas),lower('CO2m')))
    rates=@rates_CO2m;
    [nsp,species]=rates(0,P,[],area_volume,Tg,1);
    %Set initial relative densities
    n=[1;1;1;1;1;1;1;1;1000;1e7];
    
elseif (strcmp(lower(gas),lower('CO2m3')))
    rates=@rates_CO2m3;
    mg=7.307e-26;%kg，CO2分子质量；
    [nsp,species]=rates(0,P,[],area_volume,Tg,1,me,mg,e,S_abs,dt);
    %Set initial relative densities
    n=[1;1;1;1;1;1;1;1;1000;1e7];
    
elseif (strcmp(lower(gas),lower('CR_Ar')))
    rates=@cr_Ar;
    [nsp,species]=rates(0,P,[],area_volume,Tg,1);
    %Set initial relative densities
    n=ones(nsp,1)*ne;
    n(1)=3.54e13*P*1e9;%m-3    
    n(2:30)=0;
    n(31)=ne;
    n(32)=ne;
    n=n/ne;
    n(2:end)=[  
  2.5658342e+002
  2.0647860e-001
  4.3782594e+001
  5.5330014e-002
  3.2338558e-001
  5.8204280e-002
  7.3695097e-002
  1.5461018e-001
  3.0252550e-001
  1.1410669e-001
  2.0994363e-001
  1.3157183e-001
  2.0700317e-003
  1.8675975e-002
  4.5884439e-002
  9.4350130e-001
  1.2389130e-001
  3.9540199e-002
  1.2516645e-001
  4.4765584e-003
  1.5968685e-001
  1.8671126e-003
  8.9495283e-001
  1.6767913e-001
  2.6374535e-003
  8.9239697e-005
  2.6249563e-003
  6.3177122e+000
  6.1786283e+000
  1.0000000e+000
  1.0000000e+000];
n(1)=n(1)-sum(n(2:30));
else
    errordlg(sprintf('No model available for %s',gas),'Error');
    return;
end
% Find electron species
e_isp=find(strcmp(lower(species),'e')); 
if(isempty(e_isp))
e_isp=find(strcmp(lower(species),'electron')); 
end
if(isempty(e_isp))
e_isp=find(strcmp(lower(species),'electrons')); 
end
if(isempty(e_isp))
    errordlg('No electrons found in the species list. Make sure one spicies is called "e" or "electron"','ERROR');
    return;
end

MAX_HIST=100;
MAX_ITER=1e6;
MIN_ITER=1e3;

%Species
n=ne*n;%m-3
n_old=0;

%初始电子能量

Te=3;%eV

%history diagnostics
index=1;ctrl=1;jump=1;t=0;

%loop control
cont=0;

%% Create gui
fig=figure(1);
subplot(2,1,1)
xlabel('Time');
ylabel('Density');
subplot(2,1,2)
xlabel('Time');
ylabel('Te(eV)');
clf(fig);
scs=get(0,'ScreenSize');
width=.7*scs(3);
height=.8*scs(4);
set(fig,'position',[(scs(3)-width)/2 (scs(4)-height)/2 width height]);
%set(fig,'resize','off','windowstyle','modal','numbertitle','off');
set(fig,'resize','off');
bckcolor=get(fig,'color');

aux='All species';
for i=1:length(species)
    aux=[aux species(i)];
end
h_list=uicontrol(fig,'Style','popup','String',aux,'value',1,'tag','listofspecies',...
    'Position',[width/2-60,10,120,15]);

h_but=uicontrol(fig,'Style','pushbutton','String','Pause','userdata',1,'tag','pausebutton',...
    'Callback','start_pause','Position',[width/2-60-140,10,120,15]);

h_but_reaction=uicontrol(fig,'Style','pushbutton','String','Reaction Rates','userdata',1,'tag','reactionbutton',...
    'Callback','display_rates','Position',[width/2+60+20,10,120,15]);
set(fig,'Toolbar','figure');

h_but_end=uicontrol(fig,'Style','pushbutton','String','Finish','userdata',0,'tag','endbutton',...
    'Callback','finish','Position',[width/2+60+20+120,10,120,15]);
set(fig,'Toolbar','figure');

%% ITERATIONS

while(~get(h_but_end,'userdata'))
    waitfor(h_but,'string','Pause'); % Pause requested, wait for user 
    cont=cont+1;
    %Visualization
    if(rem(cont,1)==0 && cont>1)
        figure(1)
        subplot(2,1,1)
        selection=get(h_list,'value');
        if(selection>1)
            semilogy(hist_t,hist_n(:,selection-1),colour(selection-1));
            title(species(selection-1));
        else
            for i=1:length(species)
                semilogy(hist_t,hist_n(:,i),colour(i));
                hold on
            end
            title('All species');
            legend(gca,species,'location','EastOutside');
            hold off
        end

        xlabel('Time');
        ylabel('Density');
        subplot(2,1,2)
        plot(hist_t,hist_Te);
        xlabel('Time');
        ylabel('Te(eV)');

        data = guihandles(1); % initialize data to contain handles
        data.nsp = nsp;
        data.species = species;
        data.hist_t = hist_t;
        data.hist_n = hist_n;
        data.hist_Te= hist_Te;
        data.X = X;
        data.R = R;
        data.lambda=lambda;
        guidata(1, data);  % store the structure
        pause(0.2);
    end

    

    %Calculate reaction rates
    [nsp,species,X,R,lambda,dTe,B,Y]=rates(Te,P,n,area_volume,Tg,0,me,mg,e,S_abs,dt);
    
    t=t+dt;
    
    %Update Te
    Te=Te+dTe;
    
    if  Te<Tg 
        Te=Tg;
    end
    
    %Update density
    n_old=n;
    dn=dt*(sum(X')');
    ind=find(dn./n>2);
    dn(ind)=2*n(ind);
    ind=find(-dn./n>.5);
    dn(ind)=-.5*n(ind);
    
    n=n+dn;
    ind=find(n<0); n(ind)=1e-3;

    %Save histories
    if(ctrl==1)
        hist_n(index,:)=n';
        hist_t(index)=t;
        hist_Te(index)=Te;
        index=index+1;
        ctrl=jump;
        if(index>MAX_HIST)
            hist_t=hist_t(1:2:end);
            hist_Te=hist_Te(1:2:end);
            hist_n=hist_n(1:2:end,:);
            index=MAX_HIST/2+1;
            jump=jump*2;
        end
    else
        ctrl=ctrl-1;
    end
end

display_rates;
set(h_but_end,'enable','off');
set(h_but_reaction,'enable','off');
set(h_but,'enable','off');
set(h_list,'enable','off');
