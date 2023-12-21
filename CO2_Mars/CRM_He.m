close all
clear all
%% INPUT PARAMETERS SECTION

global Ratesumexini nsp species sim nall dt closeness betweenness B number_of_reactions X Rall

tic
Tg=0.02586;%eV

me=9.1e-31;%kg,电子质量

e=1.6e-19;%eV,kB

S_absm=1e6;%W/m3,能量密度

%判断是否纳入振动态
Vibration=1;

%初始振动态温度
Tv=Tg;

%判断是否进行简化
sim=0;
%通过简化除去的粒子种类
if sim==1
    spc_sim=10;
else
    spc_sim=0;
end

%Parallel plate reactor
r=0.01;
area=pi*r^2;%m2
gap=0.002;%m

gas='HeO'; %Valid Options: HeO
ne=1e14;%m-3
dt=1e-11; %sec 
%% INITIALIZATION

%Geometry
volume=area*gap;%m3
area_volume=(2*area+2*pi*r*gap)/volume;%m-1

%Gas type

if (strcmp(lower(gas),lower('HeO')))   %无M的参与
    rates=@rates_HeO;
    matrix=@HeO_matrix;
    %Set initial relative densities
    [nsp,species,K,number_of_reactions,mass,aux,auxv,Kv,number_of_reaction,number_of_vreaction]=rates_HeO0(Tg,Vibration);
    if Vibration==1
        n=[1;2.15e11;1;2.5e10;1;1;1;1;1;1;1;1];
        nv=[1e-8;1e-8;1e-8;1e-8;1e-8];%%%%%%%%%
    else
        n=[1;2.15e11;1;2.5e10;1;1;1;1;1;1;1;1];
        nv=[0;0;0;0;0];
    end
    if sim==1
        n=[];
        nv=[];
    end
        
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

MAX_HIST=1000;
% MAX_ITER=1e6;
% MIN_ITER=1e3;

%Species
n=ne*n;
nv=ne*nv;%m-3
n_old=0;
nv_old=0;

%初始电子能量

Te=0.5;

%history diagnostics
index=1;ctrl=1;jump=1;t=0;

%loop control
cont=0;
contchange=5000;

%determine random line color
% c=numel(n);
randcolor=rand(nsp,3);
    
%% Create gui
fig=figure(1);
subplot(3,1,1)
xlabel('Time');
ylabel('Density');
subplot(3,1,2)
xlabel('Time');
ylabel('Te(eV)');
subplot(3,1,3)
xlabel('Time');
ylabel('单位时间内的能量情况（eV）');
clf(fig);
scs=get(0,'ScreenSize');
width=.9*scs(3);
height=.8*scs(4);
set(fig,'position',[(scs(3)-width)/2 (scs(4)-height)/2 width height]);
%set(fig,'resize','off','windowstyle','modal','numbertitle','off');
set(fig,'resize','off');
bckcolor=get(fig,'color');

% auxx='All species';
% for i=1:length(species)
%     auxx=[auxx species(i)];
% end
auxx=species;

h_list=uicontrol(fig,'Style','popup','String',auxx,'value',1,'tag','listofspecies',...
    'Position',[width/2-60,10,120,15]);

h_but=uicontrol(fig,'Style','pushbutton','String','Pause','userdata',1,'tag','pausebutton',...
    'Callback','start_pause','Position',[width/2-60-140,10,120,15]);

h_but_reaction=uicontrol(fig,'Style','pushbutton','String','Reaction Rates','userdata',1,'tag','reactionbutton',...
    'Callback','display_rates(nsp,species,number_of_reactions)','Position',[width/2+60+20,10,120,15]);
set(fig,'Toolbar','figure');

h_but_end=uicontrol(fig,'Style','pushbutton','String','Finish','userdata',0,'tag','endbutton',...
    'Callback','finish','Position',[width/2+60+20+120,10,120,15]);
set(fig,'Toolbar','figure');

%% ITERATIONS

hist_n=zeros(MAX_HIST,length(n));
hist_nv=zeros(MAX_HIST,length(nv));
hist_nall=zeros(MAX_HIST,nsp);
hist_R=zeros(MAX_HIST,number_of_reaction);
hist_Rv=zeros(MAX_HIST,number_of_vreaction);
hist_K=zeros(MAX_HIST,number_of_reaction);
hist_Kv=zeros(MAX_HIST,number_of_vreaction);
hist_t=zeros(1,MAX_HIST);
hist_Te=zeros(1,MAX_HIST);
hist_Tv=zeros(1,MAX_HIST);
hist_dS_ape=zeros(1,MAX_HIST);
hist_dS_el=zeros(1,MAX_HIST);
hist_dS_epro=zeros(1,MAX_HIST);
hist_dS_rea=zeros(1,MAX_HIST);
hist_dS_wall=zeros(1,MAX_HIST);

while(~get(h_but_end,'userdata'))
    waitfor(h_but,'string','Pause'); % Pause requested, wait for user 
    cont=cont+1;
% if cont<=3000 || ( cont>=100000 && cont<=103000 ) || ( cont>=200000 && cont<=203000 ) || ( cont>=300000 && cont<=303000 ) || ( cont>=400000 && cont<=403000 )
%     S_abs=S_absm;
% else
%     S_abs=0;
% end
%         step=fix(cont/100);
%     if mod(fix(step),5)==0
%         S_abs=S_absm;
%     else
%         S_abs=0;
%     end

%     S_abs=S_absm*(1-heaviside(cont-contchange));
%       S_abs=S_absm*(heaviside(cont)-heaviside(cont-contchange)+heaviside(cont-1e3)-heaviside(cont-1e3-contchange)+heaviside(cont-2e3)-heaviside(cont-2e3-contchange)...
%           +heaviside(cont-3e3)-heaviside(cont-3e3-contchange)+heaviside(cont-4e3)-heaviside(cont-4e3-contchange)+heaviside(cont-5e3)-heaviside(cont-5e3-contchange)...
%           +heaviside(cont-6e3)-heaviside(cont-6e3-contchange)+heaviside(cont-7e3)-heaviside(cont-7e3-contchange)+heaviside(cont-8e3)-heaviside(cont-8e3-contchange)...
%           +heaviside(cont-9e3)-heaviside(cont-9e3-contchange)+heaviside(cont-10e3)-heaviside(cont-10e3-contchange));
%             S_abs=S_absm*(heaviside(cont)-heaviside(cont-contchange)+heaviside(cont-1e4)-heaviside(cont-1e4-contchange)+heaviside(cont-2e4)-heaviside(cont-2e4-contchange)...
%           +heaviside(cont-3e4)-heaviside(cont-3e4-contchange)+heaviside(cont-4e4)-heaviside(cont-4e4-contchange));
%       S_abs=S_absm*~mod(fix(cont/contchange),2);
  S_abs=S_absm;
    %Visualization
    if(rem(cont,2)==0 && cont>1)
        figure(1)
        subplot(3,1,1)
        selection=get(h_list,'value');
        semilogy(hist_t,hist_nall(:,selection),'color',randcolor(selection,:),'linewidth',1);
        title(species(selection));
%         if(selection>1)
%             semilogy(hist_t,hist_n(:,selection-1),'color',randcolor(selection-1,:),'linewidth',1);
%             title(species(selection-1));
%         else
%             for i=1:length(species)
%                 semilogy(hist_t,hist_n(:,i),'color',randcolor(i,:),'linewidth',1);
%                 hold on
%             end
%             title('All species');
%             %legend({'e','CO','O_2','O','O_3','O^-','O_2^-','CO_2^+','M','CO_2','C','O_2^+','CO_3^-'},'location','northeast');
%             hold off;
%         end

        xlabel('Time');
        ylabel('Density');
        
       % hist_Te(hist_Te>10)=10;%为使中后期图像便于观察变化趋势，施加强制性代码
        subplot(3,1,2)
        plot(hist_t,hist_Te,'linewidth',1);
        xlabel('Time');
        ylabel('Te(eV)');   

        subplot(3,1,3)
        plot(hist_t,hist_dS_ape,'r',hist_t,hist_dS_el,'g',hist_t,hist_dS_epro,'b',hist_t,hist_dS_rea,'m',hist_t,hist_dS_wall,'k');
        legend('单位电子能量沉积','弹性碰撞能量散失','电子密度变化引起的电子温度变化','反应化学能带来的能量变化','物质在壁上的损失带来的能量变化');
        xlabel('Time');
        ylabel('单位时间内的能量情况（eV）');

        data = guihandles(1); % initialize data to contain handles
        data.nsp = nsp;
        %data.K = K;
        %data.Kv = Kv;
        data.species = species;
        data.hist_t = hist_t;
        data.hist_n = hist_n;
        %data.hist_nv = hist_nv;
        data.hist_nall = hist_nall;
        %data.hist_Rate = hist_Rate;
        %data.hist_RateK = hist_RateK;
        data.hist_Te= hist_Te;
        data.hist_Tv= hist_Tv;
        data.hist_dS_ape= hist_dS_ape;
        data.hist_dS_el= hist_dS_el;
        data.hist_dS_epro= hist_dS_epro;
        data.hist_dS_rea= hist_dS_rea;
        data.hist_dS_wall= hist_dS_wall;
        data.X = X;
        %data.R = R;
        %data.Rv = Rv;
        guidata(1, data);  % store the structure
        pause(0.2);
    end
    
    %Calculate reaction rates
    
    [X,R,dTe,dS_ape,dS_el,dS_epro,dS_rea,dS_wall,K0,Rv,Kv0,Ralldlt]=rates(Te,n,area_volume,Tg,me,e,S_abs,dt,Vibration,K,nsp,number_of_reactions,0,mass,aux,auxv,Kv,number_of_reaction,number_of_vreaction,nv,sim);
    
    t=t+dt;
    
    K=K0;
    Kv=Kv0;
    
    %设置运行时间

    if t>=1e-3
        break
    end

    %Update Te
%   if(dTe/Te>0.2)
%        dTe=0.2*Te;
%   end
%    if(dTe/Te<-0.2)
%        dTe=-0.2*Te;
%    end
    
    Te=Te+dTe;
    Te(Te<=Tg)=Tg;
%     Te(Te>=49.95)=49.95;%供插值使用
    
    %Update density
    n_old=n;
    nv_old=nv;
    dn=dt*(sum(X')');
%     ind=find(dn./n>5);
%     dn(ind)=5*n(ind);
%     ind=find(-dn./n>.5);
%     dn(ind)=-.5*n(ind);
    
    n=n+dn(1:length(n));
    nv=nv+dn(length(n)+1:nsp);
    n(n<=1e-3)=1e-3;
    nv(nv<1e-3)=1e-3;
    nall=[n;nv];
    
    %Update Tv

%     Tv=Vibration*(-0.29/log(nv(7)/n(10)));%Tv=E1/(Kb*ln(nv1/nv0))这是K值
    Tv=Vibration*(-0.29/log(sum(nv)/(n(2)+n(4))));
    Tv(Tv<=Tg)=Tg;
    
%     if t>1e-3-5e-5-2*dt && t<1e-3-5e-5
%         index=1;ctrl=1;jump=1;
%         hist_n=zeros(MAX_HIST,length(n));
%         hist_nv=zeros(MAX_HIST,length(nv));
%         hist_nall=zeros(MAX_HIST,nsp);
%         hist_R=zeros(MAX_HIST,number_of_reaction);
%         hist_Rv=zeros(MAX_HIST,number_of_vreaction);
%         hist_K=zeros(MAX_HIST,number_of_reaction);
%         hist_Kv=zeros(MAX_HIST,number_of_vreaction);
%         hist_t=zeros(1,MAX_HIST);
%         hist_Te=zeros(1,MAX_HIST);
%         hist_Tv=zeros(1,MAX_HIST);
%         hist_dS_ape=zeros(1,MAX_HIST);
%         hist_dS_el=zeros(1,MAX_HIST);
%         hist_dS_epro=zeros(1,MAX_HIST);
%         hist_dS_rea=zeros(1,MAX_HIST);
%         hist_dS_wall=zeros(1,MAX_HIST);
%     end
    
    %Save histories
    if(ctrl==1)
        hist_n(index,:)=n';
        hist_nv(index,:)=nv';
        hist_nall(index,:)=nall';
        hist_K(index,:)=K;
        hist_Kv(index,:)=Kv;
        hist_R(index,:)=R;
        hist_Rv(index,:)=Rv;
        hist_t(index)=t;
        hist_Te(index)=Te;
        hist_Tv(index)=Tv;
        hist_dS_ape(index)=dS_ape;
        hist_dS_el(index)=dS_el;
        hist_dS_epro(index)=dS_epro;
        hist_dS_rea(index)=dS_rea;
        hist_dS_wall(index)=dS_wall;
        index=index+1;
        ctrl=jump;
        if(index>MAX_HIST)
            hist_t=hist_t(1:2:end);
            hist_Te=hist_Te(1:2:end);
            hist_Tv=hist_Tv(1:2:end);
            hist_dS_ape=hist_dS_ape(1:2:end);
            hist_dS_el=hist_dS_el(1:2:end);
            hist_dS_epro=hist_dS_epro(1:2:end);
            hist_dS_rea=hist_dS_rea(1:2:end);
            hist_dS_wall=hist_dS_wall(1:2:end);
            hist_n=hist_n(1:2:end,:);
            hist_nv=hist_nv(1:2:end,:);
            hist_K=hist_K(1:2:end,:);
            hist_Kv=hist_Kv(1:2:end,:);
            hist_R=hist_R(1:2:end,:);
            hist_Rv=hist_Rv(1:2:end,:);
            index=MAX_HIST/2+1;
            hist_t(index:MAX_HIST)=0;
            hist_Te(index:MAX_HIST)=0;
            hist_Tv(index:MAX_HIST)=0;
            hist_dS_ape(index:MAX_HIST)=0;
            hist_dS_el(index:MAX_HIST)=0;
            hist_dS_epro(index:MAX_HIST)=0;
            hist_dS_rea(index:MAX_HIST)=0;
            hist_dS_wall(index:MAX_HIST)=0;
            hist_n(index:MAX_HIST,:)=0;
            hist_nv(index:MAX_HIST,:)=0;
            hist_K(index:MAX_HIST,:)=0;
            hist_Kv(index:MAX_HIST,:)=0;
            hist_R(index:MAX_HIST,:)=0;
            hist_Rv(index:MAX_HIST,:)=0;
            hist_nall=[hist_n hist_nv];
            jump=jump*2;
        end

    else
        ctrl=ctrl-1;
    end
    
end

[B,G,Y,Bsim,Gsim]=matrix(X,nsp,number_of_reactions,Ralldlt,sim,Vibration,n,number_of_reaction);%G为拓扑图  BG列损行生
[clo,num_of_clo]=doclose(Gsim,nsp,spc_sim);%clo表示拓扑图中粒子间最短路径的长度,num_of_clo表示有多少种这样的最短路径
closeness=(nsp-spc_sim-1)./sum(clo);%临近中心性
[Bet]=dobetween(clo,nsp,num_of_clo,spc_sim);%Bet表示每两个粒子间的最短路径中，经过待求粒子的最短路径数目
for m=1:nsp-spc_sim
    BET{m,1}=Bet{m,1}./(num_of_clo+eye(nsp-spc_sim));%计算相对值，过程中使num_of_clo对角线的0部分强制置1
    betweenness(m,1)=sum(sum(BET{m,1}))./((nsp-spc_sim-1)*(nsp-spc_sim-2)/2);%相间中心性
end
[clustering]=docluster(nsp,clo,spc_sim);%紧密系数
[reaction]=checkreaction(Gsim,nsp,spc_sim);%各粒子的参与反应(拓扑分析)
Rall=[R Rv];Ratesum=zeros(nsp,nsp);Ratemax=zeros(nsp,nsp);
for i=1:nsp
    for j=1:nsp
        Ratesum(i,j)=sum(X(i,B{i,j}));
        %Ratesum(i,j)=sum(Rall(B{i,j}));
        if ~isempty(max(Rall(B{i,j})))
            Ratemax(i,j)=max(X(i,B{i,j}));
            Ratemax_reaction{i,j}=find(X(i,:)==Ratemax(i,j));
            %Ratemax(i,j)=max(Rall(B{i,j}));
        end
    end
end
degree=0;
for i=1:nsp
    for j=1:nsp
        degree=degree+length(G{i,j});
    end
    indegree(i)=degree;
    degree=0;
end
for i=1:nsp
    for j=1:nsp
        degree=degree+length(G{j,i});
    end
    outdegree(i)=degree;
    degree=0;
end
degree=indegree+outdegree;
Ratesumexini=Ratesum*1/10*dt./nall';%Crutial!! 处理dt值保证各边均大于0
%display_rates(nsp,species,number_of_reactions);
set(h_but_end,'enable','off');
set(h_but_reaction,'enable','off');
set(h_but,'enable','off');
set(h_list,'enable','off');
toc