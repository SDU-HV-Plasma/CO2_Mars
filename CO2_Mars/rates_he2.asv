function [nsp,species,X,R,lambda,dTe,B,Y]=rates_he2(Te,P,n,area_volume,Tg,only_e_balance,me,mg,e,S_abs,dt)

% R1:  e + He => He^* +e                K(1)=4.2e-9*Te^0.31*exp(-19.8/Te)*1e-6 
% R2:  e + He^* => e + He               K(2)=1.99e-10*Te^0.31*1e-6;
% R3:  e + He => 2e + He^+              K(3)=1.5e-9*Te^0.68*exp(-24.6/Te)*1e-6;
% R4:  e + He^* => 2e + He^+            K(4)=1.28e-7*Te^0.6*exp(-4.87/Te)*1e-6;
% R5:  e + He2^* => 2e + He2^+          K(5)=9.75e-10*Te^0.71*exp(-3.4/Te)*1e-6;
% R6:  e + He2^+ => He^* + He           K(6)=5.002e-9*Te^-.5*1e-6;
% R7:  2He^* => e + He^+ + He           K(7)=2.7e-10*1e-6;
% R8:  He^* + 2He => He2^* + He         K(8)=1.3e-34*1e-12;
% R9:  He^+ + 2He => He2^+ + He         K(9)=1e-31*1e-12; 
% R10: 2He2^* => e + He2^+ + 2He        K(10)=1.5e-9*1e-6;
% R11: e + He2^+ => 2He                 K(11)=1e-9*1e-6;
% R12: e + He2^+ + He => 3He            K(12)=2e-27*1e-12;
% R13: e + He => e + He                 K(13)=查表

%Species in present reaction set

species{1}='e';
species{2}='He^*';
species{3}='He^+';
species{4}='He_2';
species{5}='He_2^*';
species{6}='He_2^+';
species{7}='He';
nsp=length(species);

ng=3.54e13*P*1e9;%m-3

if(Te==0)
    return;
end

ng=3.54e13*P*1e9;%m-3

%Gas phase reactions
number_of_reactions=13;
K(1:number_of_reactions)=0;
R(1:number_of_reactions)=0;
X(nsp,number_of_reactions+1)=0;
%First: Reactions involving generation/loss of electrons
%should be changed accroding to the experimental observation.
K(3)=1.5e-9*Te^0.68*exp(-24.6/Te)*1e-6;  de(3)=24.6;          %K(3)=2.584e-12*Te^0.68*exp(-2.854092e5/Te)*1e-6;  de(3)=2.854092e5;
K(4)=1.28e-7*Te^0.6*exp(-4.87/Te)*1e-6;      de(4)=4.87; %K(4)=4.661e-10*Te^0.6*exp(-5.546e4/Te)*1e-6;      de(4)=5.546e4;
K(5)=9.75e-10*Te^0.71*exp(-3.4/Te)*1e-6;     de(5)=3.4; %K(5)=1.268e-12*Te^0.71*exp(-3.945e4/Te)*1e-6;     de(5)=3.945e4;
K(6)=5.002e-9*Te^-.5*1e-6;                        de(6)=0; %K(6)=5.386e-7*Te^-.5*1e-6;                        de(6)=0;
K(7)=2.7e-10*1e-6;                                de(7)=-15; %K(7)=2.7e-10*1e-6;                                de(7)=-15*11605;
K(10)=1.5e-9*1e-6;
K(11)=1e-9*1e-6;
K(12)=2e-27*1e-12;
R(3)=K(3)*n(1)*n(7);
R(4)=K(4)*n(1)*n(2);
R(5)=K(5)*n(1)*n(5);
R(6)=K(6)*n(1)*n(6);
R(7)=K(7)*n(2)*n(2);
R(10)=K(10)*n(5)*n(5);
R(11)=K(11)*n(1)*n(6);
R(12)=K(12)*n(1)*n(6)*n(7);
%Second: Rest of reactions
if(~only_e_balance)
    K(1)=4.2e-9*Te^0.31*exp(-19.8/Te)*1e-6;          de(1)=19.8; %K(1)=2.308e-10*Te^0.31*exp(-2.297e5/Te)*1e-6;     de(1)=2.297e5;
    K(2)=1.99e-10*Te^0.31*1e-6;                      de(2)=-19.8; %K(2)=1.099e-11*Te^0.31*1e-6;                      de(2)=-2.297e5;
    K(8)=1.3e-34*1e-12;                              de(8)=0; % Three body's reaction K(8)=1.3e-33*1e-6;                                de(8)=0;
    K(9)=1e-31*1e-12;                                de(9)=0; %K(9)=1e-31*1e-6;                                  de(9)=0;
    K(13)=lookups('reaction-rate-cross_section-e_H2elastic.lut',Te);%e与H2发生弹性碰撞的速率系数
    R(1)=K(1)*n(1)*n(7);
    R(2)=K(2)*n(1)*n(2);
    R(8)=K(8)*n(2)*n(7)*n(7);
    R(9)=K(9)*n(3)*n(7)*n(7);
end

%Wall losses
ub=sqrt(1.6e-19*Te/mg);%m/s % ub difference between electron, atom and molecule?

W(2)=-area_volume*0.25*n(2)*sqrt(8*1.6e-19*Tg/pi/mg);
W(3)=area_volume*n(3)*ub;
W(4)=0;
W(5)=-area_volume*0.25*n(5)*sqrt(8*1.6e-19*Tg/pi/mg/2);
W(6)=area_volume*n(6)*ub/sqrt(2);
W(1)=W(3)+W(6); 
W(7)=0;

% electrons
isp=1;
X(isp,3)=R(3);
X(isp,4)=R(4);
X(isp,5)=R(5);
X(isp,6)=-R(6);
X(isp,7)=R(7);
X(isp,10)=R(10);
X(isp,11)=-R(11);
X(isp,12)=-R(12);

if(~only_e_balance)
    % He*: species 2
    isp=2;
    X(isp,1)=R(1);
    X(isp,2)=-R(2);
    X(isp,4)=-R(4);
    X(isp,6)=R(6);
    X(isp,7)=-2*R(7);
    X(isp,8)=-R(8);

    % He+: species 3
    isp=3;
    X(isp,3)=R(3);
    X(isp,4)=R(4);
    X(isp,7)=R(7);
    X(isp,9)=-R(9);

    % He2: species 4
    isp=4;
    X(isp)=0;

    %He2*: species 5
    isp=5;
    X(isp,5)=-R(5);
    X(isp,8)=R(8);
    X(isp,10)=-2*R(10);

    % He2+: species 6
    isp=6;
    X(isp,5)=R(5);
    X(isp,6)=-R(6);
    X(isp,9)=R(9);
    X(isp,10)=R(10);
    X(isp,11)=-R(11);
    X(isp,12)=-R(12);

    % He: species 7
    isp=7;
    X(isp,1)=-R(1);
    X(isp,2)=R(2);
    X(isp,3)=-R(3);
    X(isp,6)=R(6);
    X(isp,7)=R(7);
    X(isp,8)=-R(8);
    X(isp,9)=-R(9);
    X(isp,10)=2*R(10);
    X(isp,11)=2*R(11);
    X(isp,12)=2*R(12);
end

% Wall Loss
X(:,(length(R)+1))=W';
%Map reactions to wavelengths
%This can be use to generate synthetic spectra
lambda=zeros(1,length(R));
%Example: Reaction 4 radiates at 850nm
lambda(4)= 850; %the value maybe is not correct here.

%能量变化
%ER1=19.82eV  ER3=24.59eV  ER4=4.7eV

dTe=2/(3*e*n(1))*S_abs*dt...
    -2*me/mg*K(13)*n(7)*(Te-Tg)*dt...
    -Te/n(1).*(sum(X(1,:)))*dt...
    -2/3*(19.82*K(1)*n(7)+24.59*K(3)*n(7)+4.7*K(4)*n(2))*dt...
    -10/(3*n(1))*Te*(-W(1))*dt;

%建立有向关系
%建立相互关系矩阵
i=1;
while i<nsp+1
    for j=1:nsp
        B{i,j}=intersect(find(X(i,:)>0),find(X(j,:)<0));
    end
    i=i+1;
end
%对该矩阵进行修正
B{2,1}=[1,B{2,1}];
B{3,1}=[3,4,B{3,1}];
B{6,1}=[5,B{6,1}];
B{7,1}=[2,B{7,1}];
%建立有向矩阵
i=1;
while i<nsp+1
    for j=1:nsp
        Y1(i,j)=sum(X(i,B{i,j}));
    end
    Y2(i)=sum(X(i,X(i,:)>0));
    i=i+1;
end
Y=Y1./Y2';

end

