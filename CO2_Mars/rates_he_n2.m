function [nsp,species,X,R,lambda]=rates_he(Te,P,n,area_volume,Tg,only_e_balance)

species{1}='e';
species{2}='He^*';
species{3}='He^+';
species{4}='He_2';
species{5}='He_2^*';
species{6}='He_2^+';
species{7}='N_2^+';
species{8}='N';
species{9}='N^+';
nsp=length(species);
Mhe=6.6969e-27;
MN=2.3e-26;


if(Te==0)
    return;
end

ng=3.54e13*P*1e9;%m-3
nN=ng*1e-7; %#ok<NOPRT>

%Gas phase reactions
number_of_reactions=20;
K(1:number_of_reactions)=0;
R(1:number_of_reactions)=0;
X(nsp,number_of_reactions+1)=0;
%First: Reactions involving generation/loss of electrons
K(3)=1.5e-9*Te^0.68*exp(-24.6/Te)*1e-6;  de(3)=24.6;          %K(3)=2.584e-12*Te^0.68*exp(-2.854092e5/Te)*1e-6;  de(3)=2.854092e5;
K(4)=1.28e-7*Te^0.6*exp(-4.87/Te)*1e-6;      de(4)=4.87; %K(4)=4.661e-10*Te^0.6*exp(-5.546e4/Te)*1e-6;      de(4)=5.546e4;
K(5)=9.75e-10*Te^0.71*exp(-3.4/Te)*1e-6;     de(5)=3.4; %K(5)=1.268e-12*Te^0.71*exp(-3.945e4/Te)*1e-6;     de(5)=3.945e4;
K(6)=5.002e-9*Te^-.5*1e-6;                        de(6)=0; %K(6)=5.386e-7*Te^-.5*1e-6;                        de(6)=0;
K(7)=2.9e-10*1e-6;                                de(7)=-15; %K(7)=2.7e-10*1e-6;                                de(7)=-15*11605;
K(10)=7.0*1e-11*1e-6;
K(11)=7.0*1e-11*1e-6;
K(16)=31.65e-30*Te^-0.8*1e-12;
K(17)=2.358e-8*Te^-0.5*1e-6;
K(19)=8.401e-5*exp(-14.5/Te)*1e-6; de(19)=14.5;
K(20)=2.71e-8*Te^-0.3*exp(-15.6/Te)*1e-6; de(20)=15.6;

R(3)=K(3)*n(1)*ng;
R(4)=K(4)*n(1)*n(2);
R(5)=K(5)*n(1)*n(5);
R(6)=K(6)*n(1)*n(6);
R(7)=K(7)*n(2)*n(2);
R(10)=K(10)*n(2)*n(8);
R(11)=K(11)*n(5)*n(8);
R(16)=K(16)*n(1)*n(1)*n(7);
R(17)=K(17)*n(1)*n(7);
R(19)=K(19)*n(1)*n(8);
R(20)=K(20)*n(1)*n(8);
%Second: Rest of reactions
if(~only_e_balance)
    K(1)=4.2e-9*Te^0.31*exp(-19.8/Te)*1e-6;     de(1)=19.8;       %K(1)=2.308e-10*Te^0.31*exp(-2.297e5/Te)*1e-6;     de(1)=2.297e5;
    K(2)=1.99e-10*Te^0.31*1e-6;                      de(2)=-19.8; %K(2)=1.099e-11*Te^0.31*1e-6;                      de(2)=-2.297e5;
    K(8)=1.66e-34*1e-12;                                de(8)=0; % Three body's reaction K(8)=1.3e-33*1e-6;                                de(8)=0;
    K(9)=1e-31*1e-12;                                  de(9)=0; %K(9)=1e-31*1e-6;                                  de(9)=0;
    K(12)=5.0*1e-10*1e-6;
    K(13)=7.0*1e-10*1e-6;
    K(14)=5.0*1e-10*1e-6;
    K(15)=7.0*1e-10*1e-6;
    K(18)=2.7975e-9*Te^-0.7*exp(-9.75/Te)*1e-6;  de(18)=9.757;

    R(1)=K(1)*n(1)*ng;
    R(2)=K(2)*n(1)*n(2);
    R(8)=K(8)*n(2)*ng^2;
    R(9)=K(9)*n(3)*ng^2;
    R(12)=K(12)*n(3)*n(8);
    R(13)=K(13)*n(3)*n(8);
    R(14)=K(14)*n(6)*n(8);
    R(15)=K(15)*n(6)*n(8);
    R(18)=K(18)*n(1)*n(8);
end
%Wall losses
ub=sqrt(1.6e-19*Te/Mhe);%m/s % ub difference between electron, atom and molecule?
uN=sqrt(1.6e-19*Te/MN);

W(2)=-area_volume*0.25*n(2)*sqrt(8*1.6e-19*Tg/pi/Mhe);
W(3)=-area_volume*n(3)*ub;
W(4)=0;
W(5)=-area_volume*0.25*n(5)*sqrt(8*1.6e-19*Tg/pi/Mhe/2);
W(6)=-area_volume*n(6)*ub/sqrt(2);
W(7)=-area_volume*n(7)*uN/sqrt(2);
W(8)=-area_volume*0.25*n(8)*sqrt(8*1.6e-19*Tg/pi/MN);
W(9)=-area_volume*n(9)*uN;
W(1)=(W(3)+W(6)+W(7)+W(9)); %??

% electrons
isp=1;
X(isp,3)=R(3);
X(isp,4)=R(4);
X(isp,5)=R(5);
X(isp,6)=-R(6);
X(isp,7)=R(7);
X(isp,10)=R(10);
X(isp,11)=R(11);
X(isp,16)=-R(16);
X(isp,17)=-R(17);
X(isp,19)=R(19);
X(isp,20)=R(20);

if(~only_e_balance)
% He*: species 2
isp=2;
X(isp,1)=R(1);
X(isp,2)=-R(2);
X(isp,4)=-R(4);
X(isp,6)=R(6);
X(isp,7)=-R(7);
X(isp,8)=-R(8);
X(isp,10)=-R(10);

% He+: species 3
isp=3;
X(isp,3)=R(3);
X(isp,4)=R(4);
X(isp,7)=R(7);
X(isp,9)=-R(9);
X(isp,12)=-R(12);
X(isp,13)=-R(13);

% He2: species 4
isp=4;
X(isp)=0;

%He2*: species 5
isp=5;
X(isp,5)=-R(5);
X(isp,8)=R(8);
X(isp,11)=-R(11);

% He2+: species 6
isp=6;
X(isp,5)=R(5);
X(isp,6)=-R(6);
X(isp,9)=R(9);
X(isp,14)=-R(14);
X(isp,15)=-R(15);

% N2+: species 7
isp=7;
X(isp,10)=R(10);
X(isp,11)=R(11);
X(isp,12)=R(12);
X(isp,14)=R(14);
X(isp,16)=-R(16);
X(isp,17)=-R(17);
X(isp,20)=R(20);

% N: species 8
isp=8;
X(isp,13)=R(13);
X(isp,15)=R(15);
X(isp,17)=2*R(17);
X(isp,18)=2*R(18);
X(isp,19)=-R(19);

% N+: species 9
isp=9;
X(isp,13)=R(13);
X(isp,15)=R(15);
X(isp,19)=R(19);
end

% Wall Loss
X(:,(length(R)+1))=W(:);
%This can be use to generate synthetic spectra
lambda=zeros(1,length(R));
%Example: Reaction 4 radiates at 850nm
lambda(4)= 850; %the value maybe is not correct here.