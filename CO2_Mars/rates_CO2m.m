function [nsp,species,X,R,lambda]=rates_CO2m(Te,P,n,area_volume,Tg,only_e_balance)

% R1:  e + CO2 => CO2^+ + 2e       k= 5.4e-17 
% R2:  e + CO2 => CO + O + e       k= 5.8e-17
% R3:  e + CO2 => CO + O^-         k= 7.0e-18
% R4:  e + O3  => O + O2 + e       k= 2.0e-15
% R5:  e + O2  => 2O + e           k= 2.0e-15
% R6:  e + O2  => O + O^-          k= 4.0e-17
% R7:  e + O2 + M => O2^- + M      k= 3.0e-42
% R8:  O^- + CO => CO2 + e         k= 5.5e-16
% R9:  O^- + O2 => O3 + e          k= 1.0e-18
% R10: O^- + O3 => 2O2 + e         k= 3.0e-16
% R11: e + CO2^+ => CO + O         k= 6.5e-13
% R12: O2^- + CO2^+ => CO + O2 + O k= 6.0e-13
% R13: O + O + M => O2 + M         k= 5.2e-47*exp(900/Tg)
% R14: O + O2 + M => O3 + M        k= 4.5e-46*(Tg/298)^-2.70
% R15: O + O3 => 2O2               k= 8.0e-18*exp(-17.13/Tg)
% R16: O + CO + M => CO2 + M       k= 1.7e-45*exp(-1510/Tg)
% R17: O3 + M => O2 + O + M        k= 4.1e-16*exp(-11430/Tg)

species{1}='e';         charge(1)=-1;       mass(1)=9.1e-31;%kg 
species{2}='CO';        charge(2)=0;        mass(2)=4.65e-26;
species{3}='O_2';       charge(3)=0;        mass(3)=5.3e-26;
species{4}='O';         charge(4)=0;        mass(4)=2.66e-26;
species{5}='O_3';       charge(5)=0;        mass(5)=7.97e-26;
species{6}='O^-';       charge(6)=-1;       mass(6)=2.66e-26;
species{7}='O_2^-';     charge(7)=-1;       mass(7)=5.3e-26;
species{8}='CO_2^+';    charge(8)=1;        mass(8)=7.31e-26;
species{9}='M';         charge(9)=0;        
species{10}='CO_2';     charge(10)=0;       mass(10)=7.31e-26;
nsp=length(species);

if(Te==0)
    return;
end

ng=3.54e13*P*1e9;%m-3
 
%Gas phase reactions
number_of_reactions=17;
K(1:number_of_reactions)=0;
R(1:number_of_reactions)=0;
X(nsp,number_of_reactions+1)=0;
%First: Reactions involving generation/loss of electrons
K(1)=5.4e-17;
K(3)=7.0e-18;
K(6)=4.0e-17;
K(7)=3.0e-42;
K(8)=5.5e-16;
K(9)=1.0e-18;
K(10)=3.0e-16;
K(11)=6.5e-13;
R(1)=K(1)*n(1)*n(10);
R(3)=K(3)*n(1)*n(10);
R(6)=K(6)*n(1)*n(3);
R(7)=K(7)*n(1)*n(3)*n(9);
R(8)=K(8)*n(6)*n(2);
R(9)=K(9)*n(6)*n(3);
R(10)=K(10)*n(6)*n(5);
R(11)=K(11)*n(1)*n(8);
%Second: Rest of reactions
if(~only_e_balance)
    K(2)=5.8e-17;
    K(4)=2.0e-15;
    K(5)=2.0e-15;
    K(12)=6.0e-13;
    K(13)=5.2e-47*exp(0.07755/Tg);
    K(14)=4.5e-46*(Tg/0.02568)^(-2.70);
    K(15)=8.0e-18*exp(-0.00148/Tg);
    K(16)=1.7e-45*exp(-0.13011/Tg);
    K(17)=4.1e-16*exp(-0.98492/Tg);
    R(2)=K(2)*n(1)*n(10);
    R(4)=K(4)*n(1)*n(5);
    R(5)=K(5)*n(1)*n(3);
    R(12)=K(12)*n(7)*n(8);
    R(13)=K(13)*n(4)*n(9)*n(4);
    R(14)=K(14)*n(4)*n(3)*n(9);
    R(15)=K(15)*n(4)*n(5);
    R(16)=K(16)*n(4)*n(2)*n(9);
    R(17)=K(17)*n(5)*n(9);
end
%Wall losses
W(1)=0;
W(2)=0;
W(3)=0;
W(4)=0;
W(5)=0;
W(9)=0;
W(10)=0;
W(6)=-area_volume*n(6)*sqrt(1.6e-19*Te/mass(6));
W(7)=-area_volume*n(7)*sqrt(1.6e-19*Te/mass(7));
W(8)=-area_volume*n(8)*sqrt(1.6e-19*Te/mass(8));
W(1)=W(8)-W(6)-W(7);
% e: species 1
isp=1;
X(isp,1)=R(1);
X(isp,3)=-R(3);
X(isp,6)=-R(6);
X(isp,7)=-R(7);
X(isp,8)=R(8);
X(isp,9)=R(9);
X(isp,10)=R(10);
X(isp,11)=-R(11);

if(~only_e_balance)
    % CO: species 2
    isp=2;
    X(isp,2)=R(2);
    X(isp,3)=R(3);
    X(isp,8)=-R(8);
    X(isp,11)=R(11);
    X(isp,12)=R(12);
    X(isp,16)=-R(16);

    % O2: species 3
    isp=3;
    X(isp,4)=R(4);
    X(isp,5)=-R(5);
    X(isp,6)=-R(6);
    X(isp,7)=-R(7);
    X(isp,9)=-R(9);
    X(isp,10)=R(10)+R(10);
    X(isp,12)=R(12);
    X(isp,13)=R(13);
    X(isp,14)=-R(14);
    X(isp,15)=R(15)+R(15);
    X(isp,17)=R(17);
    
    % O: species 4
    isp=4;
    X(isp,2)=R(2);
    X(isp,4)=R(4);
    X(isp,5)=R(5)+R(5);
    X(isp,6)=R(6);
    X(isp,11)=R(11);
    X(isp,12)=R(12);
    X(isp,13)=-R(13)-R(13);
    X(isp,14)=-R(14);
    X(isp,15)=-R(15);
    X(isp,16)=-R(16);
    X(isp,17)=R(17);
    
    % O3: species 5
    isp=5;
    X(isp,4)=-R(4);
    X(isp,9)=R(9);
    X(isp,10)=-R(10);
    X(isp,14)=R(14);
    X(isp,15)=-R(15);
    X(isp,17)=-R(17);
    
     % O—: species 6
    isp=6;
    X(isp,3)=R(3);
    X(isp,6)=R(6);
    X(isp,8)=-R(8);
    X(isp,9)=-R(9);
    X(isp,10)=-R(10);
    
    % O2-: species 7
    isp=7;
    X(isp,7)=R(7);
    X(isp,12)=-R(12);
    
    % CO2+; species 8
    isp=8;
    X(isp,1)=R(1);
    X(isp,11)=-R(11);
    X(isp,12)=-R(12);
    
    % M; species 9
    isp=9;
    X(isp)=0;
    
   % CO2; species 10
    isp=10;
    X(isp,1)=-R(1);
    X(isp,2)=-R(2);
    X(isp,3)=-R(3);
    X(isp,8)=R(8);
    X(isp,16)=R(16);
    
end

%Wall losses

X(:,(length(R)+1))=W';

%Map reactions to wavelengths
%This can be use to generate synthetic spectra
lambda=zeros(1,length(R));
%Example: Reaction 4 radiates at 850nm
lambda(4)= 850; %此处非准确数值
    
end

