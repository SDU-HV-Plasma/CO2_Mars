function [nsp,species,X,R,lambda] = rates_CO2eg(Te,P,n,area_volume,Tg,only_e_balance)

% R1:   e + CO_2   => 2e   + CO_2^+        k = 查表
% R2:   O + CO_2   => O_2  + CO            k = 2.80e-11*exp(-26500/Tg)
% R3: O_2 + CO     => CO_2 + O             k = 4.20e-12*exp(-24000/Tg)
% R4:   e + CO_2^+ => CO_2                 k = 2.15e-03*Te*Te*(-0.50)/Tg 

species{1}='e';         charge(1)=-1;       mass(1)=9.1e-31;%kg 
species{2}='CO';        charge(2)=0;        mass(2)=4.65e-26;
species{3}='O_2';       charge(3)=0;        mass(3)=5.3e-26;
species{4}='O';         charge(4)=0;        mass(4)=2.66e-26;
species{5}='CO_2^+';    charge(5)=1;        mass(5)=7.31e-26;
nsp=length(species);

if(Te==0)
    return;
end

ng=3.54e13*P*1e9;%m-3

%Gas phase reactions
number_of_reactions=4;
K(1:number_of_reactions)=0;
R(1:number_of_reactions)=0;
X(nsp,number_of_reactions+1)=0;
%First: Reactions involving generation/loss of electrons
K(1)=lookup('rate-coefficients_R8.lut',Te);
K(4)=2.15e-03*Te*Te*(-0.50)/Tg;
R(1)=K(1)*n(1)*ng;
R(4)=K(4)*n(1)*n(5);
%Second: Rest of reactions
if(~only_e_balance)
    K(2)=2.80e-11*exp(-26500/Tg);
    K(3)=4.20e-12*exp(-24000/Tg);
    R(2)=K(2)*n(4)*ng;
    R(3)=K(3)*n(3)*n(2);
end

%Wall losses
for loop=2:4
    W(loop)=-area_volume*0.25*n(loop)*sqrt(8*1.6e-19*Tg/pi/mass(loop));
end
W(5)=-area_volume*n(5)*sqrt(1.6e-19*Te/mass(5));
W(1)=W(5);

% e: species 1
isp=1;
X(isp,1)=R(1);
X(isp,4)=-R(4);

if(~only_e_balance)
    % CO: species 2
    isp=2;
    X(isp,2)=R(2);
    X(isp,3)=-R(3);
    
     % O2: species 3
    isp=3;
    X(isp,2)=R(2);
    X(isp,3)=-R(3);
    
    % O: species 4
    isp=4;
    X(isp,2)=-R(2);
    X(isp,3)=R(3);

    % CO2+; species 5
    isp=5;
    X(isp,1)=R(1);
    X(isp,4)=-R(4);
end

%Wall losses
X(:,(length(R)+1))=W(:);

%Map reactions to wavelengths
%This can be use to generate synthetic spectra
lambda=zeros(1,length(R));
%Example: Reaction 4 radiates at 850nm
lambda(4)= 850; %此处非准确数值

end

