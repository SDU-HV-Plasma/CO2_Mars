function [nsp,species,X,R,lambda]=rates_CO2(Te,P,n,area_volume,Tg,only_e_balance)

% R1:         e  + CO_2   => O^-    + CO            k = 查表
% R2:         CO + O^-    => CO_2   + e             k = 5.5e-16
% R3:         CO + CO_3^- => 2CO_2  + e             k = 5e-19
% R4: CO_2 + O^- + CO     => CO_3^- + CO            k = 1.5e-40
% R5:     CO_3^- + O      => CO_2   + O_2^-         k = 8e-17
% R6:          O + O_2^-  => O_2    + O^-           k = 1.5e-16
% R7:          e + O_2    => O      + O^-           k = 查表
% R8:          e + CO_2   => 2e     + CO_2^+        k = 查表
% R9:     CO_3^- + CO_2^+ => 2CO_2  + O             k = 5e-13

species{1}='e';         charge(1)=-1;       mass(1)=9.1e-31;%kg 
species{2}='CO';        charge(2)=0;        mass(2)=4.65e-26;
species{3}='O_2';       charge(3)=0;        mass(3)=5.3e-26;
species{4}='O';         charge(4)=0;        mass(4)=2.66e-26;
species{5}='O^-';       charge(5)=-1;       mass(5)=2.66e-26;
species{6}='O_2^-';     charge(6)=-1;       mass(6)=5.3e-26;
species{7}='CO_3^-';    charge(7)=-1;       mass(7)=10e-26;
species{8}='CO_2^+';    charge(8)=1;        mass(8)=7.31e-26;
nsp=length(species);

if(Te==0)
    return;
end

ng=3.54e13*P*1e9;%m-3

%Gas phase reactions
number_of_reactions=8;
K(1:number_of_reactions)=0;
R(1:number_of_reactions)=0;
X(nsp,number_of_reactions+1)=0;
%First: Reactions involving generation/loss of electrons  
K(1)=lookups('rate-coefficients_R1.lut',Te);
K(2)=5.5e-16;      
K(3)=5e-19;        
K(7)=lookups('rate-coefficients_R7.lut',Te);
K(8)=lookups('rate-coefficients_R8.lut',Te);
R(1)=K(1)*n(1)*ng;
R(2)=K(2)*n(2)*n(5);
R(3)=K(3)*n(2)*n(7);
R(7)=K(7)*n(1)*n(3);
R(8)=K(8)*n(1)*ng;
%Second: Rest of reactions
if(~only_e_balance)
    K(4)=1.5e-40; 
    K(5)=8e-17;   
    K(6)=1.5e-16; 
    K(9)=5e-13;   
    R(4)=K(4)*n(2)*n(5)*ng;
    R(5)=K(5)*n(4)*n(7);
    R(6)=K(6)*n(4)*n(6);
    R(9)=K(9)*n(7)*n(8);
end

%Wall losses  
W(1)=0;
W(2)=0;
W(3)=0;
W(4)=0;
for loop=5:8
    W(loop)=-area_volume*n(loop)*sqrt(1.6e-19*Te/mass(loop));
end

W(1)=W(8)-sum(W(5:7));

% e: species 1
isp=1;
X(isp,1)=-R(1);
X(isp,2)=R(2);
X(isp,3)=R(3);
X(isp,7)=-R(7);
X(isp,8)=R(8);

if(~only_e_balance)
    % CO: species 2
    isp=2;
    X(isp,1)=R(1);
    X(isp,2)=-R(2);
    X(isp,3)=-R(3);
    
    % O2: species 3
    isp=3;
    X(isp,6)=R(6);
    X(isp,7)=-R(7);
    
    % O: species 4
    isp=4;
    X(isp,5)=-R(5);
    X(isp,6)=-R(6);
    X(isp,7)=R(7);
    X(isp,9)=R(9);
    
    % O—: species 5
    isp=5;
    X(isp,1)=R(1);
    X(isp,2)=-R(2);
    X(isp,4)=-R(4);
    X(isp,6)=R(6);
    X(isp,7)=R(7);
    
    % O2-: species 6
    isp=6;
    X(isp,6)=R(6);
    X(isp,7)=-R(7);
    
    % CO3-: species 7
    isp=7;
    X(isp,3)=-R(3);
    X(isp,4)=R(4);
    X(isp,5)=-R(5);
    X(isp,9)=-R(9);
    
    % CO2+; species 8
    isp=8;
    X(isp,8)=R(8);
    X(isp,9)=-R(9);
end

%Wall losses

X(:,(length(R)+1))=W;

%Map reactions to wavelengths
%This can be use to generate synthetic spectra
lambda=zeros(1,length(R));
%Example: Reaction 4 radiates at 850nm
lambda(4)= 850; %此处非准确数值
end

