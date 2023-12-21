function [X,R,dTe,dS_ape,dS_el,dS_epro,dS_rea,dS_wall,K0,Rv,Kv0,Ralldlt]=rates_HeO(Te,n,area_volume,Tg,me,e,S_abs,dt,Vibration,K,nsp,number_of_reactions,gas_heat,mass,aux,auxv,Kv,number_of_reaction,number_of_vreaction,nv,sim)

% R1:  e + He => He^+ + 2e          k= 2.5e-18*Te^(0.68)*exp(-24.6/Te)    ER=24.6eV  
% R2:  e + He2^+ => 2He             k= 1.0e-14
% R3:  e + He2^+ + He => 3He        k= 2.0e-39
% R4:  e + O2 => O^- + O            k= 8.8e-17*exp(-4.4/Te)               ER=4.4eV
% R5:  e + O2^+ => 2O               k= 4.8e-13
% R6:  e + O^- => O + 2e            k= 2.0e-13*exp(-5.5/Te)               ER=5.5eV
% R7:  e + O2 => O^- + O^+ + 2e     k= 7.1e-17*Te^(0.5)*exp(-17/Te)       ER=17eV
% R8:  e + O2 => O + O^+ + 2e       k= 5.3e-16*Te^(0.9)*exp(-20/Te)       ER=20eV
% R9:  e + O2 => 2O + e             k= 4.2e-15*exp(-5.6/Te)               ER=5.6eV
% R10: e + O3 => O2^- + O           k= 1.0e-15
% R11: e + O3 => O^- + O2           k= 9.3e-16/(Te^(0.62))
% R12: e + O3 => O + O2 + e         k= 1.0e-14*(300/(Tg*11600))^0.5
% R13: e + O => O^+ + 2e            k= 9.0e-15*Te^(0.7)*exp(-13.6/Te)     ER=13.6eV
% R14: e + O + O2 => O^- + O2       k= 1.0e-43
% R15: e + O2 + He => O2^- + He     k= 1.0e-43
% R16: O + O + He => O2 + He        k= 1.3e-45*(3000/(Tg*11600))*exp(-170/(Tg*16000))
% R17: O + He + O2 => He + O3       k= 6.27e-46
% R18: O + O + O2 => O3 + O         k= 3.4e-46*exp(345/(Tg*11600))
% R19: O + O2 + O2 => O3 + O2       k= 6.9e-46*(300/(Tg*11600))^1.25
% R20: O + O + O2 => 2O2            k= 2.56e-46*(300/(Tg*11600))^0.63
% R21: O + O + O => O2 + O          k= 9.21e-46*(300/(Tg*11600))^0.63
% R22: O^- + O2 => O3 + e           k= 5.0e-18*(300/(Tg*11600))^(-0.5)
% R23: O^- + O3 => O3^- + O         k= 5.3e-16
% R24: O^- + O2^+ => 3O             k= 1.0e-13
% R25: O^- + O2^+ => O + O2         k= 2.0e-13*(200/(Tg*11600))^0.5
% R26: O^- + O^+ => 2O              k= 2.0e-13*(300/(Tg*11600))^0.5
% R27: O2^- + O^+ => O2 + O         k= 2.0e-13*(300/(Tg*11600))
% R28: O2^- + O2^+ => 2O2           k= 2.0e-13*(300/(Tg*11600))
% R29: O2^- + O3 => O3^- + O2       k= 6.0e-16*(300/(Tg*11600))^-0.5
% R30: O2 + O^+ => O2^+ + O         k= 2.0e-17*(300/(Tg*11600))^0.5
% R31: O3 + O => 2O2                k= 1.5e-17*exp(-2250/(Tg*11600))
% R32: O3^- + O => O2^- + O2        k= 3.2e-16
% R33: O3^- + O2^+ => 2O + O3       k= 2.0e-13
% R34: He + O3 => He + O2 + O       k= 2.28e-32
% R35: He^+ + 2He => He2^+ + He     k= 1.0e-43
% R36: e + He => e + He             k= 查表
% R37: e + O2 => e + O2             k= 查表
% 振动激发态部分
% Rv1: e + He => He^* + e           k= 2.3e-16*Te^(0.31)*exp(-19.8/Te)        ER=19.8eV
% Rv2: e + He^* => e + He           k= 2.0e-16
% Rv3: e + He^* => He^+ + 2e        k= 1.28e-13*Te^(0.6)*exp(-4.78/Te)        ER=4.78eV
% Rv4: e + He2^* => 2He + e         k= 3.8e-17
% Rv5: e + He2^* => He2^+ + 2e      k= 9.75e-16*Te^(0.71)*exp(-3.4/Te)        ER=3.4eV
% Rv6: He2^* + He => 3He            k= 1.0e-14
% Rv7: e + He2^+ => He + He^*       k= 5.0e-15*Te^(-0.5)
% Rv8: e + O2 => O + O(1D) + e      k= 5.0e-14*exp(-8.4/Te)                   ER=8.4eV
% Rv9: e + O2(Δ) => e + O2         k= 5.6e-15*exp(-2.2/Te)                   ER=2.2eV
% Rv10:e + O2 => O2(Δ) + e         k= 1.7e-15*exp(-3.1/Te)                   ER=3.1eV
% Rv11:e + O(1D) => O^+ + 2e        k= 9.0e-15*Te^(0.5)*exp(-12.6/Te)         ER=12.6eV
% Rv12:e + O(1D) => e + O           k= 8.0e-15
% Rv13:e + O => O(1D) + e           k= 4.2e-15*exp(-2.25/Te)                  ER=2.25eV
% Rv14:O^- + O2(Δ) => O2^- + O     k= 1.0e-16
% Rv15:O^- + O2(Δ) => O3 + e       k= 3.0e-16
% Rv16:O(1D) + O2 => O + O2         k= 4.8e-18*exp(-67/(Tg*11600))
% Rv17:O(1D) + O2 => O + O2(Σ)     k= 2.0e-17
% Rv18:O(1D) + O => 2O              k= 8.0e-18
% Rv19:O(1D) + O3 => O2(Δ) + O2    k= 2.7e-16
% Rv20:O(1D) + O3 => 2O2(Δ)        k= 2.5e-16
% Rv21:O(1D) + O3 => 2O + O2        k= 1.2e-16
% Rv22:O(1D) + O3 => 2O2            k= 1.2e-16
% Rv23:O2(Δ) + He => O2 + He       k= 8.0e-27*(Tg*11600/300)^0.5
% Rv24:O2(Δ) + O3 => 2O2 + O(1D)   k= 1.0e-17
% Rv25:O2(Σ) + O3 => O + 2O2       k= 1.5e-17
% Rv26:O2(Σ) + O => O2(Δ) + O     k= 3.0e-18
% Rv27:O2 + O2(Δ) => 2O2           k= 2.2e-24*(Tg*11600/300)^0.8
% Rv28:He^* + O2 => O2^+ + He + e   k= 2.54e-16*(Tg*11600/300)^0.5
% Rv29:He^* + 2He => He2^* + He     k= 1.6e-44
% Rv30:He^* + He^* => He2^+ + e     k= 1.5e-15

Ralldlt=0;
R(1:number_of_reaction)=0;
Rv=zeros(1,number_of_vreaction);
X=zeros(nsp,number_of_reactions+1);
Xgv=zeros(length(n),number_of_vreaction);
Xvv=zeros(length(nv),number_of_vreaction);

Ngg=n*n';
Ngv=n*nv';
Nvv=nv*nv';

%% Function of Te
%基态存在电子得失
K(1)=2.5e-18*Te^(0.68)*exp(-24.6/Te);
K(4)=8.8e-17*exp(-4.4/Te);
K(6)=2.0e-13*exp(-5.5/Te);
K(7)=7.1e-17*Te^(0.5)*exp(-17/Te);
K(8)=5.3e-16*Te^(0.9)*exp(-20/Te);
K(11)=9.3e-16/(Te^(0.62));
K(13)=9.0e-15*Te^(0.7)*exp(-13.6/Te);

%基态无电子得失
K(9)=4.2e-15*exp(-5.6/Te);
K(36)=interp1(aux{36}(:,1),aux{36}(:,2),Te);
K(37)=interp1(aux{37}(:,1),aux{37}(:,2),Te);

% 振动态
if Vibration==1
    Kv(1)=2.3e-16*Te^(0.31)*exp(-19.8/Te);
    Kv(3)=1.28e-13*Te^(0.6)*exp(-4.78/Te);
    Kv(5)=9.75e-18*Te^(0.71)*exp(-3.4/Te);
    Kv(7)=5.0e-15*Te^(-0.5);
    Kv(8)=5.0e-14*exp(-8.4/Te);
    Kv(9)=5.6e-15*exp(-2.2/Te);
    Kv(10)=1.7e-15*exp(-3.1/Te);
    Kv(11)=9.0e-15*Te^(0.5)*exp(-12.6/Te);
    Kv(13)=4.2e-15*exp(-2.25/Te);
end 

%% Function without Te

if gas_heat==1
    K(2)=1.0e-14;
    K(3)=2.0e-39;
    K(5)=4.8e-13;
    K(10)=1.0e-15;
    K(14)=1.0e-43;
    K(15)=1.0e-43;
    
    %基态无电子得失
    K(12)=1.0e-14*(300/(Tg*11600))^0.5;
    K(16)=1.3e-45*(3000/(Tg*11600))*exp(-170/(Tg*11600));
    K(17)=6.27e-46;
    K(18)=3.4e-46*exp(345/(Tg*11600));
    K(19)=6.9e-46*(300/(Tg*11600))^1.25;
    K(20)=2.56e-46*(300/(Tg*11600))^0.63;
    K(21)=9.21e-46*(300/(Tg*11600))^0.63;
    K(22)=5.0e-18*(300/(Tg*11600))^(-0.5);
    K(23)=5.3e-16;
    K(24)=1.0e-13;
    K(25)=2.0e-13*(200/(Tg*11600))^0.5;
    K(27)=2.0e-13*(300/(Tg*11600));
    K(26)=2.0e-13*(300/(Tg*11600))^0.5;
    K(28)=2.0e-13*(300/(Tg*11600));
    K(29)=6.0e-16*(300/(Tg*11600))^-0.5;
    K(30)=2.0e-17*(300/(Tg*11600))^0.5;
    K(31)=1.5e-17*exp(-2250/(Tg*11600));
    K(32)=3.2e-16;
    K(33)=2.0e-13;
    K(34)=2.28e-32;
    K(35)=1.0e-43;
    
    %振动态
    if Vibration==1
        Kv(2)=2.0e-16;
        Kv(4)=3.8e-15;
        Kv(6)=1.0e-12;
        Kv(12)=8.0e-15;
        Kv(14)=1.0e-16;
        Kv(15)=3.0e-16;
        Kv(16)=4.8e-18*exp(-67/(Tg*11600));
        Kv(17)=2.0e-17;
        Kv(18)=8.0e-18;
        Kv(19)=2.7e-16;
        Kv(20)=2.5e-16;
        Kv(21)=1.2e-16;
        Kv(22)=1.2e-16;
        Kv(23)=8.0e-27*(Tg*11600/300)^0.5;
        Kv(24)=1.0e-17;
        Kv(25)=1.5e-17;
        Kv(26)=3.0e-18;
        Kv(27)=2.2e-24*(Tg*11600/300)^0.8;
        Kv(28)=2.54e-16*(Tg*11600/300)^0.5;
        Kv(29)=1.6e-44;
        Kv(30)=1.5e-15;
    end   
end
    if sim==1
        Ralldlt=[];
        K(intersect([1:length(K)],Ralldlt))=0;
        Kv(intersect([1:length(Kv)],Ralldlt-length(K)))=0;
    end
    K0=K;
    Kv0=Kv;

    %% dn
    % 基态存在电子得失
    R(1)=K(1)*Ngg(1,2);
    R(2)=K(2)*Ngg(1,10);
    R(3)=K(3)*Ngg(1,10)*n(2);
    R(4)=K(4)*Ngg(1,4);
    R(5)=K(5)*Ngg(1,12);
    R(6)=K(6)*Ngg(1,6);
    R(7)=K(7)*Ngg(1,4);
    R(8)=K(8)*Ngg(1,4);
    R(10)=K(10)*Ngg(1,5);
    R(11)=K(11)*Ngg(1,5);
    R(13)=K(13)*Ngg(1,3);
    R(14)=K(14)*Ngg(1,3)*n(4);
    R(15)=K(15)*Ngg(1,4)*n(2);
    R(22)=K(22)*Ngg(4,6);
    
    % 基态无电子得失
    R(9)=K(9)*Ngg(1,4);
    R(12)=K(12)*Ngg(1,5);
    R(16)=K(16)*Ngg(3,3)*n(2);
    R(17)=K(17)*Ngg(3,4)*n(2);
    R(18)=K(18)*Ngg(3,4)*n(3);
    R(19)=K(19)*Ngg(3,4)*n(4);
    R(20)=K(20)*Ngg(3,3)*n(4);
    R(21)=K(21)*Ngg(3,3)*n(3);
    R(23)=K(23)*Ngg(5,6);
    R(24)=K(24)*Ngg(6,12);
    R(25)=K(25)*Ngg(6,12);
    R(26)=K(26)*Ngg(6,11);
    R(27)=K(27)*Ngg(7,11);
    R(28)=K(28)*Ngg(7,12);
    R(29)=K(29)*Ngg(5,7);
    R(30)=K(30)*Ngg(4,11);
    R(31)=K(31)*Ngg(3,5);
    R(32)=K(32)*Ngg(3,8);
    R(33)=K(33)*Ngg(8,12);
    R(34)=K(34)*Ngg(2,5);
    R(35)=K(35)*Ngg(2,9)*n(2);
    R(36)=K(36)*Ngg(1,2);
    R(37)=K(37)*Ngg(1,4);
    
    % 振动态
if Vibration==1
    Rv(1)=Kv(1)*Ngg(1,2);
    Rv(2)=Kv(2)*Ngv(1,1);
    Rv(3)=Kv(3)*Ngv(1,1);
    Rv(4)=Kv(4)*Ngv(1,2);
    Rv(5)=Kv(5)*Ngv(1,2);
    Rv(6)=Kv(6)*Ngv(2,2);
    Rv(7)=Kv(7)*Ngg(1,10);
    Rv(8)=Kv(8)*Ngg(1,4);
    Rv(9)=Kv(9)*Ngv(1,4);
    Rv(10)=Kv(10)*Ngg(1,4);
    Rv(11)=Kv(11)*Ngv(1,3);
    Rv(12)=Kv(12)*Ngv(1,3);
    Rv(13)=Kv(13)*Ngg(1,3);
    Rv(14)=Kv(14)*Ngv(6,4);
    Rv(15)=Kv(15)*Ngv(6,4);
    Rv(16)=Kv(16)*Ngv(4,3);
    Rv(17)=Kv(17)*Ngv(4,3);
    Rv(18)=Kv(18)*Ngv(3,3);
    Rv(19)=Kv(19)*Ngv(5,3); 
    Rv(20)=Kv(20)*Ngv(5,3);
    Rv(21)=Kv(21)*Ngv(5,3);
    Rv(22)=Kv(22)*Ngv(5,3);
    Rv(23)=Kv(23)*Ngv(2,4);
    Rv(24)=Kv(24)*Ngv(5,4);
    Rv(25)=Kv(25)*Ngv(5,5);
    Rv(26)=Kv(26)*Ngv(3,5);
    Rv(27)=Kv(27)*Ngv(4,4);
    Rv(28)=Kv(28)*Ngv(4,1);
    Rv(29)=Kv(29)*Ngv(2,1)*n(2);
    Rv(30)=Kv(30)*Nvv(1,1);
end
%% Wall losses
W(1:nsp)=0;
W(9)=-1/10*area_volume*n(9)*sqrt(1.6e-19*Te/mass(9));
W(10)=-1/10*area_volume*n(10)*sqrt(1.6e-19*Te/mass(10));
W(11)=-1/10*area_volume*n(11)*sqrt(1.6e-19*Te/mass(11));
W(12)=-1/10*area_volume*n(12)*sqrt(1.6e-19*Te/mass(12));
W(2)=-W(9)-2*W(10);
W(4)=-W(12)-0.01*W(11);
W(3)=-0.98*W(11);
W(1)=W(9)+W(10)+W(11)+W(12);
Wv=zeros(1,nsp-length(n));
if Vibration==1
    Wv(1)=-1/2*area_volume*nv(1)*sqrt(8*1.6e-19*Tg/(pi*mass(2)));
    Wv(2)=-1/2*area_volume*nv(2)*sqrt(8*1.6e-19*Tg/(2*pi*mass(2)));
    Wv(3)=-1/2*area_volume*nv(3)*sqrt(8*1.6e-19*Tg/(pi*mass(3)));
    Wv(4)=-1/2*area_volume*nv(4)*sqrt(8*1.6e-19*Tg/(pi*mass(4)))*5e-5;
    Wv(5)=-1/2*area_volume*nv(5)*sqrt(8*1.6e-19*Tg/(pi*mass(4)))*2e-2;
end
% if sim==1
%     Wv([2,6,12,13,18])=0;
%     W(1,length(n)+1:nsp)=Wv;
% end
W(length(n)+1:nsp)=Wv;
W(4)=W(4)-0.0002*Wv(3)-0.5*Wv(4)-0.8*Wv(5);
W(2)=W(2)-Wv(1)-Wv(2);
if W(1)<=-0.1*n(1)/dt
    W(2:length(n))=-W(2:length(n))*0.1*n(1)/dt/W(1);
    W(1)=-0.1*n(1)/dt;
end

% e: species 1
isp=1;
X(isp,1)=R(1);
X(isp,2)=-R(2);
X(isp,3)=-R(3);
X(isp,4)=-R(4);
X(isp,5)=-R(5);
X(isp,6)=R(6);
X(isp,7)=R(7);
X(isp,8)=R(8);
X(isp,10)=-R(10);
X(isp,11)=-R(11);
X(isp,13)=R(13);
X(isp,14)=-R(14);
X(isp,15)=-R(15);
X(isp,22)=R(22);
Xgv(isp,3)=Vibration.*Rv(3);
Xgv(isp,5)=Vibration.*Rv(5);
Xgv(isp,7)=-Vibration.*Rv(7);
Xgv(isp,11)=Vibration.*Rv(11);
Xgv(isp,15)=Vibration.*Rv(15);
Xgv(isp,28)=Vibration.*Rv(28);
Xgv(isp,30)=Vibration.*Rv(30);

% He: species 2
isp=2;
X(isp,1)=-R(1);
X(isp,2)=2*R(2);
X(isp,3)=2*R(3);
X(isp,35)=-R(35);
Xgv(isp,1)=-Vibration.*Rv(1);
Xgv(isp,2)=Vibration.*Rv(2);
Xgv(isp,4)=2*Vibration.*Rv(4);
Xgv(isp,6)=2*Vibration.*Rv(6);
Xgv(isp,7)=Vibration.*Rv(7);
Xgv(isp,28)=Vibration.*Rv(28);
Xgv(isp,29)=-Vibration.*Rv(29);

% O: species 3
isp=3;
X(isp,4)=R(4);
X(isp,5)=2*R(5);
X(isp,6)=R(6);
X(isp,8)=R(8);
X(isp,9)=2*R(9);
X(isp,10)=R(10);
X(isp,12)=R(12);
X(isp,13)=-R(13);
X(isp,14)=-R(14);
X(isp,16)=-2*R(16);
X(isp,17)=-R(17);
X(isp,18)=-R(18);
X(isp,19)=-R(19);
X(isp,20)=-2*R(20);
X(isp,21)=-2*R(21);
X(isp,23)=R(23);
X(isp,24)=3*R(24);
X(isp,25)=R(25);
X(isp,26)=2*R(26);
X(isp,27)=R(27);
X(isp,30)=R(30);
X(isp,31)=-R(31);
X(isp,32)=-R(32);
X(isp,33)=2*R(33);
X(isp,34)=R(34);
Xgv(isp,8)=Vibration.*Rv(8);
Xgv(isp,12)=Vibration.*Rv(12);
Xgv(isp,13)=-Vibration.*Rv(13);
Xgv(isp,14)=Vibration.*Rv(14);
Xgv(isp,16)=Vibration.*Rv(16);
Xgv(isp,17)=Vibration.*Rv(17);
Xgv(isp,18)=Vibration.*Rv(18);
Xgv(isp,21)=2*Vibration.*Rv(21);
Xgv(isp,25)=Vibration.*Rv(25);

% O2: species 4
isp=4;
X(isp,4)=-R(4);
X(isp,7)=-R(7);
X(isp,8)=-R(8);
X(isp,9)=-R(9);
X(isp,11)=R(11);
X(isp,12)=R(12);
X(isp,15)=-R(15);
X(isp,16)=R(16);
X(isp,17)=-R(17);
X(isp,18)=-R(18);
X(isp,19)=-R(19);
X(isp,20)=R(20);
X(isp,21)=R(21);
X(isp,22)=-R(22);
X(isp,25)=R(25);
X(isp,27)=R(27);
X(isp,28)=2*R(28);
X(isp,29)=R(29);
X(isp,30)=-R(30);
X(isp,31)=2*R(31);
X(isp,32)=R(32);
X(isp,34)=R(34);
Xgv(isp,8)=-Vibration.*Rv(8);
Xgv(isp,9)=Vibration.*Rv(9);
Xgv(isp,10)=-Vibration.*Rv(10);
Xgv(isp,17)=-Vibration.*Rv(17);
Xgv(isp,19)=Vibration.*Rv(19);
Xgv(isp,21)=Vibration.*Rv(21);
Xgv(isp,22)=2*Vibration.*Rv(22);
Xgv(isp,23)=Vibration.*Rv(23);
Xgv(isp,24)=2*Vibration.*Rv(24);
Xgv(isp,25)=2*Vibration.*Rv(25);
Xgv(isp,27)=Vibration.*Rv(27);
Xgv(isp,28)=-Vibration.*Rv(28);

% O3: species 5
isp=5;
X(isp,10)=-R(10);
X(isp,11)=-R(11);
X(isp,12)=-R(12);
X(isp,17)=R(17);
X(isp,18)=R(18);
X(isp,19)=R(19);
X(isp,22)=R(22);
X(isp,23)=-R(23);
X(isp,29)=-R(29);
X(isp,31)=-R(31);
X(isp,33)=R(33);
X(isp,34)=-R(34);
Xgv(isp,15)=Vibration.*Rv(15);
Xgv(isp,19)=-Vibration.*Rv(19);
Xgv(isp,20)=-Vibration.*Rv(20);
Xgv(isp,21)=-Vibration.*Rv(21);
Xgv(isp,22)=-Vibration.*Rv(22);
Xgv(isp,24)=-Vibration.*Rv(24);
Xgv(isp,25)=-Vibration.*Rv(25);

% O-: species 6
isp=6;
X(isp,4)=-R(4);
X(isp,6)=-R(6);
X(isp,7)=R(7);
X(isp,11)=R(11);
X(isp,14)=R(14);
X(isp,22)=-R(22);
X(isp,23)=-R(23);
X(isp,24)=-R(24);
X(isp,25)=-R(25);
X(isp,26)=-R(26);
Xgv(isp,14)=-Vibration.*Rv(14);
Xgv(isp,15)=-Vibration.*Rv(15);

% O2-: species 7
isp=7;
X(isp,10)=R(10);
X(isp,15)=R(15);
X(isp,27)=-R(27);
X(isp,28)=-R(28);
X(isp,29)=-R(29);
X(isp,32)=R(32);
Xgv(isp,14)=Vibration.*Rv(14);

% O3-: species 8
isp=8;
X(isp,23)=R(23);
X(isp,29)=R(29);
X(isp,32)=-R(32);
X(isp,33)=-R(33);

% He+: species 9
isp=9;
X(isp,1)=R(1);
X(isp,35)=-R(35);
Xgv(isp,3)=Vibration.*Rv(3);

% He2+: species 10
isp=10;
X(isp,2)=-R(2);
X(isp,3)=-R(3);
X(isp,35)=R(35);
Xgv(isp,5)=Vibration.*Rv(5);
Xgv(isp,7)=-Vibration.*Rv(7);
Xgv(isp,30)=Vibration.*Rv(30);

% O+: species 11
isp=11;
X(isp,7)=R(7);
X(isp,8)=R(8);
X(isp,13)=R(13);
X(isp,26)=-R(26);
X(isp,27)=-R(27);
X(isp,30)=-R(30);
Xgv(isp,11)=Vibration.*Rv(11);

% O2+: species 12
isp=12;
X(isp,5)=-R(5);
X(isp,24)=-R(24);
X(isp,25)=-R(25);
X(isp,28)=-R(28);
X(isp,30)=R(30);
X(isp,33)=-R(33);
Xgv(isp,28)=Vibration.*Rv(28);
   
if Vibration==1
    % He*: species 13, vspecies 1
    isp=1;
    Xvv(isp,1)=Rv(1);
    Xvv(isp,2)=-Rv(2);
    Xvv(isp,3)=-Rv(3);
    Xvv(isp,7)=Rv(7);
    Xvv(isp,28)=-Rv(28);
    Xvv(isp,29)=-Rv(29);
    Xvv(isp,30)=-2*Rv(30);
    
    % He2*: species 14, vspecies 2
    isp=2;
    Xvv(isp,4)=-Rv(4);
    Xvv(isp,5)=-Rv(5);
    Xvv(isp,6)=-Rv(6);
    Xvv(isp,29)=Rv(29);
    
    % O(1D): species 15, vspecies 3
    isp=3;
    Xvv(isp,8)=Rv(8);
    Xvv(isp,11)=-Rv(11);
    Xvv(isp,12)=-Rv(12);
    Xvv(isp,13)=Rv(13);
    Xvv(isp,16)=-Rv(16);
    Xvv(isp,17)=-Rv(17);
    Xvv(isp,18)=-Rv(18);
    Xvv(isp,19)=-Rv(19);
    Xvv(isp,20)=-Rv(20);
    Xvv(isp,21)=-Rv(21);
    Xvv(isp,22)=-Rv(22);
    Xvv(isp,24)=Rv(24);
   
    % O2(delta): species 16, vspecies 4
    isp=4;
    Xvv(isp,9)=-Rv(9);
    Xvv(isp,10)=Rv(10);
    Xvv(isp,14)=-Rv(14);
    Xvv(isp,15)=-Rv(15);
    Xvv(isp,19)=Rv(19);
    Xvv(isp,20)=2*Rv(20);
    Xvv(isp,23)=-Rv(23);
    Xvv(isp,24)=-Rv(24);
    Xvv(isp,26)=Rv(26);
    Xvv(isp,27)=-Rv(27);
    
    % O2(sigma): species 17, vspecies 5
    isp=5;
    Xvv(isp,17)=Rv(17);
    Xvv(isp,25)=-Rv(25);
    Xvv(isp,26)=-Rv(26);
end

X(1:length(n),number_of_reaction+1:number_of_reactions)=Xgv;
X(length(n)+1:nsp,number_of_reaction+1:number_of_reactions)=Xvv;

X(:,number_of_reactions+1)=W';

if n(1)<=1e6 & sum(X(1,:))<0
    inde=X(1,1:number_of_reaction)<0;
    X(:,inde)=0;
end
p(2:length(n)) = n(2:end)<=1e-3 & sum(X(2:length(n),:),2)<0;
p(length(n)+1:nsp) = nv(1:end)<=1e-3 & sum(X(length(n)+1:nsp,:),2)<0;
[~,ind]=find(p'.*X<0);
X(:,ind)=0;

% 能量变化

dS_ape=2/(3*e*n(1))*S_abs*dt;
dS_el=-2*me/mass(2)*K(36)*n(2)*(Te-Tg)*dt-2*me/mass(4)*K(37)*n(4)*(Te-Tg)*dt;
dS_epro=-Te/n(1).*(sum(X(1,:)))*dt;
dS_rea=-2/3*(24.6*K(1)*n(2)+4.4*K(4)*n(4)+5.5*K(6)*n(6)+17*K(7)*n(4)+20*K(8)*n(4)+5.6*K(9)*n(4)+13.6*K(13)*n(3)+...+
    Vibration.*(19.8*Kv(1)*n(2)+4.78*Kv(3)*nv(1)+3.4*Kv(5)*nv(2)+8.4*Kv(8)*n(4)+2.2*Kv(9)*nv(4)+3.1*Kv(10)*n(4)+12.6*Kv(11)*nv(3)+2.25*Kv(13)*n(3)))*dt;
dS_wall=-2/(3*n(1))*(-W(1))*(Te/2*log(mass(2)/(2*pi*me))+5/2*Te)*dt;%(Te/2*log(mg/(2*pi*me))+5/2*Te)为1个电子-离子对的能量

dTe=dS_ape+dS_el+dS_epro+dS_rea+dS_wall;

end

