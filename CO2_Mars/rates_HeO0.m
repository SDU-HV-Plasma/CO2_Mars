function [nsp,species,K,number_of_reactions,mass,aux,auxv,Kv,number_of_reaction,number_of_vreaction]=rates_HeO0(Tg,Vibration)

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

global species 

species{1}='e';              mass(1)=9.1e-31;
species{2}='He';             mass(2)=6.697e-27;
species{3}='O';              mass(3)=2.66e-26;
species{4}='O_2';            mass(4)=5.3e-26;
species{5}='O_3';            mass(5)=7.97e-26;
species{6}='O^-';            mass(6)=2.66e-26;
species{7}='O_2^-';          mass(7)=5.3e-26;
species{8}='O_3^-';          mass(8)=7.97e-26;
species{9}='He^+';           mass(9)=6.697e-27;
species{10}='He_2^+';        mass(10)=1.34e-26;
species{11}='O^+';           mass(11)=2.66e-26;
species{12}='O_2^+';         mass(12)=5.3e-26;
species{13}='He^*';          mass(13)=6.697e-27;
species{14}='He_2^*';        mass(14)=1.34e-26;
species{15}='O(1D)';         mass(15)=2.66e-26;
species{16}='O2(Δ)';        mass(16)=5.3e-26;
species{17}='O2(Σ)';        mass(17)=5.3e-26;

nsp=length(species);
number_of_reaction=37;
number_of_vreaction=30;
number_of_reactions = number_of_reaction + number_of_vreaction;
K(1:number_of_reaction)=0;
Kv=zeros(1,number_of_vreaction);
aux=[];
auxv=[];

aux{36}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-e_Heelastic.lut');
aux{37}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-e_O2elastic.lut');

%% Function without Te
% 基态存在电子得失
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
    Kv(4)=3.8e-17;
    Kv(6)=1.0e-14;
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

