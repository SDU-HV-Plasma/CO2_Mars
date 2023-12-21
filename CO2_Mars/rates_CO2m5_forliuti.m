function [nsp,species,X,R,lambda,dTe,B,Y,dS_ape,dS_el,dS_epro,dS_rea,dS_wall,number_of_reactions]=rates_CO2m5_forliuti(Te,P,n,area_volume,Tg,only_e_balance,me,mg,e,S_abs,dt)

% R1:  e + CO2 => CO2^+ + 2e          k= 查表                                                ER=13.7eV
% R2:  e + CO2 => CO + O + e          k= 查表                   (激发解离)                   ER=11.46eV
% R3:  e + CO2 => CO + O^-            k= 查表                                            （附着反应无ER）
% R4:  e + O3  => O + O2 + e          k= 9.0e-16                (激发解离)           (应该存在ER，但未找到数据)
% R5:  e + O2  => 2O + e              k= 查表            (生成的O一般而言属于激发态)          ER=6.12eV
% R6:  e + O2  => O + O^-             k= 查表                                             (附着反应无ER)
% R7:  O^- + CO => CO2 + e            k= 5.5e-16
% R8:  O^- + O3 => 2O2 + e            k= 3.0e-16
% R9:  e + CO2^+ => CO + O            k= 2.0e-11*Te^(-0.5)*(Tg*1.16e4)^(-1)
% R10: O2^- + CO2^+ => CO + O2 + O    k= 6.0e-13
% R11: O + O2 + CO2 => O3 + CO2       k= 1.7e-42*(Tg*1.16e4)^(-1.2)                       (O+O2+M,M为CO2)
% R12: O^- + O2 => O + e + O2         k= 6.9e-16                                           (O- + M ,M为O2)
% R13: e + CO2 => e + CO2             k= 查表                                            （弹性碰撞反应ER=0）
% R14: O + O2^- => O2 + O^-           k= 3.31e-16
% R15: O^- + CO2 => O + CO2 + e       k= 4.0e-18                                          （O- + M ,M为CO2）
% R16: O^- + O => O2 + e              k= 2.3e-16
% R17: e + CO => C + O^-              k= 查表                                                 ER=9.2eV
% R18: e + CO => e + C + O            k= 查表                                                 ER=11.1eV
% R19: e + CO2^+ => C + O2            k= 3.94e-13*Te^(-0.4)
% R20: CO2 + C => 2CO                 k= 1.0e-21
% R21: C + O2 => O + CO               k= 3.0e-17
% R22: e + O2 => 2e + O2^+            k= 查表                                               ER=12.06eV
% R23: e + O2^+ => 2O                 k= 6.0e-13*(1/Te)^0.5*(1/(Tg*1.16e4))^0.5
% R24: O^- + O2^+ => O2 + O           k= 2.6e-14*(300/(Tg*1.16e4))^0.44
% R25: O^- + O2^+ => 3O               k= 4.2e-13*(300/(Tg*1.16e4))^0.44
% R26: O2^+ + O2^- => 2O2             k= 2.01e-13*(300/(Tg*1.16e4))^0.5
% R27: O2^+ + O2^- => O2 + 2O         k= 4.2e-13
% R28: O^- + CO2 + CO2 => CO3^- + CO2 k= 9.0e-41                                       （O- + CO2 + M,M为CO2）
% R29: CO3^- + CO => 2CO2 + e         k= 5.0e-19
% R30: CO3^- + CO2^+ => 2CO2 + O      k= 5.0e-13                                       （生成的CO2为vb振动态)
% R31: O + CO3^- => CO2 + O2^-        k= 8.0e-17
% R32: CO2 + e => O2^+ + C + 2e       k= K(1)*1/13                                            ER未知
% R33: O3 + e => O2^+ + O + 2e        k= 3.2e-17*(Te*1.16e4)^(0.5)*(1+0.15*Te)*exp(-12.93/Te)
% R34: O3 + e => O2 + O^-             k= 查表                                             （附着反应无ER）
% R35: CO2^+ + O => O2^+ + CO         k= 0.63*2.6e-16
% R36: CO2^+ + O2 => O2^+ + CO2       k= 6.4e-17
% R37: O3 + e => O2^- + O             k= 查表                                             （附着反应无ER）
% R38: CO3^- + O2^+ => CO2 + O2 + O   k= 3.0e-13
% R39: O2^- + O => O3 + e             k= 3.3e-16
% R40: CO2 + O^- + CO => CO3^- + CO   k= 1.5e-40                                       (O- + CO2 + M,M为CO)
% R41: CO2 + O^- + O2 => CO3^- + O2   k= 3.1e-40                                       （O- + CO2 + M,M为O2）
% R42: O + O + CO2 => O2 + CO2        k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))（O+O+M，其中M为CO2）
% R43: O^- + CO2^+ => O + CO2         k= 1.0e-13
% R44: e + CO2 => 2e + O + CO^+       k= 查表                                                ER=19.5eV
% R45: e + CO => 2e + CO^+            k= 查表                                                ER=14.01eV
% R46: CO2 + CO^+ => CO2^+ + CO       k= 1.0e-15
% R47: O2^+ + C => CO^+ + O           k= 5.2e-17
% R48: e + O3 + O2 => O3^- + O2       k= 4.6e-40                                          (e+O3+M,M为O2)
% R49: O3^- + CO2 => CO3^- + O2       k= 5.5e-16
% R50: O^- + O3 => O3^- + O           k= 8.0e-16
% R51: O2^- + O3 => O3^- + O2         k= 4.0e-16
% R52: O3^- + O => O2^- + O2          k= 2.5e-16
% R53: e + CO2 => 2e + CO + O^+       k= 查表                                                ER=19.1eV
% R54: e + O => 2e + O^+              k= 查表                                                ER=13.6eV
% R55: O^+ + CO2 => O2^+ + CO         k= 9.4e-16
% R56: O^+ + CO2 => O + CO2^+         k= 4.5e-16
% R57: CO2^+ + O => CO2 + O^+         k= 9.62e-17
% R58: e + CO => 2e + C^+ + O         k= 查表                                                   ER=22eV
% R59: e + C => 2e + C^+              k= 查表                                                  ER=11.2eV
% R60: C^+ + CO2 => CO^+ + CO         k= 1.1e-15
% R61: O2^+ + C => C^+ + O2           k= 5.2e-17
% R62: O2^- + CO2 + O2 => CO4^- + O2  k= 4.7e-41                                         (O2- + CO2 + M,M为O2)
% R63: CO4^- + CO2^+ => 2CO2 + O2     k= 5.0e-13
% R64: O2^+ + CO4^- => CO2 + 2O2      k= 3.0e-13
% R65: CO4^- + O => CO3^- + O2        k= 1.1e-16
% R66: CO4^- + O => CO2 + O2 + O^-    k= 1.4e-17
% R67: CO4^- + O => CO2 + O3^-        k= 1.4e-16
% R68: e + CO2 => 2e + O2 + C^+       k= 查表                                                  ER=27.8eV
% R69: e + CO => 2e + C + O^+         k= 查表                                                   ER=25eV
% R70: e + 2O2 => O2^- + O2           k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4)) （三体反应无法寻找ER）(e+O2+M,其中M为O2)
% R71: e + O2 => 2e + O + O^+         k= 查表                                                  ER=19.5eV
% R72: e + CO4^+ => CO2 + O2          k= 1.61e-13*Te^(-0.5)
% R73: O2^+ + CO2 + O2 => CO4^+ + O2  k= 2.3e-41                                   (O2+ + CO2 + M,M取O2，但k系数不准确)
% R74: e + CO^+ => C + O              k= 3.46e-14*Te^(-0.48)
% R75: e + O + O2 => O^- + O2         k= 1.0e-43*exp(300/(Tg*1.16e4))                     (e+O+M，其中M为O2)
% R76: O2^- + O2 => 2O2 + e           k= 2.7e-16*((Tg*1.16e4)/300)^(0.5)*exp(-5590/(Tg*1.16e4)) (O2^-+M,其中M为O2)
% R77: CO2 + CO2 => CO + O + CO2      k= 3.91e-16*exp(-49430/(Tg*1.16e4))                   (CO2+M,M为CO2） 
% R78: CO2 + O => CO + O2             k= 2.8e-17*exp(-26500/(Tg*1.16e4))
% R79: CO + O + CO => CO2 + CO        k= 6.54e-45                                           (O+CO+M,M为CO)
% R80: O2 + CO => CO2 + O             k= 4.2e-18*exp(-24000/(Tg*1.16e4))
% R81: O + C + CO2 => CO + CO2        k= 2.14e-41*(Tg/300*1.16e4)^(-3.08)*exp(-2114/(Tg*1.16e4))(O + C + M,M取CO2，但系数k不准确)
% R82: O2 + CO2 => 2O + CO2           k= 3.0e-12*(1/((Tg*1.16e4)))*exp(-59380/(Tg*1.16e4))(O2+M,M取CO2，但系数k不准确)

global B species 

species{1}='e';         charge(1)=-1;       mass(1)=9.1e-31;%kg 
species{2}='CO';        charge(2)=0;        mass(2)=4.65e-26;
species{3}='O_2';       charge(3)=0;        mass(3)=5.3e-26;
species{4}='O';         charge(4)=0;        mass(4)=2.66e-26;
species{5}='O_3';       charge(5)=0;        mass(5)=7.97e-26;
species{6}='O^-';       charge(6)=-1;       mass(6)=2.66e-26;
species{7}='O_2^-';     charge(7)=-1;       mass(7)=5.3e-26;
species{8}='CO_2^+';    charge(8)=1;        mass(8)=7.31e-26;   
species{9}='CO_4^-';    charge(9)=-1;       mass(9)=1.261e-25;
species{10}='CO_2';     charge(10)=0;       mass(10)=7.31e-26;
species{11}='C';        charge(11)=0;       mass(11)=1.99e-26;
species{12}='O_2^+';    charge(12)=1;       mass(12)=5.3e-26;
species{13}='CO_3^-';   charge(13)=-1;      mass(13)=9.97e-26;
species{14}='CO^+';     charge(14)=1;       mass(14)=4.65e-26;
species{15}='O_3^-';    charge(15)=-1;      mass(15)=7.97e-26;
species{16}='O^+';      charge(16)=1;       mass(16)=2.66e-26;
species{17}='C^+';      charge(17)=1;       mass(17)=1.99e-26;
species{18}='CO_4^+';   charge(18)=1;       mass(18)=1.261e-25;    

nsp=length(species);

if(Te==0)
    return;
end

ng=3.54e13*P*1e9;%m-3
 
%Gas phase reactions
number_of_reactions=82;
K(1:number_of_reactions)=0;
R(1:number_of_reactions)=0;
X(nsp,number_of_reactions+1)=0;

%First: Reactions involving generation/loss of electrons      
K(1)=lookups('reaction-rate-cross_section-e_CO2ionization.lut',Te);
K(3)=lookups('reaction-rate-cross_section-e_CO2attachment.lut',Te);
K(6)=lookups('reaction-rate-cross_section-e_O2attachment.lut',Te);
K(7)=5.5e-16;
K(8)=3.0e-16;
K(9)=2.0e-11*Te^(-0.5)*(Tg*1.16e4)^(-1);
K(12)=6.9e-16;
K(15)=4.0e-18;
K(16)=2.3e-16;
K(17)=lookups('reaction-rate-cross_section-e_COattachment.lut',Te);
K(19)=3.94e-13*Te^(-0.4);
K(22)=lookups('reaction-rate-cross_section-e_O2ionization.lut',Te);
K(23)=6.0e-13*(1/Te)^0.5*(1/(Tg*1.16e4))^0.5;
K(29)=5.0e-19;
K(32)=K(1)*1/13;
K(33)=3.2e-17*(Te*1.16e4)^(0.5)*(1+0.15*Te)*exp(-12.93/Te);
K(34)=lookups('reaction-rate-cross_section-e_O3attachment_O-.lut',Te);
K(37)=lookups('reaction-rate-cross_section-e_O3attachment_O2-.lut',Te);
K(39)=3.3e-16;
K(44)=lookups('reaction-rate-cross_section-e_CO2ionization_CO+.lut',Te);
K(45)=lookups('reaction-rate-cross_section-e_COionization.lut',Te);
K(48)=4.6e-40;
K(53)=lookups('reaction-rate-cross_section-e_CO2ionization_O+.lut',Te);
K(54)=lookups('reaction-rate-cross_section-e_Oionization.lut',Te);
K(58)=lookups('reaction-rate-cross_section-e_COionization_C+.lut',Te);
K(59)=lookups('reaction-rate-cross_section-e_Cionization.lut',Te);
K(68)=lookups('reaction-rate-cross_section-e_CO2ionization_C+.lut',Te);
K(69)=lookups('reaction-rate-cross_section-e_COionization_O+.lut',Te);
K(70)=2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4));
K(71)=lookups('reaction-rate-cross_section-e_O2ionization_O+.lut',Te);
K(72)=1.61e-13*Te^(-0.5);
K(74)=3.46e-14*Te^(-0.48);
K(75)=1.0e-43*exp(300/(Tg*1.16e4));
K(76)=2.7e-16*((Tg*1.16e4)/300)^(0.5)*exp(-5590/(Tg*1.16e4));
R(1)=K(1)*n(1)*n(10);
R(3)=K(3)*n(1)*n(10);
R(6)=K(6)*n(1)*n(3);
R(7)=K(7)*n(6)*n(2);
R(8)=K(8)*n(6)*n(5);
R(9)=K(9)*n(1)*n(8);
R(12)=K(12)*n(3)*n(6);
R(15)=K(15)*n(6)*n(10);
R(16)=K(16)*n(6)*n(4);
R(17)=K(17)*n(1)*n(2);
R(19)=K(19)*n(1)*n(8);
R(22)=K(22)*n(1)*n(3);
R(23)=K(23)*n(1)*n(12);
R(29)=K(29)*n(2)*n(13);
R(32)=K(32)*n(1)*n(10);
R(33)=K(33)*n(1)*n(5);
R(34)=K(34)*n(1)*n(5);
R(37)=K(37)*n(1)*n(5);
R(39)=K(39)*n(4)*n(7);
R(44)=K(44)*n(1)*n(10);
R(45)=K(45)*n(1)*n(2);
R(48)=K(48)*n(1)*n(5)*n(3);
R(53)=K(53)*n(1)*n(10);
R(54)=K(54)*n(1)*n(4);
R(58)=K(58)*n(1)*n(2);
R(59)=K(59)*n(1)*n(11);
R(68)=K(68)*n(1)*n(10);
R(69)=K(69)*n(1)*n(2);
R(70)=K(70)*n(1)*n(3)*n(3);
R(71)=K(71)*n(1)*n(3);
R(72)=K(72)*n(1)*n(18);
R(74)=K(74)*n(1)*n(14);
R(75)=K(75)*n(1)*n(3)*n(4);
R(76)=K(76)*n(3)*n(7);

%Second: Rest of reactions
if(~only_e_balance)
    K(2)=lookups('reaction-rate-cross_section-e_CO2dissociation.lut',Te);
    K(4)=9.0e-16;
    K(5)=lookups('reaction-rate-cross_section-e_O2dissociation.lut',Te);
    K(10)=6.0e-13;
    K(11)=1.7e-42*(Tg*1.16e4)^(-1.2);
    K(13)=lookups('reaction-rate-cross_section-e_CO2elastic.lut',Te);%e与CO2发生弹性碰撞的速率系数
    K(14)=3.31e-16;
    K(18)=lookups('reaction-rate-cross_section-e_COexcitation.lut',Te);
    K(20)=1.0e-21;
    K(21)=3.0e-17;
    K(24)=2.6e-14*(300/(Tg*1.16e4))^0.44;
    K(25)=4.2e-13*(300/(Tg*1.16e4))^0.44;
    K(26)=2.01e-13*(300/(Tg*1.16e4))^0.5;
    K(27)=4.2e-13;
    K(28)=9.0e-41;
    K(30)=5.0e-13;
    K(31)=8.0e-17;
    K(35)=0.63*2.6e-16;
    K(36)=6.4e-17;
    K(38)=3.0e-13;
    K(40)=1.5e-40;
    K(41)=3.1e-40;
    K(42)=3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4));
    K(43)=1.0e-13;
    K(46)=1.0e-15;
    K(47)=5.2e-17;
    K(49)=5.5e-16;
    K(50)=8.0e-16;
    K(51)=4.0e-16;
    K(52)=2.5e-16;
    K(55)=9.4e-16;
    K(56)=4.5e-16;
    K(57)=9.62e-17;
    K(60)=1.1e-15;
    K(61)=5.2e-17;
    K(62)=4.7e-41;
    K(63)=5.0e-13;
    K(64)=3.0e-13;
    K(65)=1.1e-16;
    K(66)=1.4e-17;
    K(67)=1.4e-16;
    K(73)=2.3e-41;
    K(77)=3.91e-16*exp(-49430/(Tg*1.16e4));
    K(78)=2.8e-17*exp(-26500/(Tg*1.16e4));
    K(79)=6.54e-45;
    K(80)=4.2e-18*exp(-24000/(Tg*1.16e4));
    K(81)=2.14e-41*(Tg/300*1.16e4)^(-3.08)*exp(-2114/(Tg*1.16e4));
    K(82)=3.0e-12*(1/((Tg*1.16e4)))*exp(-59380/(Tg*1.16e4));
    R(2)=K(2)*n(1)*n(10);
    R(4)=K(4)*n(1)*n(5);
    R(5)=K(5)*n(1)*n(3);
    R(10)=K(10)*n(7)*n(8);
    R(11)=K(11)*n(4)*n(3)*n(10);
    R(13)=K(13)*n(1)*n(10);
    R(14)=K(14)*n(4)*n(7);
    R(18)=K(18)*n(1)*n(2);
    R(20)=K(20)*n(10)*n(11);
    R(21)=K(21)*n(11)*n(3);
    R(24)=K(24)*n(6)*n(12);
    R(25)=K(25)*n(6)*n(12);
    R(26)=K(26)*n(7)*n(12);
    R(27)=K(27)*n(7)*n(12);
    R(28)=K(28)*n(6)*n(10)*n(10);
    R(30)=K(30)*n(8)*n(13);
    R(31)=K(31)*n(4)*n(13);
    R(35)=K(35)*n(4)*n(8);
    R(36)=K(36)*n(3)*n(8);
    R(38)=K(38)*n(12)*n(13);
    R(40)=K(40)*n(2)*n(6)*n(10);
    R(41)=K(41)*n(3)*n(6)*n(10);
    R(42)=K(42)*n(4)*n(4)*n(10);
    R(43)=K(43)*n(6)*n(8);
    R(46)=K(46)*n(10)*n(14);
    R(47)=K(47)*n(11)*n(12);
    R(49)=K(49)*n(10)*n(15);
    R(50)=K(50)*n(5)*n(6);
    R(51)=K(51)*n(5)*n(7);
    R(52)=K(52)*n(4)*n(15);
    R(55)=K(55)*n(10)*n(16);
    R(56)=K(56)*n(10)*n(16);
    R(57)=K(57)*n(4)*n(8);
    R(60)=K(60)*n(10)*n(17);
    R(61)=K(61)*n(11)*n(12);
    R(62)=K(62)*n(7)*n(3)*n(10);
    R(63)=K(63)*n(8)*n(9);
    R(64)=K(64)*n(12)*n(9);
    R(65)=K(65)*n(4)*n(9);
    R(66)=K(66)*n(4)*n(9);
    R(67)=K(67)*n(4)*n(9);
    R(73)=K(73)*n(3)*n(10)*n(12);
    R(77)=K(77)*n(10)*n(10);
    R(78)=K(78)*n(4)*n(10);
    R(79)=K(79)*n(2)*n(2)*n(4);
    R(80)=K(80)*n(2)*n(3);
    R(81)=K(81)*n(4)*n(10)*n(11);
    R(82)=K(82)*n(3)*n(10);
end

%Wall losses
W(2)=0;
W(3)=0;
W(4)=0;
W(5)=0;
W(9)=0;
W(10)=0;
W(11)=0;
W(6)=0;
W(7)=0; 
W(13)=0;
W(15)=0;
W(8)=-1/10*area_volume*n(8)*sqrt(1.6e-19*Te/mass(8));
W(12)=-1/10*area_volume*n(12)*sqrt(1.6e-19*Te/mass(12));
W(14)=-1/10*area_volume*n(14)*sqrt(1.6e-19*Te/mass(14));
W(16)=-1/10*area_volume*n(16)*sqrt(1.6e-19*Te/mass(16));
W(17)=-1/10*area_volume*n(17)*sqrt(1.6e-19*Te/mass(17));
W(18)=-1/10*area_volume*n(18)*sqrt(1.6e-19*Te/mass(18));
W(1)=W(8)+W(12)+W(14)+W(16)+W(17)+W(18);

% e: species 1
isp=1;
X(isp,1)=R(1);
X(isp,3)=-R(3);
X(isp,6)=-R(6);
X(isp,7)=R(7);
X(isp,8)=R(8);
X(isp,9)=-R(9);
X(isp,12)=R(12);
X(isp,15)=R(15);
X(isp,16)=R(16);
X(isp,17)=-R(17);
X(isp,19)=-R(19);
X(isp,22)=R(22);
X(isp,23)=-R(23);
X(isp,29)=R(29);
X(isp,32)=R(32);
X(isp,33)=R(33);
X(isp,34)=-R(34);
X(isp,37)=-R(37);
X(isp,39)=R(39);
X(isp,44)=R(44);
X(isp,45)=R(45);
X(isp,48)=-R(48);
X(isp,53)=R(53);
X(isp,54)=R(54);
X(isp,58)=R(58);
X(isp,59)=R(59);
X(isp,68)=R(68);
X(isp,69)=R(69);
X(isp,70)=-R(70);
X(isp,71)=R(71);
X(isp,72)=-R(72);
X(isp,74)=-R(74);
X(isp,75)=-R(75);
X(isp,76)=R(76);

if(~only_e_balance)
    % CO: species 2
    isp=2;
    X(isp,2)=R(2);
    X(isp,3)=R(3);
    X(isp,7)=-R(7);
    X(isp,9)=R(9);
    X(isp,10)=R(10);
    X(isp,17)=-R(17);
    X(isp,18)=-R(18);
    X(isp,20)=2*R(20);
    X(isp,21)=R(21);
    X(isp,29)=-R(29);
    X(isp,35)=R(35);
    X(isp,45)=-R(45);
    X(isp,46)=R(46);
    X(isp,53)=R(53);
    X(isp,55)=R(55);
    X(isp,58)=-R(58);
    X(isp,60)=R(60);
    X(isp,69)=-R(69);
    X(isp,77)=R(77);
    X(isp,78)=R(78);
    X(isp,79)=-R(79);
    X(isp,80)=-R(80);
    X(isp,81)=R(81);
    
    % O2: species 3
    isp=3;
    X(isp,4)=R(4);
    X(isp,5)=-R(5);
    X(isp,6)=-R(6);
    X(isp,8)=R(8)+R(8);
    X(isp,10)=R(10);
    X(isp,11)=-R(11);
    X(isp,14)=R(14);
    X(isp,16)=R(16);
    X(isp,19)=R(19);
    X(isp,21)=-R(21);
    X(isp,22)=-R(22);
    X(isp,24)=R(24);
    X(isp,26)=2*R(26);
    X(isp,27)=R(27);
    X(isp,34)=R(34);
    X(isp,36)=-R(36);
    X(isp,38)=R(38);
    X(isp,42)=R(42);
    X(isp,49)=R(49);
    X(isp,51)=R(51);
    X(isp,52)=R(52);
    X(isp,61)=R(61);
    X(isp,63)=R(63);
    X(isp,64)=2*R(64);
    X(isp,65)=R(65);
    X(isp,66)=R(66);
    X(isp,68)=R(68);
    X(isp,70)=-R(70);
    X(isp,71)=-R(71);
    X(isp,72)=R(72);
    X(isp,76)=R(76);
    X(isp,78)=R(78);
    X(isp,80)=-R(80);
    X(isp,82)=-R(82);
    
    % O: species 4
    isp=4;
    X(isp,2)=R(2);
    X(isp,4)=R(4);
    X(isp,5)=R(5)+R(5);
    X(isp,6)=R(6);
    X(isp,9)=R(9);
    X(isp,10)=R(10);
    X(isp,11)=-R(11);
    X(isp,12)=R(12);
    X(isp,14)=-R(14);
    X(isp,15)=R(15);
    X(isp,16)=-R(16);
    X(isp,18)=R(18);
    X(isp,21)=R(21);
    X(isp,23)=2*R(23);
    X(isp,24)=R(24);
    X(isp,25)=3*R(25);
    X(isp,27)=2*R(27);
    X(isp,30)=R(30);
    X(isp,31)=-R(31);
    X(isp,33)=R(33);
    X(isp,35)=-R(35);
    X(isp,37)=R(37);
    X(isp,38)=R(38);
    X(isp,39)=-R(39);
    X(isp,42)=-2*R(42);
    X(isp,43)=R(43);
    X(isp,44)=R(44);
    X(isp,47)=R(47);
    X(isp,50)=R(50);
    X(isp,52)=-R(52);
    X(isp,54)=-R(54);
    X(isp,56)=R(56);
    X(isp,57)=-R(57);
    X(isp,58)=R(58);
    X(isp,65)=-R(65);
    X(isp,66)=-R(66);
    X(isp,67)=-R(67);
    X(isp,71)=R(71);
    X(isp,74)=R(74);
    X(isp,75)=-R(75);
    X(isp,77)=R(77);
    X(isp,78)=-R(78);
    X(isp,79)=-R(79);
    X(isp,80)=R(80);
    X(isp,81)=-R(81);
    X(isp,82)=2*R(82);
    
    % O3: species 5
    isp=5;
    X(isp,4)=-R(4);
    X(isp,8)=-R(8);
    X(isp,11)=R(11);
    X(isp,33)=-R(33);
    X(isp,34)=-R(34);
    X(isp,37)=-R(37);
    X(isp,39)=R(39);
    X(isp,48)=-R(48);
    X(isp,50)=-R(50);
    X(isp,51)=-R(51);
    
    % O—: species 6
    isp=6;
    X(isp,3)=R(3);
    X(isp,6)=R(6);
    X(isp,7)=-R(7);
    X(isp,8)=-R(8);
    X(isp,12)=-R(12);
    X(isp,14)=R(14);
    X(isp,15)=-R(15);
    X(isp,16)=-R(16);
    X(isp,17)=R(17);
    X(isp,24)=-R(24);
    X(isp,25)=-R(25);
    X(isp,28)=-R(28);
    X(isp,34)=R(34);
    X(isp,40)=-R(40);
    X(isp,41)=-R(41);
    X(isp,43)=-R(43);
    X(isp,50)=-R(50);
    X(isp,66)=R(66);
    X(isp,75)=R(75);
     
    % O2-: species 7
    isp=7;
    X(isp,10)=-R(10);
    X(isp,14)=-R(14);
    X(isp,26)=-R(26);
    X(isp,27)=-R(27);
    X(isp,31)=R(31);
    X(isp,37)=R(37);
    X(isp,39)=-R(39);
    X(isp,51)=-R(51);
    X(isp,52)=R(52);
    X(isp,62)=-R(62);
    X(isp,70)=R(70);
    X(isp,76)=-R(76);
    
    % CO2+; species 8
    isp=8;
    X(isp,1)=R(1);
    X(isp,9)=-R(9);
    X(isp,10)=-R(10);
    X(isp,19)=-R(19);
    X(isp,30)=-R(30);
    X(isp,35)=-R(35);
    X(isp,36)=-R(36);
    X(isp,43)=-R(43);
    X(isp,46)=R(46);
    X(isp,56)=R(56);
    X(isp,57)=-R(57);
    X(isp,63)=-R(63);
    
    % CO4-; species 9
    isp=9;
    X(isp,62)=R(62);
    X(isp,63)=-R(63);
    X(isp,64)=-R(64);
    X(isp,65)=-R(65);
    X(isp,66)=-R(66);
    X(isp,67)=-R(67);
    
    % CO2; species 10
    isp=10;
    X(isp,1)=-R(1);
    X(isp,2)=-R(2);
    X(isp,3)=-R(3);
    X(isp,7)=R(7);
    X(isp,20)=-R(20);
    X(isp,28)=-R(28);
    X(isp,29)=2*R(29);
    X(isp,30)=2*R(30);
    X(isp,31)=R(31);
    X(isp,32)=-R(32);
    X(isp,36)=R(36);
    X(isp,38)=R(38);
    X(isp,40)=-R(40);
    X(isp,41)=-R(41);
    X(isp,43)=R(43);
    X(isp,44)=-R(44);
    X(isp,46)=-R(46);
    X(isp,49)=-R(49);
    X(isp,53)=-R(53);
    X(isp,55)=-R(55);
    X(isp,56)=-R(56);
    X(isp,57)=R(57);
    X(isp,60)=-R(60);
    X(isp,62)=-R(62);
    X(isp,63)=2*R(63);
    X(isp,64)=R(64);
    X(isp,66)=R(66);
    X(isp,67)=R(67);
    X(isp,68)=-R(68);
    X(isp,72)=R(72);
    X(isp,73)=-R(73);
    X(isp,77)=-R(77);
    X(isp,78)=-R(78);
    X(isp,79)=R(79);
    X(isp,80)=R(80);
    
    % C: species 11
    isp=11;
    X(isp,17)=R(17);
    X(isp,18)=R(18);
    X(isp,19)=R(19);
    X(isp,20)=-R(20);
    X(isp,21)=-R(21);
    X(isp,32)=R(32);
    X(isp,47)=-R(47);
    X(isp,59)=-R(59);
    X(isp,61)=-R(61);
    X(isp,69)=R(69);
    X(isp,74)=R(74);
    X(isp,81)=-R(81);
    
    % O2+: species 12
    isp=12;
    X(isp,22)=R(22);
    X(isp,23)=-R(23);
    X(isp,24)=-R(24);
    X(isp,25)=-R(25);
    X(isp,26)=-R(26);
    X(isp,27)=-R(27);
    X(isp,32)=R(32);
    X(isp,33)=R(33);
    X(isp,35)=R(35);
    X(isp,36)=R(36);
    X(isp,38)=-R(38);
    X(isp,47)=-R(47);
    X(isp,55)=R(55);
    X(isp,61)=-R(61);
    X(isp,64)=-R(64);
    X(isp,73)=-R(73);
    
    % CO3-: species 13
    isp=13;
    X(isp,28)=R(28);
    X(isp,29)=-R(29);
    X(isp,30)=-R(30);
    X(isp,31)=-R(31);
    X(isp,38)=-R(38);
    X(isp,40)=R(40);
    X(isp,41)=R(41);
    X(isp,49)=R(49);
    X(isp,65)=R(65);
    
    % CO+: species 14
    isp=14;
    X(isp,44)=R(44);
    X(isp,45)=R(45);
    X(isp,46)=-R(46);
    X(isp,47)=R(47);
    X(isp,60)=R(60);
    X(isp,74)=-R(74);
    
    % O3-: species 15
    isp=15;
    X(isp,48)=R(48);
    X(isp,49)=-R(49);
    X(isp,50)=R(50);
    X(isp,51)=R(51);
    X(isp,52)=-R(52);
    X(isp,67)=R(67);

    % O+: species 16
    isp=16;
    X(isp,53)=R(53);
    X(isp,54)=R(54);
    X(isp,55)=-R(55);
    X(isp,56)=-R(56);
    X(isp,57)=R(57);
    X(isp,69)=R(69);
    X(isp,71)=R(71);

    % C+: species 17
    isp=17;
    X(isp,58)=R(58);
    X(isp,59)=R(59);
    X(isp,60)=-R(60);
    X(isp,61)=R(61);
    X(isp,68)=R(68);
    
    % CO4+: species 18
    isp=18;
    X(isp,72)=-R(72);
    X(isp,73)=R(73);
end

%Wall losses

X(:,(length(R)+1))=W';

%Map reactions to wavelengths
%This can be use to generate synthetic spectra
lambda=zeros(1,length(R));
%Example: Reaction 4 radiates at 850nm
lambda(4)= 850; %此处非准确数值

% 能量变化

dS_ape=2/(3*e*n(1))*S_abs*dt;
dS_el=-2*me/mg*K(13)*n(10)*(Te-Tg)*dt;
dS_epro=-Te/n(1).*(sum(X(1,:)))*dt;
dS_rea=-2/3*(13.7*K(1)*n(10)+11.46*K(2)*n(10)+6.12*K(5)*n(3)+9.2*K(17)*n(2)+11.1*K(18)*n(2)+12.06*K(22)*n(3)+19.5*K(44)*n(10)+14.01*K(45)*n(2)+19.1*K(53)*n(10)+13.6*K(54)*n(4)+22*K(58)*n(2)+11.2*K(59)*n(11)+27.8*K(68)*n(10)+25*K(69)*n(2)+19.5*K(71)*n(3))*dt;
dS_wall=-2/(3*n(1))*(-W(1))*(Te/2*log(mg/(2*pi*me))+5/2*Te)*dt;%(Te/2*log(mg/(2*pi*me))+5/2*Te)为1个电子-离子对的能量
dTe=dS_ape+dS_el+dS_epro+dS_rea+dS_wall;

% 建立有向关系
% 建立相互关系矩阵
i=1;
while i<nsp+1
    for j=1:nsp
        B{i,j}=intersect(find(X(i,:)>0),find(X(j,:)<0));
    end
    i=i+1;
end

% % 在该矩阵中对有电子参与碰撞但电子数目增加或维持不变的反应进行添加
B{2,1}=[2,53,B{2,1}];
B{3,1}=[4,68,B{3,1}];
B{4,1}=[2,4,5,18,33,44,58,71,B{4,1}];
B{8,1}=[1,B{8,1}];
B{11,1}=[18,32,69,B{11,1}];
B{12,1}=[22,32,33,B{12,1}];
B{14,1}=[44,45,B{14,1}];
B{16,1}=[53,54,69,71,B{16,1}];
B{17,1}=[58,59,68,B{17,1}];
B{13,2}=[40,B{13,2}];
B{1,3}=[12,76,B{1,3}];
B{4,3}=[12,B{4,3}];
B{6,3}=[75,B{6,3}];
B{9,3}=[62,B{9,3}];
B{13,3}=[41,B{13,3}];
B{15,3}=[48,B{15,3}];
B{18,3}=[73,B{18,3}];
B{1,10}=[15,B{1,10}];
B{2,10}=[81,B{2,10}];
B{3,10}=[42,B{3,10}];
B{4,10}=[15,82,B{4,10}];
B{5,10}=[11,B{5,10}];

% 建立有向矩阵
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

