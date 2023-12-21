function [nsp,species,X,R,lambda,dTe,B,Y,dS_ape,dS_el,dS_epro,dS_rea,dS_wall,number_of_reactions]=rates_CO2m5(Te,P,n,area_volume,Tg,only_e_balance,me,mg,e,S_abs,dt)

% R1:  e + CO2 => CO2^+ + 2e          k= 查表                                                ER=13.7eV
% R2:  e + CO2 => CO + O + e          k= 查表            (激发解离)                          ER=11.46eV
% R3:  e + CO2 => CO + O^-            k= 查表                                            （附着反应无ER）
% R4:  e + O3  => O + O2 + e          k= 1.0e-14*(0.02586/Tg)^(0.5)  (激发解离)      (应该存在ER，但未找到数据)
% R5:  e + O2  => 2O + e              k= 查表            (生成的O一般而言属于激发态)          ER=6.12eV
% R6:  e + O2  => O + O^-             k= 查表                                             (附着反应无ER)
% R7:  e + O2 + M => O2^- + M         k= 3.0e-42                                      （三体反应无法寻找ER）
% R8:  O^- + CO => CO2 + e            k= 5.5e-16
% R9:  O^- + O2 => O3 + e             k= 1.0e-18
% R10: O^- + O3 => 2O2 + e            k= 3.0e-16
% R11: e + CO2^+ => CO + O            k= 4.2e-13*(Te/300*1.16e4)^(-0.75)
% R12: O2^- + CO2^+ => CO + O2 + O    k= 6.0e-13
% R13: O + O + M => O2 + M            k= 5.2e-47*exp(900/Tg)          (K)
% R14: O + O2 + M => O3 + M           k= 4.5e-51*(Tg/298)^-2.70       (K)
% R15: O + O3 => 2O2                  k= 8.0e-18*exp(-17.13/Tg)       (K)                                             
% R16: O + CO + M => CO2 + M          k= 1.7e-45*exp(-1510/Tg)        (K)
% R17: O3 + M => O2 + O + M           k= 4.1e-16*exp(-11430/Tg)       (K)
% R18: e + CO2 => e + CO2             k= 查表                                          （弹性碰撞反应ER=0）
% R19: e + O + M => O^ - + M          k= 1.0e-43                                       (三体反应无法寻找ER)
% R20: O2^- + M => O2 + M + e         k= 2.7e-16*(Tg/300*1.16e4)^(0.5)*exp(-5590/(Tg*1.16e4))
% R21: O + O2^- => O2 + O^-           k= 3.31e-16
% R22: O^- + M => O + M + e           k= 4.0e-18
% R23: O^- + O => O2 + e              k= 2.3e-16
% R24: CO2 + M => CO + O + M          k= 1.81e-16*exp(-49000/(Tg*1.16e4))
% R25: O + CO2 => CO + O2             k= 2.8e-17*exp(-26500/(Tg*1.16e4))
% R26: CO + O2 => CO2 + O             k= 4.2e-18*exp(-24000/(Tg*1.16e4))
% R27: e + CO => C + O^-              k= 查表                                                ER=9.2eV
% R28: e + CO => e + C + O            k= 查表                                               ER=11.1eV
% R29: e + CO2^+ => C + O2            k= 3.94e-13*Te^(-0.4)
% R30: CO2 + C => 2CO                 k= 1.0e-21
% R31: C + O2 => O + CO               k= 3.0e-17
% R32: O + C + M => M + CO            k= 2.14e-41*(Tg/300*1.16e4)^(-3.08)*exp(-2114/(Tg*1.16e4))
% R33: e + O2 => 2e + O2^+            k= 查表                                               ER=12.06eV
% R34: e + O2^+ + M => O2 + M         k= 1.0e-38 
% R35: e + O2^+ => 2O                 k= 6.0e-13*(1/Te)^0.5*(1/(Tg*1.16e4))^0.5
% R36: O^- + O2^+ => O2 + O           k= 2.6e-14*(300/(Tg*1.16e4))^0.44
% R37: O^- + O2^+ => 3O               k= 4.2e-13*exp(300/(Tg*1.16e4))^0.44
% R38: O2^+ + O2^- => 2O2             k= 2.01e-13*(300/(Tg*1.16e4))^0.5
% R39: O2^+ + O2^- => O2 + 2O         k= 4.2e-13
% R40: CO + M => O + C + M            k= 1.52e-10*((Tg*1.16e4)/298)^(-3.1)*exp(-129000/(Tg*1.16e4))
% R41: O2^- + O2 => 2O2 + e           k= 2.18e-24
% R42: O^- + CO2 + M => CO3^- + M     k= 9.0e-41
% R43: CO3^- + CO => 2CO2 + e         k= 5.0e-19
% R44: CO3^- + CO2^+ => 2CO2 + O      k= 5.0e-13
% R45: O + CO3^- => CO2 + O2^-        k= 8.0e-17
% R46: CO2 + e => O2^+ + C + 2e       k= K(1)*1/13
% R47: O3 + e => O2^+ + O + 2e        k= 3.2e-17*(Te*1.16e4)^(0.5)*(1+0.15*Te)*exp(-12.93/Te)
% R48: O3 + e => O2 + O^-             k= 8.9e-18*((Tg*1.16e4)/298)^(1.5)                  （附着反应无ER）
% R49: CO2^+ + O => O2^+ + CO         k= 0.63*2.6e-16
% R50: CO2^+ + O2 => O2^+ + CO2       k= 6.4e-17
% !!(下面数值未更正)R51: O2^- + O => e + O3             k= 3.3e-16
% R52: CO + O3 => CO2 + O2            k= 4.0e-31
% R53: O2 + M => 2O + M               k= 3.0e-12*(Tg*1.16e4)^(-1)*exp(-53980/(Tg*1.16e4))
% !!(下面数值未更正)R54: e + O3 + M => O3^- + M         k= 5.0e-40*Te^(-0.5)
% R55: O^- + O2 + M => O3^- + M       k= 1.1e-42*((Tg*1.16e4)/298)
% R56: O^- + O3 => O3^- + O           k= 8.0e-16
% R57: O2^- + O3 => O3^- + O2         k= 6.0e-16
% R58: O3^- + CO2 => CO3^- + O2       k= 5.5e-16
% R59: O3^- + M => O3 + e + M         k= 2.3e-17
% R60: 2CO2 + O^- => CO3^- + CO2      k= 1.5e-40
% R61: CO2 + CO + O^- => CO3^- + CO   k= 1.5e-40
% R62: CO2 + O2 + O^- => CO3^- + O2   k= 3.1e-40
% R63: CO2 + O2^- + M => CO4^- + M    k= 4.7e-41
% R64: CO4^- + O => CO3^- + O2        k= 1.12e-16
% R65: CO4^- + O => CO2 + O2 + O^-    k= 1.4e-17
% R66: e + CO => 2e + CO^+            k= 查表                                                ER=14.01eV
% R67: e + CO^+ => C + O              k= 3.46e-14*Te^(-0.8)
% R68: CO2 + CO^+ => CO2^+ + CO       k= 1.0e-15
% R69: CO^+ + O2 => O2^+ + CO         k= 1.2e-16
% R70: CO4^- + O3 => O3^- + CO2 + O2  k= 1.3e-16 
% R71: O3^- + O => O2^- + O2          k= 2.5e-16
% R72; CO3^- + O2 => O3^- + CO2       k= 6.0e-21

global B species

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
species{11}='C';        charge(11)=0;       mass(11)=1.99e-26;
species{12}='O_2^+';    charge(12)=1;       mass(12)=5.3e-26;
species{13}='CO_3^-';   charge(13)=-1;      mass(13)=9.97e-26;
species{14}='O_3^-';    charge(14)=-1;      mass(14)=7.97e-26;
species{15}='CO_4^-';   charge(15)=-1;      mass(15)=1.261e-25;
species{16}='CO^+';     charge(16)=1;       mass(16)=4.65e-26;
nsp=length(species);

if(Te==0)
    return;
end

ng=3.54e13*P*1e9;%m-3
 
%Gas phase reactions
number_of_reactions=72;
K(1:number_of_reactions)=0;
R(1:number_of_reactions)=0;
X(nsp,number_of_reactions+1)=0;

%First: Reactions involving generation/loss of electrons      
K(1)=lookups('reaction-rate-cross_section-e_CO2ionization.lut',Te);
K(3)=lookups('reaction-rate-cross_section-e_CO2attachment.lut',Te);
K(6)=lookups('reaction-rate-cross_section-e_O2attachment.lut',Te);
K(7)=3.0e-42;
K(8)=5.5e-16;
K(9)=1.0e-18;
K(10)=3.0e-16;
K(11)=4.2e-13*(Te/300*1.16e4)^(-0.75);
K(19)=1.0e-43;
K(20)=2.7e-16*(Tg/300*1.16e4)^(0.5)*exp(-5590/(Tg*1.16e4));
K(22)=4.0e-18;
K(23)=2.3e-16;
K(27)=lookups('reaction-rate-cross_section-e_COattachment.lut',Te);
K(29)=3.94e-13*Te^(-0.4);
K(33)=lookups('reaction-rate-cross_section-e_O2ionization.lut',Te);
K(34)=1.0e-38;
K(35)=6.0e-13*(1/Te)^0.5*(1/(Tg*1.16e4))^0.5;
K(41)=2.18e-24;
K(43)=5.0e-19;
K(46)=K(1)*1/13;
K(47)=3.2e-17*(Te*1.16e4)^(0.5)*(1+0.15*Te)*exp(-12.93/Te);
K(48)=8.9e-18*((Tg*1.16e4)/298)^(1.5);
K(51)=1.5e-16;
K(54)=4.6e-40;
K(59)=2.3e-17;
K(66)=lookups('reaction-rate-cross_section-e_COionization.lut',Te);
K(67)=3.46e-14*Te^(-0.48);
R(1)=K(1)*n(1)*n(10);
R(3)=K(3)*n(1)*n(10);
R(6)=K(6)*n(1)*n(3);
R(7)=K(7)*n(1)*n(3)*n(9);
R(8)=K(8)*n(6)*n(2);
R(9)=K(9)*n(6)*n(3);
R(10)=K(10)*n(6)*n(5);
R(11)=K(11)*n(1)*n(8);
R(19)=K(19)*n(1)*n(4)*n(9);
R(20)=K(20)*n(7)*n(9);
R(22)=K(22)*n(6)*n(9);
R(23)=K(23)*n(6)*n(4);
R(27)=K(27)*n(1)*n(2);
R(29)=K(29)*n(1)*n(8);
R(33)=K(33)*n(1)*n(3);
R(34)=K(34)*n(1)*n(9)*n(12);
R(35)=K(35)*n(1)*n(12);
R(41)=K(41)*n(3)*n(7);
R(43)=K(43)*n(2)*n(13);
R(46)=K(46)*n(1)*n(10);
R(47)=K(47)*n(1)*n(5);
R(48)=K(48)*n(1)*n(5);
R(51)=K(51)*n(4)*n(7);
R(54)=K(54)*n(1)*n(5)*n(9);
R(59)=K(59)*n(9)*n(14);
R(66)=K(66)*n(1)*n(2);
R(67)=K(67)*n(1)*n(16);

%Second: Rest of reactions
if(~only_e_balance)
    K(2)=lookups('reaction-rate-cross_section-e_CO2dissociation.lut',Te);
    K(4)=2e-15;
    K(5)=lookups('reaction-rate-cross_section-e_O2dissociation.lut',Te);
    K(12)=6.0e-13;
    K(13)=5.2e-47*exp(0.07755/Tg);
    K(14)=4.5e-46*(Tg/0.02568)^(-2.70);
    K(15)=8.0e-18*exp(-0.00148/Tg);
    K(16)=1.7e-45*exp(-0.13011/Tg);
    K(17)=4.1e-16*exp(-0.98492/Tg);
    K(18)=lookups('reaction-rate-cross_section-e_CO2elastic.lut',Te);%e与CO2发生弹性碰撞的速率系数
    K(21)=3.31e-16;
    K(24)=1.81e-16*exp(-49000/(Tg*1.16e4));
    K(25)=2.8e-17*exp(-26500/(Tg*1.16e4));
    K(26)=4.2e-18*exp(-24000/(Tg*1.16e4));
    K(28)=lookups('reaction-rate-cross_section-e_COexcitation.lut',Te);
    K(30)=1.0e-21;
    K(31)=3.0e-17;
    K(32)=2.14e-41*(Tg/300*1.16e4)^(-3.08)*exp(-2114/(Tg*1.16e4));
    K(36)=2.6e-14*(300/(Tg*1.16e4))^0.44;
    K(37)=4.2e-13*exp(300/(Tg*1.16e4))^0.44;
    K(38)=2.01e-13*(300/(Tg*1.16e4))^0.5;
    K(39)=4.2e-13;
    K(40)=1.52e-10*((Tg*1.16e4)/298)^(-3.1)*exp(-129000/(Tg*1.16e4));
    K(42)=9.0e-41;
    K(44)=5.0e-13;
    K(45)=8.0e-17;
    K(49)=0.63*2.6e-16;
    K(50)=6.4e-17;
    K(52)=4.0e-31;
    K(53)=3.0e-12*(Tg*1.16e4)^(-1)*exp(-53980/(Tg*1.16e4));
    K(55)=1.1e-42*((Tg*1.16e4)/298);
    K(56)=8.0e-16;
    K(57)=6.0e-16;
    K(58)=5.5e-16;
    K(60)=1.5e-40;
    K(61)=1.5e-40;
    K(62)=3.1e-40;
    K(63)=4.7e-41;
    K(64)=1.12e-16;
    K(65)=1.4e-17;
    K(68)=1.0e-15;
    K(69)=1.2e-16;
    K(70)=1.3e-16;
    K(71)=2.5e-16;
    K(72)=6.0e-21;
    R(2)=K(2)*n(1)*n(10);
    R(4)=K(4)*n(1)*n(5);
    R(5)=K(5)*n(1)*n(3);
    R(12)=K(12)*n(7)*n(8);
    R(13)=K(13)*n(4)*n(9)*n(4);
    R(14)=K(14)*n(4)*n(3)*n(9);
    R(15)=K(15)*n(4)*n(5);
    R(16)=K(16)*n(4)*n(2)*n(9);
    R(17)=K(17)*n(5)*n(9);
    R(18)=K(18)*n(1)*n(10);
    R(21)=K(21)*n(4)*n(7);
    R(24)=K(24)*n(9)*n(10);
    R(25)=K(25)*n(4)*n(10);
    R(26)=K(26)*n(2)*n(3);
    R(28)=K(28)*n(1)*n(2);
    R(30)=K(30)*n(10)*n(11);
    R(31)=K(31)*n(11)*n(3);
    R(32)=K(32)*n(4)*n(9)*n(11);
    R(36)=K(36)*n(6)*n(12);
    R(37)=K(37)*n(6)*n(12);
    R(38)=K(38)*n(7)*n(12);
    R(39)=K(39)*n(7)*n(12);
    R(40)=K(40)*n(2)*n(9);
    R(42)=K(42)*n(6)*n(9)*n(10);
    R(44)=K(44)*n(8)*n(13);
    R(45)=K(45)*n(4)*n(13);
    R(49)=K(49)*n(4)*n(8);
    R(50)=K(50)*n(3)*n(8);
    R(52)=K(52)*n(2)*n(5);
    R(53)=K(53)*n(3)*n(9);
    R(55)=K(55)*n(3)*n(6)*n(9);
    R(56)=K(56)*n(5)*n(6);
    R(57)=K(57)*n(5)*n(7);
    R(58)=K(58)*n(10)*n(14);
    R(60)=K(60)*n(6)*n(10)*n(10);
    R(61)=K(61)*n(2)*n(6)*n(10);
    R(62)=K(62)*n(3)*n(6)*n(10);
    R(63)=K(63)*n(7)*n(9)*n(10);
    R(64)=K(64)*n(4)*n(15);
    R(65)=K(65)*n(4)*n(15);
    R(68)=K(68)*n(10)*n(16);
    R(69)=K(69)*n(3)*n(16);
    R(70)=K(70)*n(5)*n(15);
    R(71)=K(71)*n(4)*n(14);
    R(72)=K(72)*n(3)*n(13);
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
W(14)=0;
W(15)=0;
W(16)=-1/10*area_volume*n(16)*sqrt(1.6e-19*Te/mass(16));
W(8)=-1/10*area_volume*n(8)*sqrt(1.6e-19*Te/mass(8));
W(12)=-1/10*area_volume*n(12)*sqrt(1.6e-19*Te/mass(12));
W(1)=W(8)+W(12)+W(16); %注意此处需要改成W(1)=W(8)+W(12)+W(16);!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
X(isp,19)=-R(19);
X(isp,20)=R(20);
X(isp,22)=R(22);
X(isp,23)=R(23);
X(isp,27)=-R(27);
X(isp,29)=-R(29);
X(isp,33)=R(33);
X(isp,34)=-R(34);
X(isp,35)=-R(35);
X(isp,41)=R(41);
X(isp,43)=R(43);
X(isp,46)=R(46);
X(isp,47)=R(47);
X(isp,48)=-R(48);
X(isp,51)=R(51);
X(isp,54)=-R(54);
X(isp,59)=R(59);
X(isp,66)=R(66);
X(isp,67)=-R(67);

if(~only_e_balance)
    % CO: species 2
    isp=2;
    X(isp,2)=R(2);
    X(isp,3)=R(3);
    X(isp,8)=-R(8);
    X(isp,11)=R(11);
    X(isp,12)=R(12);
    X(isp,16)=-R(16);
    X(isp,24)=R(24);
    X(isp,25)=R(25);
    X(isp,26)=-R(26);
    X(isp,27)=-R(27);
    X(isp,28)=-R(28);
    X(isp,30)=2*R(30);
    X(isp,31)=R(31);
    X(isp,32)=R(32);
    X(isp,40)=-R(40);
    X(isp,43)=-R(43);
    X(isp,49)=R(49);
    X(isp,52)=-R(52);
    X(isp,66)=-R(66);
    X(isp,68)=R(68);
    X(isp,69)=R(69);
    
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
    X(isp,20)=R(20);
    X(isp,21)=R(21);
    X(isp,23)=R(23);
    X(isp,25)=R(25);
    X(isp,26)=-R(26);
    X(isp,29)=R(29);
    X(isp,31)=-R(31);
    X(isp,33)=-R(33);
    X(isp,34)=R(34);
    X(isp,36)=R(36);
    X(isp,38)=2*R(38);
    X(isp,39)=R(39);
    X(isp,41)=R(41);
    X(isp,48)=R(48);
    X(isp,50)=-R(50);
    X(isp,52)=R(52);
    X(isp,53)=-R(53);
    X(isp,55)=-R(55);
    X(isp,57)=R(57);
    X(isp,58)=R(58);
    X(isp,64)=R(64);
    X(isp,65)=R(65);
    X(isp,69)=-R(69);
    X(isp,70)=R(70);
    X(isp,71)=R(71);
    X(isp,72)=-R(72);
    
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
    X(isp,19)=-R(19);
    X(isp,21)=-R(21);
    X(isp,22)=R(22);
    X(isp,23)=-R(23);
    X(isp,24)=R(24);
    X(isp,25)=-R(25);
    X(isp,26)=R(26);
    X(isp,28)=R(28);
    X(isp,31)=R(31);
    X(isp,32)=-R(32);
    X(isp,35)=2*R(35);
    X(isp,36)=R(36);
    X(isp,37)=3*R(37);
    X(isp,39)=2*R(39);
    X(isp,40)=R(40);
    X(isp,44)=R(44);
    X(isp,45)=-R(45);
    X(isp,47)=R(47);
    X(isp,49)=-R(49);
    X(isp,51)=-R(51);
    X(isp,53)=2*R(53);
    X(isp,56)=R(56);
    X(isp,64)=-R(64);
    X(isp,65)=-R(65);
    X(isp,67)=R(67);
    X(isp,71)=-R(71);
    
    % O3: species 5
    isp=5;
    X(isp,4)=-R(4);
    X(isp,9)=R(9);
    X(isp,10)=-R(10);
    X(isp,14)=R(14);
    X(isp,15)=-R(15);
    X(isp,17)=-R(17);
    X(isp,47)=-R(47);
    X(isp,48)=-R(48);
    X(isp,51)=R(51);
    X(isp,52)=-R(52);
    X(isp,54)=-R(54);
    X(isp,56)=-R(56);
    X(isp,57)=-R(57);
    X(isp,59)=R(59);
    X(isp,70)=-R(70);
    
    % O—: species 6
    isp=6;
    X(isp,3)=R(3);
    X(isp,6)=R(6);
    X(isp,8)=-R(8);
    X(isp,9)=-R(9);
    X(isp,10)=-R(10);
    X(isp,19)=R(19);
    X(isp,21)=R(21);
    X(isp,22)=-R(22);
    X(isp,23)=-R(23);
    X(isp,27)=R(27);
    X(isp,36)=-R(36);
    X(isp,37)=-R(37);
    X(isp,42)=-R(42);
    X(isp,48)=R(48);
    X(isp,55)=-R(55);
    X(isp,56)=-R(56);
    X(isp,60)=-R(60);
    X(isp,61)=-R(61);
    X(isp,62)=-R(62);
    X(isp,65)=R(65);
     
    % O2-: species 7
    isp=7;
    X(isp,7)=R(7);
    X(isp,12)=-R(12);
    X(isp,20)=-R(20);
    X(isp,21)=-R(21);
    X(isp,38)=-R(38);
    X(isp,39)=-R(39);
    X(isp,41)=-R(41);
    X(isp,45)=R(45);
    X(isp,51)=-R(51);
    X(isp,57)=-R(57);
    X(isp,63)=-R(63);
    X(isp,71)=R(71);
    
    % CO2+; species 8
    isp=8;
    X(isp,1)=R(1);
    X(isp,11)=-R(11);
    X(isp,12)=-R(12);
    X(isp,29)=-R(29);
    X(isp,44)=-R(44);
    X(isp,49)=-R(49);
    X(isp,50)=-R(50);
    X(isp,68)=R(68);
    
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
    X(isp,24)=-R(24);
    X(isp,25)=-R(25);
    X(isp,26)=R(26);
    X(isp,30)=-R(30);
    X(isp,42)=-R(42);
    X(isp,43)=2*R(43);
    X(isp,44)=2*R(44);
    X(isp,45)=R(45);
    X(isp,46)=-R(46);
    X(isp,50)=R(50);
    X(isp,52)=R(52);
    X(isp,58)=-R(58);
    X(isp,60)=-R(60);
    X(isp,61)=-R(61);
    X(isp,62)=-R(62);
    X(isp,63)=-R(63);
    X(isp,65)=R(65);
    X(isp,68)=-R(68);
    X(isp,70)=R(70);
    X(isp,72)=R(72);
    
    % C: species 11
    isp=11;
    X(isp,27)=R(27);
    X(isp,28)=R(28);
    X(isp,29)=R(29);
    X(isp,30)=-R(30);
    X(isp,31)=-R(31);
    X(isp,32)=-R(32);
    X(isp,40)=R(40);
    X(isp,46)=R(46);
    X(isp,67)=R(67);
    
    % O2+: species 12
    isp=12;
    X(isp,33)=R(33);
    X(isp,34)=-R(34);
    X(isp,35)=-R(35);
    X(isp,36)=-R(36);
    X(isp,37)=-R(37);
    X(isp,38)=-R(38);
    X(isp,46)=R(46);
    X(isp,47)=R(47);
    X(isp,49)=R(49);
    X(isp,50)=R(50);
    X(isp,69)=R(69);
    
    % CO3-: species 13
    isp=13;
    X(isp,42)=R(42);
    X(isp,43)=-R(43);
    X(isp,44)=-R(44);
    X(isp,45)=-R(45);
    X(isp,58)=R(58);
    X(isp,60)=R(60);
    X(isp,61)=R(61);
    X(isp,62)=R(62);
    X(isp,64)=R(64);
    X(isp,72)=-R(72);
    
    % O3-: species 14
    isp=14;
    X(isp,54)=R(54);
    X(isp,55)=R(55);
    X(isp,56)=R(56);
    X(isp,57)=R(57);
    X(isp,58)=-R(58);
    X(isp,59)=-R(59);
    X(isp,70)=R(70);
    X(isp,71)=-R(71);
    X(isp,72)=R(72);
    
    % CO4-: species 15
    isp=15;
    X(isp,63)=R(63);
    X(isp,64)=-R(64);
    X(isp,65)=-R(65);
    X(isp,70)=-R(70);
    
    % CO+: species 16
    isp=16;
    X(isp,66)=R(66);
    X(isp,67)=-R(67);
    X(isp,68)=-R(68);
    X(isp,69)=-R(69);
    
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
dS_el=-2*me/mg*K(18)*n(10)*(Te-Tg)*dt;
%dS_epro=0;
dS_epro=-Te/n(1).*(sum(X(1,:)))*dt;
%dS_rea=-2/3*(13.7*K(1)*n(10)+6.12*K(5)*n(3)+9.2*K(27)*n(2)+11.1*K(28)*n(2)+12.06*K(33)*n(3))*dt;
dS_rea=-2/3*(13.7*K(1)*n(10)+11.46*K(2)*n(10)+6.12*K(5)*n(3)+9.2*K(27)*n(2)+11.1*K(28)*n(2)+12.06*K(33)*n(3)+14.01*K(66)*n(2))*dt;
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

% 在该矩阵中对有电子参与碰撞但电子数目增加或维持不变的反应进行添加
B{2,1}=[2,B{2,1}];
B{3,1}=[4,B{3,1}];
B{4,1}=[2,4,5,28,47,B{4,1}];
B{8,1}=[1,B{8,1}];
B{11,1}=[28,46,B{11,1}];
B{12,1}=[33,46,47,B{12,1}];
B{16,1}=[66,B{16,1}];
B{13,2}=[61,B{13,2}];
B{1,3}=[41,B{1,3}];
B{13,3}=[62,B{13,3}];
% M的添加仅为findequation.m提供帮助
B{1,9}=[20,22,59];
B{2,9}=[24,32];
B{3,9}=[13,17,20,34];
B{4,9}=[17,22,24,40,53];
B{5,9}=[14,59];
B{6,9}=19;
B{7,9}=7;
B{10,9}=16;
B{11,9}=40;
B{13,9}=42;
B{14,9}=[54,55];
B{15,9}=63;

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

