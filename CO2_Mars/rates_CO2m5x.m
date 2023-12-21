function [nsp,species,X,R,lambda,dTe,B,Y,dS_ape,dS_el,dS_epro,dS_rea,dS_wall,number_of_reactions]=rates_CO2m5x(Te,P,n,area_volume,Tg,only_e_balance,me,mg,e,S_abs,dt)

%R13、14、16、17M的变化见文章Fluid modelling of CO2 dissociation in a dielectric barrier discharge（DBD liuti）
% R1:  e + CO2 => CO2^+ + 2e          k= 查表                                                ER=13.7eV
% R2:  e + CO2 => CO + O + e          k= 查表                   (激发解离)                   ER=11.46eV
% R3:  e + CO2 => CO + O^-            k= 查表                                            （附着反应无ER）
% R4:  e + O3  => O + O2 + e          k= 9.0e-16                (激发解离)      (应该存在ER，但未找到数据)
% R5:  e + O2  => 2O + e              k= 查表            (生成的O一般而言属于激发态)          ER=6.12eV
% R6:  e + O2  => O + O^-             k= 查表                                             (附着反应无ER)
% R7:  e + O2 + M => O2^- + M         k= 3.0e-42                                      （三体反应无法寻找ER）
% R8:  O^- + CO => CO2 + e            k= 5.5e-16
% R9:  O^- + O2 => O3 + e             k= 1.0e-18
% R10: O^- + O3 => 2O2 + e            k= 3.0e-16
% R11: e + CO2^+ => CO + O            k= 2.0e-11*Te^(-0.5)*(Tg*1.16e4)^(-1)
% R12: O2^- + CO2^+ => CO + O2 + O    k= 6.0e-13
% R13: O + O + M => O2 + M            k= 1.27e-44*((Tg*1.16e4)/300)^(-1)*exp(-170/(Tg*1.16e4))M可以为CO2，见60
% R14: O + O2 + M => O3 + M           k= 6.11e-46*(Tg/300)^(-2.60)                         M可以为O2或CO2
% R15: O + O3 => 2O2                  k= 3.1e-20*(Tg*1.16e4)^(0.75)*exp(-1575/(Tg*1.16e4))                                                    
% R16: O + CO + M => CO2 + M          k= 8.2e-46*exp(-1510/(Tg*1.16e4))                 (M可以为CO2，系数k不变)
% R17: O3 + M => O2 + O + M           k= 4.1e-16*exp(-11430/(Tg*1.16e4))       (速率系数是M为O2计算而来)用O2时4.1改为7.26,见64
% R18: e + CO2 => e + CO2             k= 查表                                            （弹性碰撞反应ER=0）
% R19: e + O + M => O^ - + M          k= 1.0e-43                                         (三体反应无法寻找ER)
% R20: O2^- + M => O2 + M + e         k=2.7e-16*(Tg/300*1.16e4)^(0.5)*exp(-5590/(Tg*1.16e4))(M可以为O2，见41)
% R21: O + O2^- => O2 + O^-           k= 3.31e-16
% R22: O^- + M => O + M + e           k= 4.0e-18                                      （M可以为CO2，系数k不变）
% R23: O^- + O => O2 + e              k= 2.3e-16
% R24: CO2 + M => CO + O + M          k= 1.81e-16*exp(-49000/(Tg*1.16e4))                      见R63
% R25: O + CO2 => CO + O2             k= 2.8e-17*exp(-26500/(Tg*1.16e4))
% R26: CO + O2 => CO2 + O             k= 4.2e-18*exp(-24000/(Tg*1.16e4))
% R27: e + CO => C + O^-              k= 查表                                                 ER=9.2eV
% R28: e + CO => e + C + O            k= 查表                                                 ER=11.1eV
% R29: e + CO2^+ => C + O2            k= 3.94e-13*Te^(-0.4)
% R30: CO2 + C => 2CO                 k= 1.0e-21
% R31: C + O2 => O + CO               k= 3.0e-17
% R32: O + C + M => M + CO            k= 2.14e-41*(Tg/300*1.16e4)^(-3.08)*exp(-2114/(Tg*1.16e4))
% R33: e + O2 => 2e + O2^+            k= 查表                                               ER=12.06eV
% R34: e + O2^+ + M => O2 + M         k= 1.0e-38 
% R35: e + O2^+ => 2O                 k= 6.0e-13*(1/Te)^0.5*(1/(Tg*1.16e4))^0.5
% R36: O^- + O2^+ => O2 + O           k= 2.6e-14*(300/(Tg*1.16e4))^0.44
% R37: O^- + O2^+ => 3O               k= 4.2e-13*(300/(Tg*1.16e4))^0.44
% R38: O2^+ + O2^- => 2O2             k= 2.01e-13*(300/(Tg*1.16e4))^0.5
% R39: O2^+ + O2^- => O2 + 2O         k= 4.2e-13
% R40: CO + M => O + C + M            k= 1.52e-10*((Tg*1.16e4)/298)^(-3.1)*exp(-129000/(Tg*1.16e4))
% R41: O2^- + O2 => 2O2 + e           k= 2.18e-24
% R42: O^- + CO2 + M => CO3^- + M     k= 9.0e-41                             （M可以为CO2，系数k不变，为CO和O2时见58、59）
% R43: CO3^- + CO => 2CO2 + e         k= 5.0e-19
% R44: CO3^- + CO2^+ => 2CO2 + O      k= 5.0e-13                                       （生成的CO2为vb振动态)
% R45: O + CO3^- => CO2 + O2^-        k= 8.0e-17
% R46: CO2 + e => O2^+ + C + 2e       k= K(1)*1/13                                            ER未知
% R47: O3 + e => O2^+ + O + 2e        k= 3.2e-17*(Te*1.16e4)^(0.5)*(1+0.15*Te)*exp(-12.93/Te)
% R48: O3 + e => O2 + O^-             k= 查表                                             （附着反应无ER）
% R49: CO2^+ + O => O2^+ + CO         k= 0.63*2.6e-16
% R50: CO2^+ + O2 => O2^+ + CO2       k= 6.4e-17
% R51: O3 + e => O2^- + O             k= 查表                                             （附着反应无ER）
% R52: CO3^- + O2^+ => CO2 + O2 + O   k= 3.0e-13
% R53: CO + O3 => CO2 + O2            k= 4.0e-31
% R54: O2 + M => 2O + M               k= 3.0e-12*(Tg*1.16e4)^(-1)*exp(-53980/(Tg*1.16e4))      (M为Ar)
% R55: O2^- + O => O3 + e             k= 3.3e-16
% R56: O2^- + O2^+ + M => 2O2 + M     k= 2.0e-37
% R57: O^- + O2^+ + M => O3 + M       k= 2.0e-37
% R58: CO2 + O^- + CO => CO3^- + CO   k= 1.5e-40
% R59: CO2 + O^- + O2 => CO3^- + O2   k= 3.1e-40 
% R60: O + O + CO2 => O2 + CO2        k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% R61: O^- + CO2^+ => O + CO2         k= 1.0e-13
% R62: CO + O + CO2 => 2CO2           k= 1.6e-45*exp(-1510/(Tg*1.16e4))
% R63: CO2 + CO2 => CO + O + CO2      k= 3.91e-16*exp(-49430/(Tg*1.16e4))
% R64: O2 + O3 => 2O2 + O             k= 7.26e-16*exp(-11400/(Tg*1.16e4))
% R65: O2 + O + O2 => O3 + O2         k= 8.61e-43*(Tg*1.16e4)^(-1.25)
% R66: O2 + O + O => O3 + O           k= 2.15e-40*exp(345/(Tg*1.16e4))
% R67: e + CO2 => 2e + O + CO^+       k= 查表                                                ER=19.5eV
% R68: e + CO => 2e + CO^+            k= 查表                                                ER=14.01eV
% R69: e + CO^+ => C + O              k= 3.46e-14*Te^(-0.48)
% R70: CO2 + CO^+ => CO2^+ + CO       k= 1.0e-15
% R71: O2^+ + C => CO^+ + O           k= 5.2e-17
% R72: CO^+ + O2 => O2^+ + CO         k= 1.2e-16
% R73: O^- + CO^+ => CO + O           k= 1.0e-13
% R74: CO^+ + O2^- => CO + O2         k= 2.0e-13
% R75: e + O3 + M => O3^- + M         k= 5.0e-43*Te^(-0.5)
% R76: O3^- + CO2 => CO3^- + O2       k= 5.5e-16
% R77: O^- + O2 + M => O3^- + M       k= 3.0e-40*((Tg*1.16e4)/298)^(-1)
% R78: O^- + O3 => O3^- + O           k= 8.0e-16
% R79: O2^- + O3 => O3^- + O2         k= 4.0e-16
% R80: O3^- + O3 => 3O2 + e           k= 3.0e-16
% R81: O3^- + O => O3 + O^-           k= 1.0e-19
% R82: O3^- + O => 2O2 + e            k= 1.0e-19
% R83: O3^- + O => O2^- + O2          k= 2.5e-16
% R84: O3^- + O2^+ => O2 + O3         k= 2.0e-13
% R85: O3^- + O2^+ => 2O + O3         k= 1.0e-13
% R86: O3^- + M => O3 + e + M         k= 2.3e-17                                     (这里的k就是把M按O2算的)
% R87: O3^- + O2 => O2 + O3 + e       k= 2.3e-17
% R88: e + CO2 => 2e + CO + O^+       k= 查表                                                ER=19.1eV
% R89: e + CO => 2e + C + O^+         k= 查表                                                 ER=25eV
% R90: e + O2 => 2e + O + O^+         k= 查表                                                ER=19.5eV
% R91: e + O => 2e + O^+              k= 查表                                                ER=13.6eV
% R92: e + O^+ + M => O + M           k= 1.0e-38
% R93: O^+ + CO2 => O2^+ + CO         k= 9.4e-16
% R94: O^+ + CO2 => O + CO2^+         k= 4.5e-16
% R95: O^+ + CO => O + CO^+           k= 4.9e-18*((Tg*1.16e4)/298)^(0.5)*exp(-4580/(Tg*1.16e4))
% R96: CO2^+ + O => CO2 + O^+         k= 9.62e-17
% R97: CO^+ + O => CO + O^+           k= 1.4e-16
% R98: O^+ + O2 => O + O2^+           k= 1.9e-17*((Tg*1.16e4)/298)^(-0.5)
% R99: O^+ + O3 => O2^+ + O2          k= 1.0e-16
% R100:O^+ + O + M => O2^+ + M        k= 1.0e-41
% R101:O2^- + O^+ + M => O3 + M       k= 2.0e-37
% R102:O2^- + O^+ => O + O2           k= 2.7e-13
% R103:O3^- + O^+ => O3 + O           k= 1.0e-13
% R104:O^- + O^+ => 2O                k= 4.0e-14
% R105:O^- + O^+ + M => O2 + M        k= 2.0e-37
% R106:O3 + e => O^+ + O^- + O + e    k= 8.4e-18*(Te*1.16e4)^(0.5)*(1+0.11*Te)*exp(-18.1/Te)    ER未知
% R107:e + CO2 => 2e + C^+ + O2       k= 查表                                                  ER=27.8eV
% R108:e + CO => 2e + C^+ + O         k= 查表                                                   ER=22eV
% R109:e + C => 2e + C^+              k= 查表                                                  ER=11.2eV
% R110:C^+ + CO2 => CO^+ + CO         k= 1.1e-15
% R111:C^+ + CO => CO^+ + C           k= 5.0e-19
% R112:CO^+ + C => CO + C^+           k= 1.1e-16
% R113:O2^+ + C => C^+ + O2           k= 5.2e-17
% R114:O2 + C^+ => CO + O^+           k= 6.2e-16
% R115:O2 + C^+ => CO^+ + O           k= 3.8e-16
% R116:C^+ + O^- => C + O             k= 5.0e-14
% R117:C^+ + O2^- => C + O2           k= 5.0e-14
% R118:e + CO4^+ => CO2 + O2          k= 1.61e-13*Te^(-0.5)
% R119:O2^+ + CO2 + M => CO4^+ + M    k= 2.3e-41
% R120:O2^- + CO2 + M => CO4^- + M    k= 1.0e-41
% R121:CO4^- + CO2^+ => 2CO2 + O2     k= 5.0e-13
% R122:O2^+ + CO4^- => CO2 + 2O2      k= 3.0e-13
% R123:CO4^- + O => CO3^- + O2        k= 1.1e-16
% R124:CO4^- + O => CO2 + O2 + O^-    k= 1.4e-17
% R125:CO4^- + O => CO2 + O3^-        k= 1.4e-16
% R126:CO4^- + O3 => CO2 + O3^- + O2  k= 1.3e-16
% R127:e + O4^+ => 2O2                k= 2.25e-13*Te^(-0.5)
% R128:C + CO + M => C2O + M          k= 6.5e-44                                        此处的K是M为CO2计算获得
% R129:O + C2O => 2CO                 k= 5.0e-17
% R130:O2 + C2O => CO2 + CO           k= 3.3e-19
% R131:O4^- + CO2 => O2 + CO4^-       k= 4.8e-16
% R132:O2^+ + O2 + M => O4^+ + M      k= 2.4e-42
% R133:O2^- + O2 + M => O4^- + M      k= 3.5e-43
% R134:O4^- + O => O3^- + O2          k= 4.0e-16
% R135:O4^- + O => O^- + 2O2          k= 3.0e-16
% R136:O4^+ + O => O2^+ + O3          k= 3.0e-16
% R137:O4^- + M => O2^- + O2 + M      k= 4.0e-18    
% R138:O4^+ + M => O2^+ + O2 + M      k= 1.73e-19                                         此处的k是M为O2计算而来
% R139:O5^+ + e => O2 + O3            k= 5.0e-12*((Te*1.16e4)/300)^(-0.6)
% R140:CO4^+ + O3 => O5^+ + CO2       k= 1.0e-16
% R141:O5^+ + CO2 => CO4^+ + O3       k= 1.0e-17
% R142:e + C2O2^+ => 2CO              k= 4.0e-13*Te^(-0.34) 
% R143:e + C2O3^+ => CO2 + CO         k= 5.4e-14*Te^(-0.7)
% R144:e + C2O4^+ => 2CO2             k= 2.0e-11*(Te^(-0.5))*(Tg*1.16e4)^(-1)
% R145:CO2^+ + CO2 + M => C2O4^+ + M  k= 3.0e-40
% R146:C2O3^+ + CO => CO2 + C2O2^+    k= 1.1e-15
% R147:C2O4^+ + CO => CO2 + C2O3^+    k= 9.0e-16
% R148:C2O3^+ + CO + M => CO2 + C2O2^+ + M   k= 2.6e-38
% R149:C2O4^+ + CO + M => CO2 + C2O3^+ + M   k= 4.2e-38
% R150:C2O2^+ + O2 => 2CO + O2^+      k= 5.0e-18
% R151:C2O2^+ + M => CO^+ + CO + M    k= 1.0e-18
% R152:C2O2^+ + CO3^- => CO2 + 2CO + O       k= 5.0e-13
% R153:C2O2^+ + CO4^- => CO2 + 2CO + O2      k= 5.0e-13
% R154:C2O2^+ + O2^- => 2CO + O2      k= 6.0e-13
% R155:C2O3^+ + CO3^- => 2CO2 + CO + O       k= 5.0e-13
% R156:C2O3^+ + CO4^- => 2CO2 + CO + O2      k= 5.0e-13
% R157:C2O3^+ + O2^- => CO2 + CO + O2        k= 6.0e-13
% R158:C2O4^+ + M => CO2^+ + CO2 + M  k= 1.0e-20
% R159:C2O4^+ + CO3^- => 3CO2 + O     k= 5.0e-13
% R160:C2O4^+ + CO4^- => 3CO2 + O2    k= 5.0e-13
% R161:C2O4^+ + O2^- => 2CO2 + O2     k= 6.0e-13

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
species{14}='CO^+';     charge(14)=1;       mass(14)=4.65e-26;
species{15}='O_3^-';    charge(15)=-1;      mass(15)=7.97e-26;
species{16}='O^+';      charge(16)=1;       mass(16)=2.66e-26;
species{17}='C^+';      charge(17)=1;       mass(17)=1.99e-26;
species{18}='CO_4^+';   charge(18)=1;       mass(18)=1.261e-25;
species{19}='CO_4^-';   charge(19)=-1;      mass(19)=1.261e-25;
species{20}='O_4^+';    charge(20)=1;       mass(20)=1.063e-25;
species{21}='C_2O';     charge(21)=0;       mass(21)=6.65e-26;
species{22}='O_4^-';    charge(22)=-1;      mass(22)=1.063e-25;
species{23}='O_5^+';    charge(23)=1;       mass(23)=1.329e-25;
species{24}='C_2O_2^+'; charge(24)=1;       mass(24)=9.3e-26;
species{25}='C_2O_3^+'; charge(25)=1;       mass(25)=1.196e-25;
species{26}='C_2O_4^+'; charge(26)=1;       mass(26)=1.461e-25;

nsp=length(species);

if(Te==0)
    return;
end

ng=3.54e13*P*1e9;%m-3
 
%Gas phase reactions
number_of_reactions=161;
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
K(11)=2.0e-11*Te^(-0.5)*(Tg*1.16e4)^(-1);
K(19)=1.0e-43;
K(20)=2.7e-16*((Tg*1.16e4)/300)^(0.5)*exp(-5590/(Tg*1.16e4));
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
K(48)=lookups('reaction-rate-cross_section-e_O3attachment_O-.lut',Te);
K(51)=lookups('reaction-rate-cross_section-e_O3attachment_O2-.lut',Te);
K(55)=3.3e-16;
K(67)=lookups('reaction-rate-cross_section-e_CO2ionization_CO+.lut',Te);
K(68)=lookups('reaction-rate-cross_section-e_COionization.lut',Te);
K(69)=3.46e-14*Te^(-0.48);
K(75)=5.0e-43*Te^(-0.5);
K(80)=3.0e-16;
K(82)=1.0e-19;
K(86)=2.3e-17;
K(87)=2.3e-17;
K(88)=lookups('reaction-rate-cross_section-e_CO2ionization_O+.lut',Te);
K(89)=lookups('reaction-rate-cross_section-e_COionization_O+.lut',Te);
K(90)=lookups('reaction-rate-cross_section-e_O2ionization_O+.lut',Te);
K(91)=lookups('reaction-rate-cross_section-e_Oionization.lut',Te);
K(92)=1.0e-38;
K(107)=lookups('reaction-rate-cross_section-e_CO2ionization_C+.lut',Te);
K(108)=lookups('reaction-rate-cross_section-e_COionization_C+.lut',Te);
K(109)=lookups('reaction-rate-cross_section-e_Cionization.lut',Te);
K(118)=1.61e-13*Te^(-0.5);
K(127)=2.25e-13*Te^(-0.5);
K(139)=5.0e-12*((Te*1.16e4)/300)^(-0.6);
K(142)=4.0e-13*Te^(-0.34);
K(143)=5.4e-14*Te^(-0.7);
K(144)=2.0e-11*(Te^(-0.5))*(Tg*1.16e4)^(-1);
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
R(51)=K(51)*n(1)*n(5);
R(55)=K(55)*n(4)*n(7);
R(67)=K(67)*n(1)*n(10);
R(68)=K(68)*n(1)*n(2);
R(69)=K(69)*n(1)*n(14);
R(75)=K(75)*n(1)*n(5)*n(9);
R(80)=K(80)*n(5)*n(15);
R(82)=K(82)*n(4)*n(15);
R(86)=K(86)*n(9)*n(15);
R(87)=K(87)*n(3)*n(15);
R(88)=K(88)*n(1)*n(10);
R(89)=K(89)*n(1)*n(2);
R(90)=K(90)*n(1)*n(3);
R(91)=K(91)*n(1)*n(4);
R(92)=K(92)*n(1)*n(9)*n(16);
R(107)=K(107)*n(1)*n(10);
R(108)=K(108)*n(1)*n(2);
R(109)=K(109)*n(1)*n(11);
R(118)=K(118)*n(1)*n(18);
R(127)=K(127)*n(1)*n(20);
R(139)=K(139)*n(1)*n(23);
R(142)=K(142)*n(1)*n(24);
R(143)=K(143)*n(1)*n(25);
R(144)=K(144)*n(1)*n(26);

%Second: Rest of reactions
if(~only_e_balance)
    K(2)=lookups('reaction-rate-cross_section-e_CO2dissociation.lut',Te);
    K(4)=9.0e-16;
    K(5)=lookups('reaction-rate-cross_section-e_O2dissociation.lut',Te);
    K(12)=6.0e-13;
    K(13)=1.27e-44*((Tg*1.16e4)/300)^(-1)*exp(-170/(Tg*1.16e4));
    K(14)=6.11e-46*(Tg/300)^(-2.60);
    K(15)=3.1e-20*(Tg*1.16e4)^(0.75)*exp(-1575/(Tg*1.16e4));
    K(16)=8.2e-46*exp(-1510/(Tg*1.16e4));
    K(17)=4.1e-16*exp(-11430/(Tg*1.16e4));
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
    K(37)=4.2e-13*(300/(Tg*1.16e4))^0.44;
    K(38)=2.01e-13*(300/(Tg*1.16e4))^0.5;
    K(39)=4.2e-13;
    K(40)=1.52e-10*((Tg*1.16e4)/298)^(-3.1)*exp(-129000/(Tg*1.16e4));
    K(42)=9.0e-41;
    K(44)=5.0e-13;
    K(45)=8.0e-17;
    K(49)=0.63*2.6e-16;
    K(50)=6.4e-17;
    K(52)=3.0e-13;
    K(53)=4.0e-31;
    K(54)=3.0e-12*(Tg*1.16e4)^(-1)*exp(-53980/(Tg*1.16e4));
    K(56)=2.0e-37;
    K(57)=2.0e-37;
    K(58)=1.5e-40;
    K(59)=3.1e-40;
    K(60)=3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4));
    K(61)=1.0e-13;
    K(62)=1.6e-45*exp(-1510/(Tg*1.16e4));
    K(63)=3.91e-16*exp(-49430/(Tg*1.16e4));
    K(64)=7.26e-16*exp(-11400/(Tg*1.16e4));
    K(65)=8.61e-43*(Tg*1.16e4)^(-1.25);
    K(66)=2.15e-40*exp(345/(Tg*1.16e4));
    K(70)=1.0e-15;
    K(71)=5.2e-17;
    K(72)=1.2e-16;
    K(73)=1.0e-13;
    K(74)=2.0e-13;
    K(76)=5.5e-16;
    K(77)=3.0e-40*((Tg*1.16e4)/298)^(-1);
    K(78)=8.0e-16;
    K(79)=4.0e-16;
    K(81)=1.0e-19;
    K(83)=2.5e-16;
    K(84)=2.0e-13;
    K(85)=1.0e-13;
    K(93)=9.4e-16;
    K(94)=4.5e-16;
    K(95)=4.9e-18*((Tg*1.16e4)/298)^(0.5)*exp(-4580/(Tg*1.16e4));
    K(96)=9.62e-17;
    K(97)=1.4e-16;
    K(98)=1.9e-17*((Tg*1.16e4)/298)^(-0.5);
    K(99)=1.0e-16;
    K(100)=1.0e-41;
    K(101)=2.0e-37;
    K(102)=2.7e-13;
    K(103)=1.0e-13;
    K(104)=4.0e-14;
    K(105)=2.0e-37;
    K(106)=8.4e-18*(Te*1.16e4)^(0.5)*(1+0.11*Te)*exp(-18.1/Te);
    K(110)=1.1e-15;
    K(111)=5.0e-19;
    K(112)=1.1e-16;
    K(113)=5.2e-17;
    K(114)=6.2e-16;
    K(115)=3.8e-16;
    K(116)=5.0e-14;
    K(117)=5.0e-14;
    K(119)=2.3e-41;
    K(120)=1.0e-41;
    K(121)=5.0e-13;
    K(122)=3.0e-13;
    K(123)=1.1e-16;
    K(124)=1.4e-17;
    K(125)=1.4e-16;
    K(126)=1.3e-16;
    K(128)=6.5e-44;
    K(129)=5.0e-17;
    K(130)=3.3e-19;
    K(131)=4.8e-16;
    K(132)=2.4e-42;
    K(133)=3.5e-43;
    K(134)=4.0e-16;
    K(135)=3.0e-16;
    K(136)=3.0e-16;
    K(137)=4.0e-18;
    K(138)=1.73e-19;
    K(140)=1.0e-16;
    K(141)=1.0e-17;
    K(145)=3.0e-40;
    K(146)=1.1e-15;
    K(147)=9.0e-16;
    K(148)=2.6e-38;
    K(149)=4.2e-38;
    K(150)=5.0e-18;
    K(151)=1.0e-18;
    K(152)=5.0e-13;
    K(153)=5.0e-13;
    K(154)=6.0e-13;
    K(155)=5.0e-13;
    K(156)=5.0e-13;
    K(157)=6.0e-13;
    K(158)=1.0e-20;
    K(159)=5.0e-13;
    K(160)=5.0e-13;
    K(161)=6.0e-13;
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
    R(52)=K(52)*n(12)*n(13);
    R(53)=K(53)*n(2)*n(5);
    R(54)=K(54)*n(3)*n(9);
    R(56)=K(56)*n(7)*n(9)*n(12);
    R(57)=K(57)*n(6)*n(9)*n(12);
    R(58)=K(58)*n(2)*n(6)*n(10);
    R(59)=K(59)*n(3)*n(6)*n(10);
    R(60)=K(60)*n(4)*n(4)*n(10);
    R(61)=K(61)*n(6)*n(8);
    R(62)=K(62)*n(2)*n(4)*n(10);
    R(63)=K(63)*n(10)*n(10);
    R(64)=K(64)*n(3)*n(5);
    R(65)=K(65)*n(3)*n(3)*n(4);
    R(66)=K(66)*n(3)*n(4)*n(4);
    R(70)=K(70)*n(10)*n(14);
    R(71)=K(71)*n(11)*n(12);
    R(72)=K(72)*n(3)*n(14);
    R(73)=K(73)*n(6)*n(14);
    R(74)=K(74)*n(7)*n(14);
    R(76)=K(76)*n(10)*n(15);
    R(77)=K(77)*n(3)*n(6)*n(9);
    R(78)=K(78)*n(5)*n(6);
    R(79)=K(79)*n(5)*n(7);
    R(81)=K(81)*n(4)*n(15);
    R(83)=K(83)*n(4)*n(15);
    R(84)=K(84)*n(12)*n(15);
    R(85)=K(85)*n(12)*n(15);
    R(93)=K(93)*n(10)*n(16);
    R(94)=K(94)*n(10)*n(16);
    R(95)=K(95)*n(2)*n(16);
    R(96)=K(96)*n(4)*n(8);
    R(97)=K(97)*n(4)*n(14);
    R(98)=K(98)*n(3)*n(16);
    R(99)=K(99)*n(5)*n(16);
    R(100)=K(100)*n(4)*n(9)*n(16);
    R(101)=K(101)*n(7)*n(9)*n(16);
    R(102)=K(102)*n(7)*n(16);
    R(103)=K(103)*n(15)*n(16);
    R(104)=K(104)*n(6)*n(16);
    R(105)=K(105)*n(6)*n(9)*n(16);
    R(106)=K(106)*n(1)*n(5);
    R(110)=K(110)*n(10)*n(17);
    R(111)=K(111)*n(2)*n(17);
    R(112)=K(112)*n(11)*n(14);
    R(113)=K(113)*n(11)*n(12);
    R(114)=K(114)*n(3)*n(17);
    R(115)=K(115)*n(3)*n(17);
    R(116)=K(116)*n(6)*n(17);
    R(117)=K(117)*n(7)*n(17);
    R(119)=K(119)*n(9)*n(10)*n(12);
    R(120)=K(120)*n(7)*n(9)*n(10);
    R(121)=K(121)*n(8)*n(19);
    R(122)=K(122)*n(12)*n(19);
    R(123)=K(123)*n(4)*n(19);
    R(124)=K(124)*n(4)*n(19);
    R(125)=K(125)*n(4)*n(19);
    R(126)=K(126)*n(5)*n(19);
    R(128)=K(128)*n(2)*n(9)*n(11);
    R(129)=K(129)*n(4)*n(21);
    R(130)=K(130)*n(3)*n(21);
    R(131)=K(131)*n(10)*n(22);
    R(132)=K(132)*n(3)*n(9)*n(12);
    R(133)=K(133)*n(3)*n(7)*n(9);
    R(134)=K(134)*n(4)*n(22);
    R(135)=K(135)*n(4)*n(22);
    R(136)=K(136)*n(4)*n(20);
    R(137)=K(137)*n(9)*n(22);
    R(138)=K(138)*n(9)*n(20);
    R(140)=K(140)*n(5)*n(18);
    R(141)=K(141)*n(10)*n(23);
    R(145)=K(145)*n(8)*n(9)*n(10);
    R(146)=K(146)*n(2)*n(25);
    R(147)=K(147)*n(2)*n(26);
    R(148)=K(148)*n(2)*n(9)*n(25);
    R(149)=K(149)*n(2)*n(9)*n(26);
    R(150)=K(150)*n(3)*n(24);
    R(151)=K(151)*n(9)*n(24);
    R(152)=K(152)*n(13)*n(24);
    R(153)=K(153)*n(19)*n(24);
    R(154)=K(154)*n(7)*n(24);
    R(155)=K(155)*n(13)*n(25);
    R(156)=K(156)*n(19)*n(25);
    R(157)=K(157)*n(7)*n(25);
    R(158)=K(158)*n(9)*n(26);
    R(159)=K(159)*n(13)*n(26);
    R(160)=K(160)*n(19)*n(26);
    R(161)=K(161)*n(7)*n(26);
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
W(19)=0;
W(21)=0;
W(22)=0;
W(8)=-1/10*area_volume*n(8)*sqrt(1.6e-19*Te/mass(8));
W(12)=-1/10*area_volume*n(12)*sqrt(1.6e-19*Te/mass(12));
W(14)=-1/10*area_volume*n(14)*sqrt(1.6e-19*Te/mass(14));
W(16)=-1/10*area_volume*n(16)*sqrt(1.6e-19*Te/mass(16));
W(17)=-1/10*area_volume*n(17)*sqrt(1.6e-19*Te/mass(17));
W(18)=-1/10*area_volume*n(18)*sqrt(1.6e-19*Te/mass(18));
W(20)=-1/10*area_volume*n(20)*sqrt(1.6e-19*Te/mass(20));
W(23)=-1/10*area_volume*n(23)*sqrt(1.6e-19*Te/mass(23));
W(24)=-1/10*area_volume*n(24)*sqrt(1.6e-19*Te/mass(24));
W(25)=-1/10*area_volume*n(25)*sqrt(1.6e-19*Te/mass(25));
W(26)=-1/10*area_volume*n(26)*sqrt(1.6e-19*Te/mass(26));
W(1)=W(8)+W(12)+W(14)+W(16)+W(17)+W(18)+W(20)+W(23)+W(24)+W(25)+W(26);

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
X(isp,51)=-R(51);
X(isp,55)=R(55);
X(isp,67)=R(67);
X(isp,68)=R(68);
X(isp,69)=-R(69);
X(isp,75)=-R(75);
X(isp,80)=R(80);
X(isp,82)=R(82);
X(isp,86)=R(86);
X(isp,87)=R(87);
X(isp,88)=R(88);
X(isp,89)=R(89);
X(isp,90)=R(90);
X(isp,91)=R(91);
X(isp,92)=-R(92);
X(isp,107)=R(107);
X(isp,108)=R(108);
X(isp,109)=R(109);
X(isp,118)=-R(118);
X(isp,127)=-R(127);
X(isp,139)=-R(139);
X(isp,142)=-R(142);
X(isp,143)=-R(143);
X(isp,144)=-R(144);

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
    X(isp,53)=-R(53);
    X(isp,62)=-R(62);
    X(isp,63)=R(63);
    X(isp,68)=-R(68);
    X(isp,70)=R(70);
    X(isp,72)=R(72);
    X(isp,73)=R(73);
    X(isp,74)=R(74);
    X(isp,88)=R(88);
    X(isp,89)=-R(89);
    X(isp,93)=R(93);
    X(isp,95)=-R(95);
    X(isp,97)=R(97);
    X(isp,108)=-R(108);
    X(isp,110)=R(110);
    X(isp,111)=-R(111);
    X(isp,112)=R(112);
    X(isp,114)=R(114);
    X(isp,128)=-R(128);
    X(isp,129)=2*R(129);
    X(isp,130)=R(130);
    X(isp,142)=2*R(142);
    X(isp,143)=R(143);
    X(isp,146)=-R(146);
    X(isp,147)=-R(147);
    X(isp,148)=-R(148);
    X(isp,149)=-R(149);
    X(isp,150)=2*R(150);
    X(isp,151)=R(151);
    X(isp,152)=2*R(152);
    X(isp,153)=2*R(153);
    X(isp,154)=2*R(154);
    X(isp,155)=R(155);
    X(isp,156)=R(156);
    X(isp,157)=R(157);
    
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
    X(isp,53)=R(53);
    X(isp,54)=-R(54);
    X(isp,56)=2*R(56);
    X(isp,60)=R(60);
    X(isp,64)=R(64);
    X(isp,65)=-R(65);
    X(isp,66)=-R(66);
    X(isp,72)=-R(72);
    X(isp,74)=R(74);
    X(isp,76)=R(76);
    X(isp,77)=-R(77);
    X(isp,79)=R(79);
    X(isp,80)=3*R(80);
    X(isp,82)=2*R(82);
    X(isp,83)=R(83);
    X(isp,84)=R(84);
    X(isp,90)=-R(90);
    X(isp,98)=-R(98);
    X(isp,99)=R(99);
    X(isp,102)=R(102);
    X(isp,105)=R(105);
    X(isp,107)=R(107);
    X(isp,113)=R(113);
    X(isp,114)=-R(114);
    X(isp,115)=-R(115);
    X(isp,117)=R(117);
    X(isp,118)=R(118);
    X(isp,121)=R(121);
    X(isp,122)=2*R(122);
    X(isp,123)=R(123);
    X(isp,124)=R(124);
    X(isp,126)=R(126);
    X(isp,127)=2*R(127);
    X(isp,130)=-R(130);
    X(isp,131)=R(131);
    X(isp,132)=-R(132);
    X(isp,133)=-R(133);
    X(isp,134)=R(134);
    X(isp,135)=2*R(135);
    X(isp,137)=R(137);
    X(isp,138)=R(138);
    X(isp,139)=R(139);
    X(isp,150)=-R(150);
    X(isp,153)=R(153);
    X(isp,154)=R(154);
    X(isp,156)=R(156);
    X(isp,157)=R(157);
    X(isp,160)=R(160);
    X(isp,161)=R(161);
    
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
    X(isp,51)=R(51);
    X(isp,52)=R(52);
    X(isp,54)=2*R(54);
    X(isp,55)=-R(55);
    X(isp,60)=-2*R(60);
    X(isp,61)=R(61);
    X(isp,62)=-R(62);
    X(isp,63)=R(63);
    X(isp,64)=R(64);
    X(isp,65)=-R(65);
    X(isp,66)=-R(66);
    X(isp,67)=R(67);
    X(isp,69)=R(69);
    X(isp,71)=R(71);
    X(isp,73)=R(73);
    X(isp,78)=R(78);
    X(isp,81)=-R(81);
    X(isp,82)=-R(82);
    X(isp,83)=-R(83);
    X(isp,85)=2*R(85);
    X(isp,90)=R(90);
    X(isp,91)=-R(91);
    X(isp,92)=R(92);
    X(isp,94)=R(94);
    X(isp,95)=R(95);
    X(isp,96)=-R(96);
    X(isp,97)=-R(97);
    X(isp,98)=R(98);
    X(isp,100)=-R(100);
    X(isp,102)=R(102);
    X(isp,103)=R(103);
    X(isp,104)=2*R(104);
    X(isp,106)=R(106);
    X(isp,108)=R(108);
    X(isp,115)=R(115);
    X(isp,116)=R(116);
    X(isp,123)=-R(123);
    X(isp,124)=-R(124);
    X(isp,125)=-R(125);
    X(isp,129)=-R(129);
    X(isp,134)=-R(134);
    X(isp,135)=-R(135);
    X(isp,136)=-R(136);
    X(isp,152)=R(152);
    X(isp,155)=R(155);
    X(isp,159)=R(159);
    
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
    X(isp,51)=-R(51);
    X(isp,53)=-R(53);
    X(isp,55)=R(55);
    X(isp,57)=R(57);
    X(isp,64)=-R(64);
    X(isp,65)=R(65);
    X(isp,66)=R(66);
    X(isp,75)=-R(75);
    X(isp,78)=-R(78);
    X(isp,79)=-R(79);
    X(isp,80)=-R(80);
    X(isp,81)=R(81);
    X(isp,84)=R(84);
    X(isp,85)=R(85);
    X(isp,86)=R(86);
    X(isp,87)=R(87);
    X(isp,99)=-R(99);
    X(isp,101)=R(101);
    X(isp,103)=R(103);
    X(isp,106)=-R(106);
    X(isp,126)=-R(126);
    X(isp,136)=R(136);
    X(isp,139)=R(139);
    X(isp,140)=-R(140);
    X(isp,141)=R(141);
    
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
    X(isp,57)=-R(57);
    X(isp,58)=-R(58);
    X(isp,59)=-R(59);
    X(isp,61)=-R(61);
    X(isp,73)=-R(73);
    X(isp,77)=-R(77);
    X(isp,78)=-R(78);
    X(isp,81)=R(81);
    X(isp,104)=-R(104);
    X(isp,105)=-R(105);
    X(isp,106)=R(106);
    X(isp,116)=-R(116);
    X(isp,124)=R(124);
    X(isp,135)=R(135);
     
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
    X(isp,51)=R(51);
    X(isp,55)=-R(55);
    X(isp,56)=-R(56);
    X(isp,74)=-R(74);
    X(isp,79)=-R(79);
    X(isp,83)=R(83);
    X(isp,101)=-R(101);
    X(isp,102)=-R(102);
    X(isp,117)=-R(117);
    X(isp,120)=-R(120);
    X(isp,133)=-R(133);
    X(isp,137)=R(137);
    X(isp,154)=-R(154);
    X(isp,157)=-R(157);
    X(isp,161)=-R(161);
    
    % CO2+; species 8
    isp=8;
    X(isp,1)=R(1);
    X(isp,11)=-R(11);
    X(isp,12)=-R(12);
    X(isp,29)=-R(29);
    X(isp,44)=-R(44);
    X(isp,49)=-R(49);
    X(isp,50)=-R(50);
    X(isp,61)=-R(61);
    X(isp,70)=R(70);
    X(isp,94)=R(94);
    X(isp,96)=-R(96);
    X(isp,121)=-R(121);
    X(isp,145)=-R(145);
    X(isp,158)=R(158);
    
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
    X(isp,53)=R(53);
    X(isp,58)=-R(58);
    X(isp,59)=-R(59);
    X(isp,61)=R(61);
    X(isp,62)=R(62);
    X(isp,63)=-R(63);
    X(isp,67)=-R(67);
    X(isp,70)=-R(70);
    X(isp,76)=-R(76);
    X(isp,88)=-R(88);
    X(isp,93)=-R(93);
    X(isp,94)=-R(94);
    X(isp,96)=R(96);
    X(isp,107)=-R(107);
    X(isp,110)=-R(110);
    X(isp,118)=R(118);
    X(isp,119)=-R(119);
    X(isp,120)=-R(120);
    X(isp,121)=2*R(121);
    X(isp,122)=R(122);
    X(isp,124)=R(124);
    X(isp,125)=R(125);
    X(isp,126)=R(126);
    X(isp,130)=R(130);
    X(isp,131)=-R(131);
    X(isp,140)=R(140);
    X(isp,141)=-R(141);
    X(isp,143)=R(143);
    X(isp,144)=2*R(144);
    X(isp,145)=-R(145);
    X(isp,146)=R(146);
    X(isp,147)=R(147);
    X(isp,148)=R(148);
    X(isp,149)=R(149);
    X(isp,152)=R(152);
    X(isp,153)=R(153);
    X(isp,155)=2*R(155);
    X(isp,156)=2*R(156);
    X(isp,157)=R(157);
    X(isp,158)=R(158);
    X(isp,159)=3*R(159);
    X(isp,160)=3*R(160);
    X(isp,161)=2*R(161);
    
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
    X(isp,69)=R(69);
    X(isp,71)=-R(71);
    X(isp,89)=R(89);
    X(isp,109)=-R(109);
    X(isp,111)=R(111);
    X(isp,112)=-R(112);
    X(isp,113)=-R(113);
    X(isp,116)=R(116);
    X(isp,117)=R(117);
    X(isp,128)=-R(128);
    
    % O2+: species 12
    isp=12;
    X(isp,33)=R(33);
    X(isp,34)=-R(34);
    X(isp,35)=-R(35);
    X(isp,36)=-R(36);
    X(isp,37)=-R(37);
    X(isp,38)=-R(38);
    X(isp,39)=-R(39);
    X(isp,46)=R(46);
    X(isp,47)=R(47);
    X(isp,49)=R(49);
    X(isp,50)=R(50);
    X(isp,52)=-R(52);
    X(isp,56)=-R(56);
    X(isp,57)=-R(57);
    X(isp,71)=-R(71);
    X(isp,72)=R(72);
    X(isp,84)=-R(84);
    X(isp,85)=-R(85);
    X(isp,93)=R(93);
    X(isp,98)=R(98);
    X(isp,99)=R(99);
    X(isp,100)=R(100);
    X(isp,113)=-R(113);
    X(isp,119)=-R(119);
    X(isp,122)=-R(122);
    X(isp,132)=-R(132);
    X(isp,136)=R(136);
    X(isp,138)=R(138);
    X(isp,150)=R(150);
    
    % CO3-: species 13
    isp=13;
    X(isp,42)=R(42);
    X(isp,43)=-R(43);
    X(isp,44)=-R(44);
    X(isp,45)=-R(45);
    X(isp,52)=-R(52);
    X(isp,58)=R(58);
    X(isp,59)=R(59);
    X(isp,76)=R(76);
    X(isp,123)=R(123);
    X(isp,152)=-R(152);
    X(isp,155)=-R(155);
    X(isp,159)=-R(159);
    
    % CO+: species 14
    isp=14;
    X(isp,67)=R(67);
    X(isp,68)=R(68);
    X(isp,69)=-R(69);
    X(isp,70)=-R(70);
    X(isp,71)=R(71);
    X(isp,72)=-R(72);
    X(isp,73)=-R(73);
    X(isp,74)=-R(74);
    X(isp,95)=R(95);
    X(isp,97)=-R(97);
    X(isp,110)=R(110);
    X(isp,111)=R(111);
    X(isp,112)=-R(112);
    X(isp,115)=R(115);
    X(isp,151)=R(151);
    
    % O3-: species 15
    isp=15;
    X(isp,75)=R(75);
    X(isp,76)=-R(76);
    X(isp,77)=R(77);
    X(isp,78)=R(78);
    X(isp,79)=R(79);
    X(isp,80)=-R(80);
    X(isp,81)=-R(81);
    X(isp,82)=-R(82);
    X(isp,83)=-R(83);
    X(isp,84)=-R(84);
    X(isp,85)=-R(85);
    X(isp,86)=-R(86);
    X(isp,87)=-R(87);
    X(isp,103)=-R(103);
    X(isp,125)=R(125);
    X(isp,126)=R(126);
    X(isp,134)=R(134);

    % O+: species 16
    isp=16;
    X(isp,88)=R(88);
    X(isp,89)=R(89);
    X(isp,90)=R(90);
    X(isp,91)=R(91);
    X(isp,92)=-R(92);
    X(isp,93)=-R(93);
    X(isp,94)=-R(94);
    X(isp,95)=-R(95);
    X(isp,96)=R(96);
    X(isp,97)=R(97);
    X(isp,98)=-R(98);
    X(isp,99)=-R(99);
    X(isp,100)=-R(100);
    X(isp,101)=-R(101);
    X(isp,102)=-R(102);
    X(isp,103)=-R(103);
    X(isp,104)=-R(104);
    X(isp,105)=-R(105);
    X(isp,106)=R(106);
    X(isp,114)=R(114);

    % C+: species 17
    isp=17;
    X(isp,107)=R(107);
    X(isp,108)=R(108);
    X(isp,109)=R(109);
    X(isp,110)=-R(110);
    X(isp,111)=-R(111);
    X(isp,112)=R(112);
    X(isp,113)=R(113);
    X(isp,114)=-R(114);
    X(isp,115)=-R(115);
    X(isp,116)=-R(116);
    X(isp,117)=-R(117);

    % CO4+: species 18
    isp=18;
    X(isp,118)=-R(118);
    X(isp,119)=R(119);
    X(isp,140)=-R(140);
    X(isp,141)=R(141);

    % CO4-: species 19
    isp=19;
    X(isp,120)=R(120);
    X(isp,121)=-R(121);
    X(isp,122)=-R(122);
    X(isp,123)=-R(123);
    X(isp,124)=-R(124);
    X(isp,125)=-R(125);
    X(isp,126)=-R(126);
    X(isp,131)=R(131);
    X(isp,153)=-R(153);
    X(isp,156)=-R(156);
    X(isp,160)=-R(160);

    % O4+: species 20
    isp=20;
    X(isp,127)=-R(127);
    X(isp,132)=R(132);
    X(isp,136)=-R(136);
    X(isp,138)=-R(138);

    % C2O: species 21
    isp=21;
    X(isp,128)=R(128);
    X(isp,129)=-R(129);
    X(isp,130)=-R(130);
    
    % O4-: species 22
    isp=22;
    X(isp,131)=-R(131);
    X(isp,133)=R(133);
    X(isp,134)=-R(134);
    X(isp,135)=-R(135);
    X(isp,137)=-R(137);

    % O5+: species 23
    isp=23;
    X(isp,139)=-R(139);
    X(isp,140)=R(140);
    X(isp,141)=-R(141);

    % C2O2+: species 24
    isp=24;
    X(isp,142)=-R(142);
    X(isp,146)=R(146);
    X(isp,148)=R(148);
    X(isp,150)=-R(150);
    X(isp,151)=-R(151);
    X(isp,152)=-R(152);
    X(isp,153)=-R(153);
    X(isp,154)=-R(154);

    % C2O3+: species 25
    isp=25;
    X(isp,143)=-R(143);
    X(isp,146)=-R(146);
    X(isp,147)=R(147);
    X(isp,148)=-R(148);
    X(isp,149)=R(149);
    X(isp,155)=-R(155);
    X(isp,156)=-R(156);
    X(isp,157)=-R(157);

    % C2O4+: species 26
    isp=26;
    X(isp,144)=-R(144);
    X(isp,145)=R(145);
    X(isp,147)=-R(147);
    X(isp,149)=-R(149);
    X(isp,158)=-R(158);
    X(isp,159)=-R(159);
    X(isp,160)=-R(160);
    X(isp,161)=-R(161);
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
dS_epro=-Te/n(1).*(sum(X(1,:)))*dt;
dS_rea=-2/3*(13.7*K(1)*n(10)+11.46*K(2)*n(10)+6.12*K(5)*n(3)+9.2*K(27)*n(2)+11.1*K(28)*n(2)+12.06*K(33)*n(3)+19.5*K(67)*n(10)+14.01*K(68)*n(2)+19.1*K(88)*n(10)+25*K(89)*n(2)+19.5*K(90)*n(3)+13.6*K(91)*n(4)+27.8*K(107)*n(10)+22*K(108)*n(2)+11.2*K(109)*n(11))*dt;
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
B{2,1}=[2,88,B{2,1}];
B{3,1}=[4,107,B{3,1}];
B{4,1}=[2,4,5,28,47,67,90,106,108,B{4,1}];
B{6,1}=[106,B{6,1}];
B{8,1}=[1,B{8,1}];
B{11,1}=[28,46,89,B{11,1}];
B{12,1}=[33,46,47,B{12,1}];
B{14,1}=[67,68,B{14,1}];
B{16,1}=[88,89,90,91,106,B{16,1}];
B{17,1}=[107,108,109,B{17,1}];
B{13,2}=[58,B{13,2}];
B{1,3}=[41,87,B{1,3}];
B{4,3}=[64,B{4,3}];
B{5,3}=[87,B{5,3}];
B{13,3}=[59,B{13,3}];
B{3,10}=[60,B{3,10}];
% % M的添加仅为findequation.m提供帮助
B{1,9}=[20,22,86];
B{2,9}=[24,32,151];
B{3,9}=[13,17,20,34,56,105,137,138];
B{4,9}=[17,22,24,40,54,92];
B{5,9}=[14,57,86,101];
B{6,9}=19;
B{7,9}=[7,137];
B{8,9}=158;
B{10,9}=[16,148,149,158];
B{11,9}=40;
B{12,9}=100;
B{13,9}=42;
B{14,9}=151;
B{15,9}=[75,77];
B{18,9}=119;
B{19,9}=120;
B{20,9}=[132,138];
B{21,9}=128;
B{22,9}=133;
B{24,9}=148;
B{25,9}=149;
B{26,9}=145;

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

