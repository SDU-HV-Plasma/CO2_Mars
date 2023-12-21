function [nsp,species,K,number_of_reactions,mass,aux,auxv,Kv,number_of_reaction,number_of_vreaction]=rates_CO2m5z0(Tg,Vibration)

% R1:  e + CO2 => CO2^+ + 2e          k= 查表                                                ER=13.8eV
% R2:  e + CO2 => CO + O + e          k= 查表                   (激发解离)                   ER=11.46eV
% R3:  e + CO2 => CO + O^-            k= 查表                                            （附着反应无ER）
% R4:  e + O3  => O + O2 + e          k= 9.0e-16                (激发解离)           (应该存在ER，但未找到数据)
% R5:  e + O2  => 2O + e              k= 查表            (生成的O一般而言属于激发态)            ER=6eV
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
% R32: CO2 + e => O2^+ + C + 2e       k= K(1)*1/13                               ER未知,无法对应给出CO2e的电离阈值
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
% R58: e + CO => 2e + C^+ + O         k= 查表                                                ER=22eV
% R59: e + C => 2e + C^+              k= 查表                                               ER=11.2eV
% R60: C^+ + CO2 => CO^+ + CO         k= 1.1e-15
% R61: O2^+ + C => C^+ + O2           k= 5.2e-17
% R62: O2^- + CO2 + O2 => CO4^- + O2  k= 4.7e-41                                      (O2- + CO2 + M,M为O2)
% R63: CO4^- + CO2^+ => 2CO2 + O2     k= 5.0e-13
% R64: O2^+ + CO4^- => CO2 + 2O2      k= 3.0e-13
% R65: CO4^- + O => CO3^- + O2        k= 1.1e-16
% R66: CO4^- + O => CO2 + O2 + O^-    k= 1.4e-17
% R67: CO4^- + O => CO2 + O3^-        k= 1.4e-16
% R68: O^- + C => CO + e              k= 5.0e-16
% R69: e + CO4^+ => CO2 + O2          k= 1.61e-13*Te^(-0.5)
% R70: O2^+ + CO2 + CO2 => CO4^+ + CO2   k= 2.3e-41
% R71: e + O2^+ + CO2 => O2 + CO2     k=1.0e-38                                       (e + O2+ + M,M为CO2)
% R72: O2^- + CO2 + CO2 => CO4^- + CO2   k= 1.0e-41                         (O2- + CO2 + M.M为CO2)  ------1.339606
% R73: e + O2 + CO2 => O2^- + CO2        k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% R74: e + C2O2^+ => 2CO              k= 4.0e-13*Te^(-0.34) 
% R75: e + C2O3^+ => CO2 + CO         k= 5.4e-14*Te^(-0.7)
% R76: e + C2O4^+ => 2CO2             k= 2.0e-11*(Te^(-0.5))*(Tg*1.16e4)^(-1)
% R77: CO2^+ + CO2 + CO => C2O4^+ + CO  k= 3.0e-40
% R78: C2O3^+ + CO => CO2 + C2O2^+    k= 1.1e-15
% R79: C2O4^+ + CO => CO2 + C2O3^+    k= 9.0e-16
% R80: C2O3^+ + CO + CO2 => CO2 + C2O2^+ + CO2  k= 2.6e-38
% R81: C2O4^+ + CO + CO2 => CO2 + C2O3^+ + CO2  k= 4.2e-38
% R82: C2O2^+ + O2 => 2CO + O2^+      k= 5.0e-18
% R83: C2O2^+ + CO2 => CO^+ + CO + CO2    k= 1.0e-18
% R84: C2O2^+ + CO3^- => CO2 + 2CO + O       k= 5.0e-13
% R85: C2O2^+ + CO4^- => CO2 + 2CO + O2      k= 5.0e-13
% R86: C2O2^+ + O2^- => 2CO + O2      k= 6.0e-13
% R87: C2O3^+ + CO3^- => 2CO2 + CO + O       k= 5.0e-13
% R88: C2O3^+ + CO4^- => 2CO2 + CO + O2      k= 5.0e-13
% R89: C2O3^+ + O2^- => CO2 + CO + O2        k= 6.0e-13
% R90: C2O4^+ + CO2 => CO2^+ + CO2 + CO2  k= 1.0e-20
% R91: C2O4^+ + CO3^- => 3CO2 + O     k= 5.0e-13
% R92: C2O4^+ + CO4^- => 3CO2 + O2    k= 5.0e-13
% R93: C2O4^+ + O2^- => 2CO2 + O2     k= 6.0e-13
% R94: e + N2 => e + N2               k= 查表
% R95: e + N2 => 2e + N2^+            k= 查表                                          ER=15.58eV
% R96: e + N2 => e + 2N               k= 查表                                          ER=9.753eV
% R97: e + NO => N + O^-              k= 查表    
% R98: e + NO => 2e + NO^+            k= 查表                                          ER=9.26eV
% R99: e + N2O => N2 + O^-            k= 查表   
% R100:e + N2O => 2e + N2O^+          k= 查表   2.6×10-7  cm3                                       ER=12.89eV
% R101:e + NO2 => NO2^+ + 2e          k= 2.6e-15*Te^1/2*exp(-10/Te)                     ER=10eV
% R102:e + NO2 => NO^+ + O + 2e       k= 8.1e-15*Te^1/2*exp(-12.9/Te)                  ER=12.9eV
% R103:e + NO => N + O + e            k= 7.4e-15*exp(-6.5/Te)                           ER=6.5eV
% R104:e + NO2 => NO + O + e          k= 5.6e-15*exp(-3.11/Te)                         ER=3.11eV
% R105:e + N2O => N2 + O + e          k= 1.4e-15*exp(-1.67/Te)                         ER=1.67eV
% R106:e + NO2 + N2 => NO2^- + N2     k= 1.6e-43
% R107:e + NO2 + CO2 => NO2^- + CO2   k= 1.6e-43
% R108:O^- + N => NO + e              k= 2.2e-16
% R109:O^- + N2 => N2O + e            k= 1.0e-18
% R110:O^- + NO => NO2 + e            k= 2.5e-16
% R111:O2^- + N => NO2 + e            k= 4.0e-16
% R112:NO2^- + O => NO3 + e           k= 1.0e-18
% R113:NO2^- + N => N2 + O2 + e       k= 1.0e-18
% R114:NO3^- + O => NO2 + O2 + e      k= 1.0e-18
% R115:NO3^- + N => N2 + O3 + e       k= 1.0e-18
% R116:NO^- + N2 => NO + N2 + e       k= 5.0e-16
% R117:e + NO^+ => N + O              k= 7.0e-14*Te^(-0.5)
% R118:e + N2^+ => 2N                 k= 4.6e-14*Te^(-0.39)
% R119:e + NO2^+ => N + O2            k= 7.0e-14*Te^(-0.5)
% R120:e + N2O^+ => N2 + O            k= 7.0e-14*Te^(-0.5)
% R121:CO2^+ + NO => NO^+ + CO2       k= 1.2e-16
% R122:N2^+ + O2 => O2^+ + N2         k= 7.0e-17
% R123:N2^+ + NO => NO^+ + N2         k= 4.9e-16
% R124:N2^+ + CO => CO^+ + N2         k= 7.0e-17
% R125:N2^+ + CO2 => CO2^+ + N2       k= 9.0e-16
% R126:O2^+ + N => NO^+ + O           k= 1.9e-17
% R127:O2^+ + NO => NO^+ + O2         k= 7.0e-16
% R128:O2^+ + NO2 => NO2^+ + O2       k= 6.6e-16
% R129:NO^+ + O3 => NO2^+ + O2        k= 1.0e-20
% R130:NO2^+ + NO => NO^+ + NO2       k= 2.9e-16
% R131:N2O^+ + NO => NO^+ + N2O       k= 2.9e-16
% R132:O^- + NO2 => O2^- + NO         k= 1.0e-15
% R133:O^- + NO2 => NO2^- + O         k= 1.0e-15
% R134:O^- + N2O => NO^- + NO         k= 2.0e-16
% R135:O2^- + NO2 => NO2^- + O2       k= 2.0e-15
% R136:O2^- + N2O => O3^- + N2        k= 1.0e-17
% R137:O2^- + NO3 => NO3^- + O2       k= 5.0e-16
% R138:O3^- + NO => NO3^- + O         k= 1.0e-17
% R139:O3^- + NO => NO2^- + O2        k= 2.6e-18
% R140:O3^- + NO2 => NO2^- + O3       k= 7.0e-16
% R141:O3^- + NO2 => NO3^- + O2       k= 2.8e-16
% R142:O3^- + NO3 => NO3^- + O3       k= 5.0e-16
% R143:CO3^- + NO => NO2^- + CO2      k= 9.0e-18
% R144:CO3^- + NO2 => NO3^- + CO2     k= 1.0e-16
% R145:CO4^- + NO => NO3^- + CO2      k= 5.0e-17
% R146:NO^- + O2 => O2^- + NO         k= 9.0e-16
% R147:NO2^- + O3 => NO3^- + O2       k= 1.8e-17
% R148:NO2^- + NO2 => NO3^- + NO      k= 4.0e-18
% R149:NO2^- + NO3 => NO3^- + NO2     k= 3.0e-16
% R150:NO3^- + NO => NO2^- + NO2      k= 2.0e-17
% R151:O^- + NO + N2 => NO2^- + N2    k= 3.0e-40*(Tg*11600)^(-1)
% R152:CO2^+ + NO2^- => CO + NO2 + O  k= 6.0e-13
% R153:CO2^+ + NO3^- => CO + NO3 + O  k= 5.0e-13
% R154:NO^+ + CO3^- => CO2 + NO + O   k= 6.0e-13
% R155:NO^+ + CO4^- => CO2 + NO + O2  k= 6.0e-13
% R156:NO^+ + NO2^- => NO2 + N + O    k= 5.1e-13
% R157:NO^+ + NO3^- => NO3 + N + O    k= 8.1e-13
% R158:NO^+ + O2^- => O2 + N + O      k= 5.8e-13
% R159:NO^+ + O^- => NO + O           k= 4.9e-13
% R160:O2^+ + NO2^- => NO2 + 2O       k= 4.1e-13
% R161:O2^+ + NO3^- => NO3 + 2O       k= 1.3e-13
% R162:N2^+ + CO3^- => CO2 + N2 + O   k= 3.0e-13
% R163:N2^+ + CO4^- => CO2 + N2 + O2  k= 3.0e-13
% R164:N2^+ + NO2^- => NO2 + 2N       k= 4.1e-13
% R165:N2^+ + NO3^- => NO3 + 2N       k= 1.3e-13
% R166:N2^+ + O2^- => N2 + 2O         k= 4.2e-13
% R167:N2^+ + O^- => N2 + O           k= 1.0e-13
% R168:C2O4^+ + NO2^- => 2CO2 + NO2   k= 6.0e-13
% R169:C2O4^+ + NO3^- => 2CO2 + NO2 + O  k= 5.0e-13
% R170:C2O3^+ + NO2^- => CO2 + CO + NO2  k= 6.0e-13
% R171:C2O3^+ + NO3^- => CO2 + CO + NO2 + O  k=5.0e-13
% R172:C2O2^+ + NO2^- => 2CO + NO2    k= 6.0e-13
% R173:C2O2^+ + NO3^- => 2CO + NO2 + O   k=5.0e-13
% R174:N + CO2 => NO + CO             k= 1.0e-25
% R175:CO + NO2 => CO2 + NO           k= 1.48e-16*exp(-16967/(Tg*1.16e4))
% R176:NO + O + CO2 => CO2 + NO2      k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% R177:NO^- + CO => NO + CO + e       k= 5.0e-19
% R178:NO^- + CO2 => NO + CO2 + e     k= 8.3e-18
% 振动激发态部分
% Rv1: e + CO2 => e + CO2e1           k= 查表                                                EE=7eV  
% Rv2: e + CO2 => e + CO2e2           k= 查表                                              EE=10.5eV
% Rv3: e + CO2e1 => CO2^+ + 2e        k= K(1)*0.944                                    EI=6.8eV(未必准确)
% Rv4: e + CO2e2 => CO2^+ + 2e        k= K(1)*0.099                                    EI=3.3eV(未必准确)
% Rv5: e + CO2e1 => 2e + O + CO^+     k= K(44)*3.189                                   EI=12.5eV(未必准确)
% Rv6: e + CO2e2 => 2e + O + CO^+     k= K(44)*0.735                                    EI=9eV(未必准确)
% Rv7: CO2e1 + CO^+ => CO2^+ + CO     k= K(46)*0.916
% Rv8: CO2e2 + CO^+ => CO2^+ + CO     k= K(46)*0.093
% Rv9: e + CO2e1 => 2e + CO + O^+     k= K(53)*2.988                                   EI=12.1eV(未必准确)
% Rv10: e + CO2e2 => 2e + CO + O^+     k= K(53)*0.671                                  EI=8.6eV((未必准确))
% Rv11: O^+ + CO2e1 => O + CO2^+       k= K(56)*0.916
% Rv12: O^+ + CO2e2 => O + CO2^+       k= K(56)*0.093
% Rv13: C^+ + CO2e1 => CO^+ + CO       k= K(60)*3.189
% Rv14: C^+ + CO2e2 => CO^+ + CO       k= K(60)*0.735
% Rv15: e + CO2 => e + CO2va           k= 查表                                              EV=0.083
% Rv16: e + CO2 => e + CO2vb           k= 查表                                              EV=0.167
% Rv17: e + CO2 => e + CO2vc           k= 查表                                              EV=0.252
% Rv18: e + CO2 => e + CO2vd           k= 查表                                              EV=0.339
% Rv19: CO2va + CO2 => 2CO2            k= 7.14e-14*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv20: CO2va + CO => CO2 + CO         k= 5.00e-14*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv21: CO2va + O2 => CO2 + O2         k= 5.00e-14*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv22: CO2vb + CO2 => 2CO2            k= 1.071e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv23: CO2vb + CO => CO2 + CO         k= 3.213e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv24: CO2vb + O2 => CO2 + O2         k= 3.213e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv25: CO2vb + CO2 => CO2va + CO2     k= 1.438e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv26: CO2vb + CO => CO2va + CO       k= 1.007e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv27: CO2vb + O2 => CO2va + O2       k= 1.007e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv28: CO2vc + CO2 => CO2va + CO2     k= 1.071e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv29: CO2vc + CO => CO2va + CO       k= 3.213e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv30: CO2vc + O2 => CO2va + O2       k= 3.213e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv31: e + CO2 => e + CO2v1           k= 查表                                            EV=0.291eV
% Rv32: e + CO2 => e + CO2v2           k= 查表                                             EV=0.58eV
% Rv33: e + CO2 => e + CO2v3           k= 查表                                             EV=0.86eV
% Rv34: e + CO2 => e + CO2v4           k= 查表                                             EV=1.14eV
% Rv35: e + CO2v1 => CO2^+ + 2e        k= K(1)                                             ER=13.8eV
% Rv36: e + CO2v2 => CO2^+ + 2e        k= K(1)                                             ER=13.8eV
% Rv37: e + CO2v3 => CO2^+ + 2e        k= K(1)                                             ER=13.8eV
% Rv38: e + CO2v4 => CO2^+ + 2e        k= K(1)                                             ER=13.8eV
% Rv39: e + CO2v1 => CO + O + e        k= 查表                                            ER=11.169eV
% Rv40: e + CO2v2 => CO + O + e        k= 查表                                            ER=10.88eV
% Rv41: e + CO2v3 => CO + O + e        k= 查表                                             ER=10.6eV
% Rv42: e + CO2v4 => CO + O + e        k= 查表                                            ER=10.32eV
% Rv43: e + CO2v1 => CO + O-           k= 查表                                        附着过程无ER，下同
% Rv44: e + CO2v2 => CO + O-           k= 查表
% Rv45: e + CO2v3 => CO + O-           k= 查表                    
% Rv46: e + CO2v4 => CO + O-           k= 查表
% Rv47: e + CO2v1 => O2^+ + C + 2e     k= K(1)*1/13                                   ER未知，k值未修正
% Rv48: e + CO2v2 => O2^+ + C + 2e     k= K(1)*1/13                                   ER未知，k值未修正
% Rv49: e + CO2v3 => O2^+ + C + 2e     k= K(1)*1/13                                   ER未知，k值未修正
% Rv50: e + CO2v4 => O2^+ + C + 2e     k= K(1)*1/13                                   ER未知，k值未修正
% Rv51: e + CO2v1 => 2e + O + CO^+     k= 查表                                            ER=19.21eV
% Rv52: e + CO2v2 => 2e + O + CO^+     k= 查表                                            ER=18.92eV
% Rv53: e + CO2v3 => 2e + O + CO^+     k= 查表                                            ER=18.64eV
% Rv54: e + CO2v4 => 2e + O + CO^+     k= 查表                                            ER=18.36eV
% Rv55: e + CO2v1 => 2e + CO + O^+     k= 查表                                            ER=18.81eV
% Rv56: e + CO2v2 => 2e + CO + O^+     k= 查表                                            ER=18.52eV
% Rv57: e + CO2v3 => 2e + CO + O^+     k= 查表                                            ER=18.24eV
% Rv58: e + CO2v4 => 2e + CO + O^+     k= 查表                                            ER=17.96eV
% Rv59: CO2v1 + CO2 => CO2va + CO2     k= 4.25e-7*exp(-407*(Tg*11600)^(-1/3)+824*(Tg*11600)^(-2/3))
% Rv60: CO2v1 + CO => CO2va + CO       k= Kv(59)*0.3
% Rv61: CO2v1 + O2 => CO2va + O2       k= Kv(59)*0.4
% Rv62: CO2v1 + CO2 => CO2vb + CO2     k= 8.57e-7*exp(-404*(Tg*11600)^(-1/3)+1096*(Tg*11600)^(-2/3))
% Rv63: CO2v1 + CO => CO2vb + CO       k= Kv(62)*0.3
% Rv64: CO2v1 + O2 => CO2vb + O2       k= Kv(62)*0.4
% Rv65: CO2v1 + CO2 => CO2vc + CO2     k= 1.43e-11*exp(-252*(Tg*11600)^(-1/3)+685*(Tg*11600)^(-2/3))
% Rv66: CO2v1 + CO => CO2vc + CO       k= Kv(65)*0.3;
% Rv67: CO2v1 + O2 => CO2vc + O2       k= Kv(65)*0.4;
% Rv68: CO2v2 + CO2 => CO2v1 + CO2     k= 见下文
% Rv69: CO2v2 + CO => CO2v1 + CO       k= 见下文
% Rv70: CO2v2 + O2 => CO2v1 + O2       k= 见下文
% Rv71: CO2v3 + CO2 => CO2v2 + CO2     k= 见下文
% Rv72: CO2v3 + CO => CO2v2 + CO       k= 见下文
% Rv73: CO2v3 + O2 => CO2v2 + O2       k= 见下文
% Rv74: CO2v4 + CO2 => CO2v3 + CO2     k= 见下文
% Rv75: CO2v4 + CO => CO2v3 + CO       k= 见下文
% Rv76: CO2v4 + O2 => CO2v3 + O2       k= 见下文
% Rv77: CO2v1 + CO2 => CO2va + CO2vb   k= 1.06e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3))
% Rv78: CO2v2 + CO2 => CO2v1 + CO2va   k= 2.13e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3))
% Rv79: CO2v2 + CO2 => CO2v1 + CO2vb   k= 2.13e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3))
% Rv80: CO2v3 + CO2 => CO2v2 + CO2va   k= 2.13e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3))
% Rv81: CO2v3 + CO2 => CO2v2 + CO2vb   k= 2.13e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3))
% Rv82: CO2v4 + CO2 => CO2v3 + CO2va   k= 2.13e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3))
% Rv83: CO2v4 + CO2 => CO2v3 + CO2vb   k= 2.13e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3))
% Rv84: CO2vb + CO2 => CO2va + CO2va   k= 2.157e-15*exp(-88*(Tg*11600)^(-1/3)+233*(Tg*11600)^(-2/3))
% Rv85: CO2vc + CO2 => CO2va + CO2vb   k= 5.305e-15*exp(-88.2*(Tg*11600)^(-1/3)+233*(Tg*11600)^(-2/3))
% Rv86: CO2vc + CO2va => CO2vb + CO2vb k= 5.305e-15*exp(-88*(Tg*11600)^(-1/3)+233*(Tg*11600)^(-2/3))
% Rv87: CO2va + CO2 => CO + O + CO2    k= 3.91e-16*exp(-48688/(Tg*11600))
% Rv88: CO2vb + CO2 => CO + O + CO2    k= 3.91e-16*exp(-47852/(Tg*11600))
% Rv89: CO2vc + CO2 => CO + O + CO2    k= 3.91e-16*exp(-47110/(Tg*11600))
% Rv90: CO2vd + CO2 => CO + O + CO2    k= 3.91e-16*exp(-46368/(Tg*11600))
% Rv91: CO2v1 + CO2 => CO + O + CO2    k= 3.91e-16*exp(-46739/(Tg*11600))
% Rv92: CO2v2 + CO2 => CO + O + CO2    k= 3.91e-16*exp(-44048/(Tg*11600))
% Rv93: CO2v3 + CO2 => CO + O + CO2    k= 3.91e-16*exp(-41449/(Tg*11600))
% Rv94: CO2v4 + CO2 => CO + O + CO2    k= 3.91e-16*exp(-38851/(Tg*11600))
% Rv95: CO2va + O => CO + O2           k= 2.8e-17*exp(-26036/(Tg*11600))
% Rv96: CO2vb + O => CO + O2           k= 2.8e-17*exp(-25514/(Tg*11600))
% Rv97: CO2vc + O => CO + O2           k= 2.8e-17*exp(-25050/(Tg*11600))
% Rv98: CO2vd + O => CO + O2           k= 2.8e-17*exp(-24586/(Tg*11600))
% Rv99: CO2v1 + O => CO + O2           k= 2.8e-17*exp(-24818/(Tg*11600))
% Rv100:CO2v2 + O => CO + O2           k= 2.8e-17*exp(-23136/(Tg*11600))
% Rv101:CO2v3 + O => CO + O2           k= 2.8e-17*exp(-21512/(Tg*11600))
% Rv102:CO2v4 + O => CO + O2           k= 2.8e-17*exp(-19888/(Tg*11600))
% Rv103:e + CO2va => CO + O + e        k= 查表                                            ER=11.37eV
% Rv104:e + CO2vb => CO + O + e        k= 查表                                            ER=11.293eV
% Rv105:e + CO2vc => CO + O + e        k= 查表                                            ER=11.208eV
% Rv106:e + CO2vd => CO + O + e        k= 查表                                            ER=11.121eV
% Rv107:e + CO2va => CO + O-           k= 查表
% Rv108:e + CO2vb => CO + O-           k= 查表
% Rv109:e + CO2vc => CO + O-           k= 查表
% Rv110:e + CO2vd => CO + O-           k= 查表
% Rv111:e + CO2va => CO2^+ + 2e        k= K(1)                                             ER=13.8eV
% Rv112:e + CO2vb => CO2^+ + 2e        k= K(1)                                             ER=13.8eV
% Rv113:e + CO2vc => CO2^+ + 2e        k= K(1)                                             ER=13.8eV
% Rv114:e + CO2vd => CO2^+ + 2e        k= K(1)                                             ER=13.8eV
% Rv115:e + CO2va => O2^+ + C + 2e     k= K(1)*1/13                                    ER未知，k值未修正
% Rv116:e + CO2vb => O2^+ + C + 2e     k= K(1)*1/13                                    ER未知，k值未修正
% Rv117:e + CO2vc => O2^+ + C + 2e     k= K(1)*1/13                                    ER未知，k值未修正
% Rv118:e + CO2vd => O2^+ + C + 2e     k= K(1)*1/13                                    ER未知，k值未修正
% Rv119:e + CO2va => 2e + O + CO^+     k= 查表                                            ER=19.417eV                        
% Rv120:e + CO2vb => 2e + O + CO^+     k= 查表                                            ER=19.333eV                   
% Rv121:e + CO2vc => 2e + O + CO^+     k= 查表                                            ER=19.248eV                                 
% Rv122:e + CO2vd => 2e + O + CO^+     k= 查表                                            ER=19.161eV
% Rv123:e + CO2va => 2e + CO + O^+     k= 查表                                            ER=19.017eV
% Rv124:e + CO2vb => 2e + CO + O^+     k= 查表                                            ER=18.933eV                                   
% Rv125:e + CO2vc => 2e + CO + O^+     k= 查表                                            ER=18.848eV                                 
% Rv126:e + CO2vd => 2e + CO + O^+     k= 查表                                            ER=18.761eV                                       
% Rv127:CO2e1 + C => 2CO               k= 1.0e-21
% Rv128:CO2e2 + C => 2CO               k= 1.0e-21
% Rv129:CO2va + C => 2CO               k= 1.0e-21
% Rv130:CO2vb + C => 2CO               k= 1.0e-21
% Rv131:CO2vc + C => 2CO               k= 1.0e-21
% Rv132:CO2vd + C => 2CO               k= 1.0e-21
% Rv133:CO2v1 + C => 2CO               k= 1.0e-21
% Rv134:CO2v2 + C => 2CO               k= 1.0e-21
% Rv135:CO2v3 + C => 2CO               k= 1.0e-21
% Rv136:CO2v4 + C => 2CO               k= 1.0e-21
% Rv137:O^- + CO2e1 + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv138:O^- + CO2e2 + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv139:O^- + CO2va + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv140:O^- + CO2vb + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv141:O^- + CO2vc + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv142:O^- + CO2vd + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv143:O^- + CO2v1 + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv144:O^- + CO2v2 + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv145:O^- + CO2v3 + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv146:O^- + CO2v4 + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv147:CO2e1 + O^- + CO => CO3^- + CO   k= 1.5e-40
% Rv148:CO2e2 + O^- + CO => CO3^- + CO   k= 1.5e-40
% Rv149:CO2va + O^- + CO => CO3^- + CO   k= 1.5e-40
% Rv150:CO2vb + O^- + CO => CO3^- + CO   k= 1.5e-40
% Rv151:CO2vc + O^- + CO => CO3^- + CO   k= 1.5e-40
% Rv152:CO2vd + O^- + CO => CO3^- + CO   k= 1.5e-40
% Rv153:CO2v1 + O^- + CO => CO3^- + CO   k= 1.5e-40
% Rv154:CO2v2 + O^- + CO => CO3^- + CO   k= 1.5e-40
% Rv155:CO2v3 + O^- + CO => CO3^- + CO   k= 1.5e-40
% Rv156:CO2v4 + O^- + CO => CO3^- + CO   k= 1.5e-40
% Rv157:CO2e1 + O^- + O2 => CO3^- + O2   k= 3.1e-40
% Rv158:CO2e2 + O^- + O2 => CO3^- + O2   k= 3.1e-40
% Rv159:CO2va + O^- + O2 => CO3^- + O2   k= 3.1e-40
% Rv160:CO2vb + O^- + O2 => CO3^- + O2   k= 3.1e-40
% Rv161:CO2vc + O^- + O2 => CO3^- + O2   k= 3.1e-40
% Rv162:CO2vd + O^- + O2 => CO3^- + O2   k= 3.1e-40
% Rv163:CO2v1 + O^- + O2 => CO3^- + O2   k= 3.1e-40
% Rv164:CO2v2 + O^- + O2 => CO3^- + O2   k= 3.1e-40
% Rv165:CO2v3 + O^- + O2 => CO3^- + O2   k= 3.1e-40
% Rv166:CO2v4 + O^- + O2 => CO3^- + O2   k= 3.1e-40
% Rv167:CO2va + CO^+ => CO2^+ + CO       k= 1.0e-15
% Rv168:CO2vb + CO^+ => CO2^+ + CO       k= 1.0e-15
% Rv169:CO2vc + CO^+ => CO2^+ + CO       k= 1.0e-15
% Rv170:CO2vd + CO^+ => CO2^+ + CO       k= 1.0e-15
% Rv171:CO2v1 + CO^+ => CO2^+ + CO       k= 1.0e-15
% Rv172:CO2v2 + CO^+ => CO2^+ + CO       k= 1.0e-15
% Rv173:CO2v3 + CO^+ => CO2^+ + CO       k= 1.0e-15
% Rv174:CO2v4 + CO^+ => CO2^+ + CO       k= 1.0e-15
% Rv175:O3^- + CO2e1 => CO3^- + O2       k= 5.5e-16
% Rv176:O3^- + CO2e2 => CO3^- + O2       k= 5.5e-16
% Rv177:O3^- + CO2va => CO3^- + O2       k= 5.5e-16
% Rv178:O3^- + CO2vb => CO3^- + O2       k= 5.5e-16
% Rv179:O3^- + CO2vc => CO3^- + O2       k= 5.5e-16
% Rv180:O3^- + CO2vd => CO3^- + O2       k= 5.5e-16
% Rv181:O3^- + CO2v1 => CO3^- + O2       k= 5.5e-16
% Rv182:O3^- + CO2v2 => CO3^- + O2       k= 5.5e-16
% Rv183:O3^- + CO2v3 => CO3^- + O2       k= 5.5e-16
% Rv184:O3^- + CO2v4 => CO3^- + O2       k= 5.5e-16
% Rv185:O^+ + CO2e1 => O2^+ + CO         k= 9.4e-16
% Rv186:O^+ + CO2e2 => O2^+ + CO         k= 9.4e-16
% Rv187:O^+ + CO2va => O2^+ + CO         k= 9.4e-16
% Rv188:O^+ + CO2vb => O2^+ + CO         k= 9.4e-16
% Rv189:O^+ + CO2vc => O2^+ + CO         k= 9.4e-16
% Rv190:O^+ + CO2vd => O2^+ + CO         k= 9.4e-16
% Rv191:O^+ + CO2v1 => O2^+ + CO         k= 9.4e-16
% Rv192:O^+ + CO2v2 => O2^+ + CO         k= 9.4e-16
% Rv193:O^+ + CO2v3 => O2^+ + CO         k= 9.4e-16
% Rv194:O^+ + CO2v4 => O2^+ + CO         k= 9.4e-16
% Rv195:O^+ + CO2va => O + CO2^+         k= 4.5e-16
% Rv196:O^+ + CO2vb => O + CO2^+         k= 4.5e-16
% Rv197:O^+ + CO2vc => O + CO2^+         k= 4.5e-16
% Rv198:O^+ + CO2vd => O + CO2^+         k= 4.5e-16
% Rv199:O^+ + CO2v1 => O + CO2^+         k= 4.5e-16
% Rv200:O^+ + CO2v2 => O + CO2^+         k= 4.5e-16
% Rv201:O^+ + CO2v3 => O + CO2^+         k= 4.5e-16
% Rv202:O^+ + CO2v4 => O + CO2^+         k= 4.5e-16
% Rv203:C^+ + CO2va => CO^+ + CO         k= 1.1e-15
% Rv204:C^+ + CO2vb => CO^+ + CO         k= 1.1e-15
% Rv205:C^+ + CO2vc => CO^+ + CO         k= 1.1e-15
% Rv206:C^+ + CO2vd => CO^+ + CO         k= 1.1e-15
% Rv207:C^+ + CO2v1 => CO^+ + CO         k= 1.1e-15
% Rv208:C^+ + CO2v2 => CO^+ + CO         k= 1.1e-15
% Rv209:C^+ + CO2v3 => CO^+ + CO         k= 1.1e-15
% Rv210:C^+ + CO2v4 => CO^+ + CO         k= 1.1e-15
% Rv211:O2^- + CO2e1 + O2 => CO4^- + O2  k= 4.7e-41
% Rv212:O2^- + CO2e2 + O2 => CO4^- + O2  k= 4.7e-41
% Rv213:O2^- + CO2va + O2 => CO4^- + O2  k= 4.7e-41
% Rv214:O2^- + CO2vb + O2 => CO4^- + O2  k= 4.7e-41
% Rv215:O2^- + CO2vc + O2 => CO4^- + O2  k= 4.7e-41
% Rv216:O2^- + CO2vd + O2 => CO4^- + O2  k= 4.7e-41
% Rv217:O2^- + CO2v1 + O2 => CO4^- + O2  k= 4.7e-41
% Rv218:O2^- + CO2v2 + O2 => CO4^- + O2  k= 4.7e-41
% Rv219:O2^- + CO2v3 + O2 => CO4^- + O2  k= 4.7e-41
% Rv220:O2^- + CO2v4 + O2 => CO4^- + O2  k= 4.7e-41
% Rv221:e + CO2e1 => CO + O + e          k= 查表
% Rv222:e + CO2e2 => CO + O + e          k= 查表
% Rv223:e + CO2e1 => CO + O^-            k= 查表
% Rv224:e + CO2e2 => CO + O^-            k= 查表
% Rv225:O + O2 + CO2e1 => O3 + CO2e1     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv226:O + O2 + CO2e2 => O3 + CO2e2     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv227:O + O2 + CO2va => O3 + CO2va     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv228:O + O2 + CO2vb => O3 + CO2vb     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv229:O + O2 + CO2vc => O3 + CO2vc     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv230:O + O2 + CO2vd => O3 + CO2vd     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv231:O + O2 + CO2v1 => O3 + CO2v1     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv232:O + O2 + CO2v2 => O3 + CO2v2     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv233:O + O2 + CO2v3 => O3 + CO2v3     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv234:O + O2 + CO2v4 => O3 + CO2v4     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv235:O^- + CO2e1 => O + CO2e1 + e    k= 4.0e-18
% Rv236:O^- + CO2e2 => O + CO2e2 + e    k= 4.0e-18
% Rv237:O^- + CO2va => O + CO2va + e    k= 4.0e-18
% Rv238:O^- + CO2vb => O + CO2vb + e    k= 4.0e-18
% Rv239:O^- + CO2vc => O + CO2vc + e    k= 4.0e-18
% Rv240:O^- + CO2vd => O + CO2vd + e    k= 4.0e-18
% Rv241:O^- + CO2v1 => O + CO2v1 + e    k= 4.0e-18
% Rv242:O^- + CO2v2 => O + CO2v2 + e    k= 4.0e-18
% Rv243:O^- + CO2v3 => O + CO2v3 + e    k= 4.0e-18
% Rv244:O^- + CO2v4 => O + CO2v4 + e    k= 4.0e-18
% Rv245:O^- + CO2 + CO2e1 => CO3^- + CO2e1 k= 9.0e-41
% Rv246:O^- + CO2 + CO2e2 => CO3^- + CO2e2 k= 9.0e-41
% Rv247:O^- + CO2 + CO2va => CO3^- + CO2va k= 9.0e-41
% Rv248:O^- + CO2 + CO2vb => CO3^- + CO2vb k= 9.0e-41
% Rv249:O^- + CO2 + CO2vc => CO3^- + CO2vc k= 9.0e-41
% Rv250:O^- + CO2 + CO2vd => CO3^- + CO2vd k= 9.0e-41
% Rv251:O^- + CO2 + CO2v1 => CO3^- + CO2v1 k= 9.0e-41
% Rv252:O^- + CO2 + CO2v2 => CO3^- + CO2v2 k= 9.0e-41
% Rv253:O^- + CO2 + CO2v3 => CO3^- + CO2v3 k= 9.0e-41
% Rv254:O^- + CO2 + CO2v4 => CO3^- + CO2v4 k= 9.0e-41
% Rv255:O + O + CO2e1 => O2 + CO2e1     k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv256:O + O + CO2e2 => O2 + CO2e2     k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv257:O + O + CO2va => O2 + CO2va     k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv258:O + O + CO2vb => O2 + CO2vb     k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv259:O + O + CO2vc => O2 + CO2vc     k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv260:O + O + CO2vd => O2 + CO2vd     k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv261:O + O + CO2v1 => O2 + CO2v1     k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv262:O + O + CO2v2 => O2 + CO2v2     k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv263:O + O + CO2v3 => O2 + CO2v3     k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv264:O + O + CO2v4 => O2 + CO2v4     k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv265:O2^+ + CO2e1 + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv266:O2^+ + CO2e2 + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv267:O2^+ + CO2va + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv268:O2^+ + CO2vb + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv269:O2^+ + CO2vc + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv270:O2^+ + CO2vd + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv271:O2^+ + CO2v1 + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv272:O2^+ + CO2v2 + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv273:O2^+ + CO2v3 + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv274:O2^+ + CO2v4 + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv275:e + O2^+ + CO2e1 => O2 + CO2e1     k=1.0e-38
% Rv276:e + O2^+ + CO2e2 => O2 + CO2e2     k=1.0e-38
% Rv277:e + O2^+ + CO2va => O2 + CO2va     k=1.0e-38
% Rv278:e + O2^+ + CO2vb => O2 + CO2vb     k=1.0e-38
% Rv279:e + O2^+ + CO2vc => O2 + CO2vc     k=1.0e-38
% Rv280:e + O2^+ + CO2vd => O2 + CO2vd     k=1.0e-38
% Rv281:e + O2^+ + CO2v1 => O2 + CO2v1     k=1.0e-38
% Rv282:e + O2^+ + CO2v2 => O2 + CO2v2     k=1.0e-38
% Rv283:e + O2^+ + CO2v3 => O2 + CO2v3     k=1.0e-38
% Rv284:e + O2^+ + CO2v4 => O2 + CO2v4     k=1.0e-38
% Rv285:O2^- + CO2e1 + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv286:O2^- + CO2e2 + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv287:O2^- + CO2va + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv288:O2^- + CO2vb + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv289:O2^- + CO2vc + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv290:O2^- + CO2vd + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv291:O2^- + CO2v1 + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv292:O2^- + CO2v2 + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv293:O2^- + CO2v3 + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv294:O2^- + CO2v4 + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv295:O2^- + CO2 + CO2e1 => CO4^- + CO2e1   k= 1.0e-41
% Rv296:O2^- + CO2 + CO2e2 => CO4^- + CO2e2   k= 1.0e-41
% Rv297:O2^- + CO2 + CO2va => CO4^- + CO2va   k= 1.0e-41
% Rv298:O2^- + CO2 + CO2vb => CO4^- + CO2vb   k= 1.0e-41
% Rv299:O2^- + CO2 + CO2vc => CO4^- + CO2vc   k= 1.0e-41
% Rv300:O2^- + CO2 + CO2vd => CO4^- + CO2vd   k= 1.0e-41
% Rv301:O2^- + CO2 + CO2v1 => CO4^- + CO2v1   k= 1.0e-41
% Rv302:O2^- + CO2 + CO2v2 => CO4^- + CO2v2   k= 1.0e-41
% Rv303:O2^- + CO2 + CO2v3 => CO4^- + CO2v3   k= 1.0e-41
% Rv304:O2^- + CO2 + CO2v4 => CO4^- + CO2v4   k= 1.0e-41
% Rv305:CO2vc + CO2 => CO2vb + CO2   k= 2.897e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv306:CO2vc + CO => CO2vb + CO     k= 2.028e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv307:CO2vc + O2 => CO2vb + O2     k= 2.028e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv308:CO2vd + CO2 => CO2vb + CO2   k= 1.582e-11*exp(-272*(Tg*1.16e4)^(-1/3)+437*(Tg*1.16e4)^(-2/3))
% Rv309:CO2vd + CO => CO2vb + CO     k= 4.745e-11*exp(-272*(Tg*1.16e4)^(-1/3)+437*(Tg*1.16e4)^(-2/3))
% Rv310:CO2vd + O2 => CO2vb + O2     k= 4.745e-11*exp(-272*(Tg*1.16e4)^(-1/3)+437*(Tg*1.16e4)^(-2/3))
% Rv311:CO2vd + CO2 => CO2vc + CO2   k= 4.321e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv312:CO2vd + CO => CO2vc + CO     k= 3.025e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv313:CO2vd + O2 => CO2vc + O2     k= 3.025e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv314:CO2vd + CO2 => CO2v1 + CO2   k= 1.775e-17*exp(-108*(Tg*1.16e4)^(-1/3)+165*(Tg*1.16e4)^(-2/3))
% Rv315:CO2vd + CO => CO2v1 + CO     k= 2.13e-17*exp(-108*(Tg*1.16e4)^(-1/3)+165*(Tg*1.16e4)^(-2/3))
% Rv316:CO2vd + O2 => CO2v1 + O2     k= 2.13e-17*exp(-108*(Tg*1.16e4)^(-1/3)+165*(Tg*1.16e4)^(-2/3))
% Rv317:CO2vd + CO2 => CO2vb + CO2vb k= 2.384e-15*exp(-89*(Tg*1.16e4)^(-1/3)+234*(Tg*1.16e4)^(-2/3))
% Rv318:CO2vd + CO2 => CO2va + CO2vc k= 6.48e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3))
% Rv319:CO2vd + CO2va => CO2vb + CO2vc   k= 1.442e-14*exp(-88.2*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3))
% Rv320:CO2vd + CO2vb => CO2vc + CO2vc   k= 2.628e-14*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3))
% Rv321:CO2v1 + CO2v1 => CO2 + CO2v2 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv322:CO2v1 + CO2v2 => CO2 + CO2v3 k= 2.927e-17*exp(21.4*(Tg*1.16e4)^(-1/3)-53*(Tg*1.16e4)^(-2/3))
% Rv323:CO2v1 + CO2v3 => CO2 + CO2v4 k= 4.825e-17*exp(20.1*(Tg*1.16e4)^(-1/3)-65*(Tg*1.16e4)^(-2/3))
% Rv324:CO2v2 + CO2v2 => CO2v1 + CO2v3   k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv325:CO2v2 + CO2v3 => CO2v1 + CO2v4   k= 2.927e-17*exp(21.4*(Tg*1.16e4)^(-1/3)-53*(Tg*1.16e4)^(-2/3))
% Rv326:CO2v3 + CO2v3 => CO2v2 + CO2v4   k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv327:e + O2^+ + CO2v1a => O2 + CO2v1a k=1.0e-38
% Rv328:e + O2^+ + CO2v1b => O2 + CO2v1b k=1.0e-38
% Rv329:e + O2^+ + CO2v1c => O2 + CO2v1c k=1.0e-38
% Rv330:CO2v1 + CO2v4 => CO2v2 + CO2v3   k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv331:CO2v1 + CO2v3 => CO2v2 + CO2v2   k= Kv(330)
% Rv332:CO2v2 + CO2v4 => CO2v3 + CO2v3   k= Kv(330)
% Rv333:e + CO2v1 => e + CO2v2           k= 查表                                                EV=0.29eV
% Rv334:e + CO2v2 => e + CO2v3           k= 查表                                                EV=0.28eV
% Rv335:e + CO2v3 => e + CO2v4           k= 查表                                                EV=0.28eV
% Rv336:e + CO2 => e + CO2v1a            k= 查表                                                EV=0.339eV
% Rv337:e + CO2 => e + CO2v1b            k= 查表                                                EV=0.422eV
% Rv338:e + CO2 => e + CO2v1c            k= 查表                                                EV=0.505eV
% Rv339:e + CO2v1a => CO2^+ + 2e         k= K(1)                                                ER=13.8eV
% Rv340:e + CO2v1b => CO2^+ + 2e         k= K(1)                                                ER=13.8eV
% Rv341:e + CO2v1c => CO2^+ + 2e         k= K(1)                                                ER=13.8eV
% Rv342:e + CO2v1a => O2^+ + C + 2e      k= ?
% Rv343:e + CO2v1b => O2^+ + C + 2e      k= ?
% Rv344:e + CO2v1c => O2^+ + C + 2e      k= ?
% Rv345:e + CO2v1a => CO^+ + O + 2e      k= 查表                                                ER=19.161eV
% Rv346:e + CO2v1b => CO^+ + O + 2e      k= 查表                                                ER=19.078eV
% Rv347:e + CO2v1c => CO^+ + O + 2e      k= 查表                                                ER=18.995eV
% Rv348:e + CO2v1a => O^+ + CO + 2e      k= 查表                                                ER=18.761eV
% Rv349:e + CO2v1b => O^+ + CO + 2e      k= 查表                                                ER=18.678eV
% Rv350:e + CO2v1c => O^+ + CO + 2e      k= 查表                                                ER=18.595eV
% Rv351:e + CO2v1a => CO + O + e         k= 查表                                                ER=11.121eV
% Rv352:e + CO2v1b => CO + O + e         k= 查表                                                ER=11.038eV
% Rv353:e + CO2v1c => CO + O + e         k= 查表                                                ER=10.955eV
% Rv354:e + CO2v1a => CO + O^-           k= 查表
% Rv355:e + CO2v1b => CO + O^-           k= 查表
% Rv356:e + CO2v1c => CO + O^-           k= 查表
% Rv357:CO2v1a + CO2 => CO2vc + CO2      k= 8.568e-7*exp(-404*(Tg*1.16e4)^(-1/3)+1096*(Tg*1.16e4)^(-2/3))
% Rv358:CO2v1a + CO => CO2vc + CO        k= 2.57e-7*exp(-404*(Tg*1.16e4)^(-1/3)+1096*(Tg*1.16e4)^(-2/3))
% Rv359:CO2v1a + O2 => CO2vc + O2        k= 3.427e-7*exp(-404*(Tg*1.16e4)^(-1/3)+1096*(Tg*1.16e4)^(-2/3))
% Rv360:CO2v1a + CO2 => CO2vd + CO2      k= 1.431e-11*exp(-252*(Tg*1.16e4)^(-1/3)+685*(Tg*1.16e4)^(-2/3))
% Rv361:CO2v1a + CO => CO2vd + CO        k= 4.293e-12*exp(-252*(Tg*1.16e4)^(-1/3)+685*(Tg*1.16e4)^(-2/3))
% Rv362:CO2v1a + O2 => CO2vd + O2        k= 5.724e-12*exp(-252*(Tg*1.16e4)^(-1/3)+685*(Tg*1.16e4)^(-2/3))
% Rv363:CO2v1a + CO2 => CO2v1 + CO2      k= 7.143e-14*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv364:CO2v1a + CO => CO2v1 + CO        k= 5e-14*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv365:CO2v1a + O2 => CO2v1 + O2        k= 5e-14*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv366:CO2v1b + CO2 => CO2vc + CO2      k= 3.218e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv367:CO2v1b + CO => CO2vc + CO        k= 9.976e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv368:CO2v1b + O2 => CO2vc + O2        k= 9.976e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv369:CO2v1b + CO2 => CO2vd + CO2      k= 6.447e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv370:CO2v1b + CO => CO2vd + CO        k= 4.513e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv371:CO2v1b + O2 => CO2vd + O2        k= 4.513e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv372:CO2v1b + CO2 => CO2v1 + CO2      k= 1.071e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv373:CO2v1b + CO => CO2v1 + CO        k= 3.32e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv374:CO2v1b + O2 => CO2v1 + O2        k= 3.32e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv375:CO2v1b + CO2 => CO2v1a + CO2     k= 6.447e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv376:CO2v1b + CO => CO2v1a + CO       k= 4.513e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv377:CO2v1b + O2 => CO2v1a + O2       k= 4.513e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv378:CO2v1c + CO2 => CO2v1a + CO2     k= 1.071e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv379:CO2v1c + CO => CO2v1a + CO       k= 3.32e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv380:CO2v1c + O2 => CO2v1a + O2       k= 3.32e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv381:CO2v1c + CO2 => CO2v1b + CO2     k= 2.897e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv382:CO2v1c + CO => CO2v1b + CO       k= 2.028e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv383:CO2v1c + O2 => CO2v1b + O2       k= 2.028e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv384:CO2v2 + CO2 => CO2v1a + CO2      k= 8.568e-7*exp(-406*(Tg*1.16e4)^(-1/3)+829*(Tg*1.16e4)^(-2/3))
% Rv385:CO2v2 + CO => CO2v1a + CO        k= 2.57e-7*exp(-406*(Tg*1.16e4)^(-1/3)+829*(Tg*1.16e4)^(-2/3))
% Rv386:CO2v2 + O2 => CO2v1a + O2        k= 3.427e-7*exp(-406*(Tg*1.16e4)^(-1/3)+829*(Tg*1.16e4)^(-2/3))
% Rv387:CO2v2 + CO2 => CO2v1b + CO2      k= 1.725e-6*exp(-404*(Tg*1.16e4)^(-1/3)+1098*(Tg*1.16e4)^(-2/3))
% Rv388:CO2v2 + CO => CO2v1b + CO        k= 5.175e-7*exp(-404*(Tg*1.16e4)^(-1/3)+1098*(Tg*1.16e4)^(-2/3))
% Rv389:CO2v2 + O2 => CO2v1b + O2        k= 6.9e-7*exp(-404*(Tg*1.16e4)^(-1/3)+1098*(Tg*1.16e4)^(-2/3))
% Rv390:CO2v2 + CO2 => CO2v1c + CO2      k= 2.882e-11*exp(-253*(Tg*1.16e4)^(-1/3)+683*(Tg*1.16e4)^(-2/3))
% Rv391:CO2v2 + CO => CO2v1c + CO        k= 8.646e-12*exp(-253*(Tg*1.16e4)^(-1/3)+683*(Tg*1.16e4)^(-2/3))
% Rv392:CO2v2 + O2 => CO2v1c + O2        k= 1.153e-11*exp(-253*(Tg*1.16e4)^(-1/3)+683*(Tg*1.16e4)^(-2/3))
% Rv393:CO2v2 + CO2 => CO2v1a + CO2vb    k= 4.299e-11*exp(-241*(Tg*1.16e4)^(-1/3)+637*(Tg*1.16e4)^(-2/3))
% Rv394:CO2v2 + CO2 => CO2v1b + CO2va    k= 4.299e-11*exp(-241*(Tg*1.16e4)^(-1/3)+635*(Tg*1.16e4)^(-2/3))
% Rv395:CO2v1a + CO2 => CO2v1 + CO2va    k= 1.071e-15*exp(-87.7*(Tg*1.16e4)^(-1/3)+230*(Tg*1.16e4)^(-2/3))
% Rv396:CO2v1a + CO2va => CO2v1 + CO2vb  k= 2.157e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3))
% Rv397:CO2v1a + CO2vb => CO2v1 + CO2vc  k= 4.344e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3))
% Rv398:CO2v1a + CO2vc => CO2v1 + CO2vd  k= 6.48e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3))
% Rv399:CO2v1b + CO2 => CO2va + CO2vd    k= 9.667e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3))
% Rv400:CO2v1b + CO2va => CO2vb + CO2vd  k= 9.667e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3))
% Rv401:CO2v1b + CO2vb => CO2vc + CO2vd  k= 4.332e-14*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3))
% Rv402:CO2v1b + CO2vc => CO2vd + CO2vd  k= 5.848e-14*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3))
% Rv403:CO2v1b + CO2v1 => CO2v1a + CO2v1a  k= 2.157e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3))
% Rv404:CO2v1c + CO2 => CO2va + CO2vd    k= 2.378e-14*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3))
% Rv405:CO2v1c + CO2va => CO2vb + CO2vd  k= 5.292e-14*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3))
% Rv406:CO2v1c + CO2vb => CO2vc + CO2vd  k= 1.066e-13*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3))
% Rv407:CO2v1c + CO2vc => CO2vd + CO2vd  k= 1.438e-13*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3))
% Rv408:CO2v1c + CO2v1 => CO2v1a + CO2v1b  k= 5.305e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3))
% Rv409:CO2v1c + CO2v1a => CO2v1b + CO2v1b k= 5.305e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3))
% Rv410:CO2 + CO2va => CO + O + CO2va    k= 3.91e-16*exp(-48688/(Tg*11600))
% Rv411:CO2 + CO2vb => CO + O + CO2vb    k= 3.91e-16*exp(-47852/(Tg*11600))
% Rv412:CO2 + CO2vc => CO + O + CO2vc    k= 3.91e-16*exp(-47110/(Tg*11600))
% Rv413:CO2 + CO2vd => CO + O + CO2vd    k= 3.91e-16*exp(-46368/(Tg*11600))
% Rv414:CO2 + CO2v1 => CO + O + CO2v1    k= 3.91e-16*exp(-46739/(Tg*11600))
% Rv415:CO2 + CO2v2 => CO + O + CO2v2    k= 3.91e-16*exp(-44048/(Tg*11600))
% Rv416:CO2 + CO2v3 => CO + O + CO2v3    k= 3.91e-16*exp(-41449/(Tg*11600))
% Rv417:CO2 + CO2v4 => CO + O + CO2v4    k= 3.91e-16*exp(-38851/(Tg*11600))
% Rv418:CO2v1a + CO2 => CO + O + CO2     k= 3.91e-16*exp(-46284/(Tg*11600))
% Rv419:CO2v1b + CO2 => CO + O + CO2     k= 3.91e-16*exp(-45514/(Tg*11600))
% Rv420:CO2v1c + CO2 => CO + O + CO2     k= 3.91e-16*exp(-44744/(Tg*11600))
% Rv421:CO2 + CO2v1a => CO + O + CO2v1a  k= 3.91e-16*exp(-46284/(Tg*11600))
% Rv422:CO2 + CO2v1b => CO + O + CO2v1b  k= 3.91e-16*exp(-45514/(Tg*11600))
% Rv423:CO2 + CO2v1c => CO + O + CO2v1c  k= 3.91e-16*exp(-44744/(Tg*11600))
% Rv424:CO2v1a + O => CO + O2            k= 2.8e-17*exp(-24534/(Tg*11600))
% Rv425:CO2v1b + O => CO + O2            k= 2.8e-17*exp(-24052/(Tg*11600))
% Rv426:CO2v1c + O => CO + O2            k= 2.8e-17*exp(-23571/(Tg*11600))
% Rv427:CO2v1a + C => 2CO                k= 1.0e-21
% Rv428:CO2v1b + C => 2CO                k= 1.0e-21
% Rv429:CO2v1c + C => 2CO                k= 1.0e-21
% Rv430:O^- + CO2v1a + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv431:O^- + CO2v1b + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv432:O^- + CO2v1c + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv433:O^- + CO2v1a + CO => CO3^- + CO   k= 1.5e-40
% Rv434:O^- + CO2v1b + CO => CO3^- + CO   k= 1.5e-40
% Rv435:O^- + CO2v1c + CO => CO3^- + CO   k= 1.5e-40
% Rv436:O^- + CO2v1a + O2 => CO3^- + O2   k= 3.1e-40
% Rv437:O^- + CO2v1b + O2 => CO3^- + O2   k= 3.1e-40
% Rv438:O^- + CO2v1c + O2 => CO3^- + O2   k= 3.1e-40
% Rv439:O^- + CO2 + CO2v1a => CO3^- + CO2v1a k= 9.0e-41
% Rv440:O^- + CO2 + CO2v1b => CO3^- + CO2v1b k= 9.0e-41
% Rv441:O^- + CO2 + CO2v1c => CO3^- + CO2v1c k= 9.0e-41
% Rv442:CO2v1a + CO^+ => CO2^+ + CO       k= 1.0e-15
% Rv443:CO2v1b + CO^+ => CO2^+ + CO       k= 1.0e-15
% Rv444:CO2v1c + CO^+ => CO2^+ + CO       k= 1.0e-15
% Rv445:O3^- + CO2v1a => CO3^- + O2       k= 5.5e-16
% Rv446:O3^- + CO2v1b => CO3^- + O2       k= 5.5e-16
% Rv447:O3^- + CO2v1c => CO3^- + O2       k= 5.5e-16
% Rv448:O^+ + CO2v1a => O2^+ + CO         k= 9.4e-16
% Rv449:O^+ + CO2v1b => O2^+ + CO         k= 9.4e-16
% Rv450:O^+ + CO2v1c => O2^+ + CO         k= 9.4e-16
% Rv451:O^+ + CO2v1a => O + CO2^+         k= 4.5e-16
% Rv452:O^+ + CO2v1b => O + CO2^+         k= 4.5e-16
% Rv453:O^+ + CO2v1c => O + CO2^+         k= 4.5e-16
% Rv454:C^+ + CO2v1a => CO^+ + CO         k= 1.1e-15
% Rv455:C^+ + CO2v1b => CO^+ + CO         k= 1.1e-15
% Rv456:C^+ + CO2v1c => CO^+ + CO         k= 1.1e-15
% Rv457:O2^+ + CO2v1a + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv458:O2^+ + CO2v1b + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv459:O2^+ + CO2v1c + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv460:O2^- + CO2v1a + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv461:O2^- + CO2v1b + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv462:O2^- + CO2v1c + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv463:O2^- + CO2v1a + O2 => CO4^- + O2  k= 4.7e-41
% Rv464:O2^- + CO2v1b + O2 => CO4^- + O2  k= 4.7e-41
% Rv465:O2^- + CO2v1c + O2 => CO4^- + O2  k= 4.7e-41
% Rv466:O2^- + CO2 + CO2v1a => CO4^- + CO2v1a   k= 1.0e-41
% Rv467:O2^- + CO2 + CO2v1b => CO4^- + CO2v1b   k= 1.0e-41
% Rv468:O2^- + CO2 + CO2v1c => CO4^- + CO2v1c   k= 1.0e-41
% Rv469:O + O2 + CO2v1a => O3 + CO2v1a     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv470:O + O2 + CO2v1b => O3 + CO2v1b     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv471:O + O2 + CO2v1c => O3 + CO2v1c     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv472:O^- + CO2v1a => O + CO2v1a + e    k= 4.0e-18
% Rv473:O^- + CO2v1b => O + CO2v1b + e    k= 4.0e-18
% Rv474:O^- + CO2v1c => O + CO2v1c + e    k= 4.0e-18
% Rv475:O + O + CO2v1a => O2 + CO2v1a     k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv476:O + O + CO2v1b => O2 + CO2v1b     k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv477:O + O + CO2v1c => O2 + CO2v1c     k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv478:e + O2 + CO2e1 => O2^- + CO2e1      k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv479:e + O2 + CO2e2 => O2^- + CO2e2      k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv480:e + O2 + CO2va => O2^- + CO2va      k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv481:e + O2 + CO2vb => O2^- + CO2vb      k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv482:e + O2 + CO2vc => O2^- + CO2vc      k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv483:e + O2 + CO2vd => O2^- + CO2vd      k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv484:e + O2 + CO2v1 => O2^- + CO2v1      k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv485:e + O2 + CO2v2 => O2^- + CO2v2      k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv486:e + O2 + CO2v3 => O2^- + CO2v3      k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv487:e + O2 + CO2v4 => O2^- + CO2v4      k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv488:e + O2 + CO2v1a => O2^- + CO2v1a    k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv489:e + O2 + CO2v1b => O2^- + CO2v1b    k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv490:e + O2 + CO2v1c => O2^- + CO2v1c    k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv491:e + CO2 => e + CO2v5           k= 查表                                            EV=1.43eV
% Rv492:e + CO2 => e + CO2v6           k= 查表                                            EV=1.70eV
% Rv493:e + CO2 => e + CO2v7           k= 查表                                            EV=1.97eV
% Rv494:e + CO2 => e + CO2v8           k= 查表                                            EV=2.24eV
% Rv495:e + CO2v4 => e + CO2v5         k= 查表                                            EV=0.29eV
% Rv496:e + CO2v5 => e + CO2v6         k= 查表                                            EV=0.27eV
% Rv497:e + CO2v6 => e + CO2v7         k= 查表                                            EV=0.27eV
% Rv498:e + CO2v7 => e + CO2v8         k= 查表                                            EV=0.27eV
% Rv499:e + CO2v5 => CO2^+ + 2e        k= K(1)                                             ER=13.8eV
% Rv500:e + CO2v6 => CO2^+ + 2e        k= K(1)                                             ER=13.8eV
% Rv501:e + CO2v7 => CO2^+ + 2e        k= K(1)                                             ER=13.8eV
% Rv502:e + CO2v8 => CO2^+ + 2e        k= K(1)                                             ER=13.8eV
% Rv503:e + CO2v5 => O2^+ + C + 2e     k= ?                                    
% Rv504:e + CO2v6 => O2^+ + C + 2e     k= ?                                   
% Rv505:e + CO2v7 => O2^+ + C + 2e     k= ?                                    
% Rv506:e + CO2v8 => O2^+ + C + 2e     k= ?                                    
% Rv507:e + CO2v5 => 2e + O + CO^+     k= 查表                                            ER=18.07eV                        
% Rv508:e + CO2v6 => 2e + O + CO^+     k= 查表                                            ER=17.8eV                   
% Rv509:e + CO2v7 => 2e + O + CO^+     k= 查表                                            ER=17.53eV                                 
% Rv510:e + CO2v8 => 2e + O + CO^+     k= 查表                                            ER=17.26eV
% Rv511:e + CO2v5 => 2e + CO + O^+     k= 查表                                            ER=17.67eV
% Rv512:e + CO2v6 => 2e + CO + O^+     k= 查表                                            ER=17.4eV
% Rv513:e + CO2v7 => 2e + CO + O^+     k= 查表                                            ER=17.13eV
% Rv514:e + CO2v8 => 2e + CO + O^+     k= 查表                                            ER=16.86eV
% Rv515:e + CO2v5 => CO + O + e        k= 查表                                            ER=10.03eV
% Rv516:e + CO2v6 => CO + O + e        k= 查表                                            ER=9.76eV
% Rv517:e + CO2v7 => CO + O + e        k= 查表                                            ER=9.49eV
% Rv518:e + CO2v8 => CO + O + e        k= 查表                                            ER=9.22eV
% Rv519:e + CO2v5 => CO + O-           k= 查表                                        附着过程无ER，下同
% Rv520:e + CO2v6 => CO + O-           k= 查表
% Rv521:e + CO2v7 => CO + O-           k= 查表                    
% Rv522:e + CO2v8 => CO + O-           k= 查表
% Rv523:CO2v5 + CO2 => CO2v4 + CO2     k= 见下文
% Rv524:CO2v5 + CO => CO2v4 + CO       k= 见下文
% Rv525:CO2v5 + O2 => CO2v4 + O2       k= 见下文
% Rv526:CO2v6 + CO2 => CO2v5 + CO2     k= 见下文
% Rv527:CO2v6 + CO => CO2v5 + CO       k= 见下文
% Rv528:CO2v6 + O2 => CO2v5 + O2       k= 见下文
% Rv529:CO2v7 + CO2 => CO2v6 + CO2     k= 见下文
% Rv530:CO2v7 + CO => CO2v6 + CO       k= 见下文
% Rv531:CO2v7 + O2 => CO2v6 + O2       k= 见下文
% Rv532:CO2v8 + CO2 => CO2v7 + CO2     k= 见下文
% Rv533:CO2v8 + CO => CO2v7 + CO       k= 见下文
% Rv534:CO2v8 + O2 => CO2v7 + O2       k= 见下文
% Rv535:CO2v5 + CO2 => CO2v4 + CO2va   k= 2.13e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3))
% Rv536:CO2v5 + CO2 => CO2v4 + CO2vb   k= 2.13e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3))
% Rv537:CO2v6 + CO2 => CO2v5 + CO2va   k= 2.13e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3))
% Rv538:CO2v6 + CO2 => CO2v5 + CO2vb   k= 2.13e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3))
% Rv539:CO2v7 + CO2 => CO2v6 + CO2va   k= 2.13e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3))
% Rv540:CO2v7 + CO2 => CO2v6 + CO2vb   k= 2.13e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3))
% Rv541:CO2v8 + CO2 => CO2v7 + CO2va   k= 2.13e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3))
% Rv542:CO2v8 + CO2 => CO2v7 + CO2vb   k= 2.13e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3))
% Rv543:CO2v1 + CO2v4 => CO2 + CO2v5   k= 7.199e-17*exp(18.4*(Tg*1.16e4)^(-1/3)-76*(Tg*1.16e4)^(-2/3))
% Rv544:CO2v2 + CO2v4 => CO2v1 + CO2v5 k= 4.825e-17*exp(20.1*(Tg*1.16e4)^(-1/3)-65*(Tg*1.16e4)^(-2/3))
% Rv545:CO2v2 + CO2v5 => CO2v1 + CO2v6 k= 7.199e-17*exp(18.4*(Tg*1.16e4)^(-1/3)-76*(Tg*1.16e4)^(-2/3))
% Rv546:CO2v3 + CO2v4 => CO2v2 + CO2v5 k= 2.927e-17*exp(21.4*(Tg*1.16e4)^(-1/3)-53*(Tg*1.16e4)^(-2/3))
% Rv547:CO2v3 + CO2v5 => CO2v2 + CO2v6 k= 4.825e-17*exp(20.1*(Tg*1.16e4)^(-1/3)-65*(Tg*1.16e4)^(-2/3))
% Rv548:CO2v3 + CO2v6 => CO2v2 + CO2v7 k= 7.199e-17*exp(18.4*(Tg*1.16e4)^(-1/3)-76*(Tg*1.16e4)^(-2/3))
% Rv549:CO2v4 + CO2v4 => CO2v3 + CO2v5 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv550:CO2v4 + CO2v5 => CO2v3 + CO2v6 k= 2.927e-17*exp(21.4*(Tg*1.16e4)^(-1/3)-53*(Tg*1.16e4)^(-2/3))
% Rv551:CO2v4 + CO2v6 => CO2v3 + CO2v7 k= 4.825e-17*exp(20.1*(Tg*1.16e4)^(-1/3)-65*(Tg*1.16e4)^(-2/3))
% Rv552:CO2v4 + CO2v7 => CO2v3 + CO2v8 k= 7.199e-17*exp(18.4*(Tg*1.16e4)^(-1/3)-76*(Tg*1.16e4)^(-2/3))
% Rv553:CO2v5 + CO2v5 => CO2v4 + CO2v6 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv554:CO2v5 + CO2v6 => CO2v4 + CO2v7 k= 2.927e-17*exp(21.4*(Tg*1.16e4)^(-1/3)-53*(Tg*1.16e4)^(-2/3))
% Rv555:CO2v5 + CO2v7 => CO2v4 + CO2v8 k= 4.825e-17*exp(20.1*(Tg*1.16e4)^(-1/3)-65*(Tg*1.16e4)^(-2/3))
% Rv556:CO2v6 + CO2v6 => CO2v5 + CO2v7 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv557:CO2v6 + CO2v7 => CO2v5 + CO2v8 k= 2.927e-17*exp(21.4*(Tg*1.16e4)^(-1/3)-53*(Tg*1.16e4)^(-2/3))
% Rv558:CO2v7 + CO2v7 => CO2v6 + CO2v8 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv559:CO2v1 + CO2v5 => CO2v2 + CO2v4 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv560:CO2v1 + CO2v6 => CO2v2 + CO2v5 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv561:CO2v1 + CO2v7 => CO2v2 + CO2v6 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv562:CO2v1 + CO2v8 => CO2v2 + CO2v7 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv563:CO2v2 + CO2v5 => CO2v3 + CO2v4 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv564:CO2v2 + CO2v6 => CO2v3 + CO2v5 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv565:CO2v2 + CO2v7 => CO2v3 + CO2v6 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv566:CO2v2 + CO2v8 => CO2v3 + CO2v7 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv567:CO2v3 + CO2v5 => CO2v4 + CO2v4 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv568:CO2v3 + CO2v6 => CO2v4 + CO2v5 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv569:CO2v3 + CO2v7 => CO2v4 + CO2v6 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv570:CO2v3 + CO2v8 => CO2v4 + CO2v7 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv571:CO2v4 + CO2v6 => CO2v5 + CO2v5 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv572:CO2v4 + CO2v7 => CO2v5 + CO2v6 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv573:CO2v4 + CO2v8 => CO2v5 + CO2v7 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv574:CO2v5 + CO2v7 => CO2v6 + CO2v6 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv575:CO2v5 + CO2v8 => CO2v6 + CO2v7 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv576:CO2v6 + CO2v8 => CO2v7 + CO2v7 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv577:CO2 + CO2v5 => CO + O + CO2v5    k= 3.91e-16*exp(-36160/(Tg*11600))
% Rv578:CO2 + CO2v6 => CO + O + CO2v6    k= 3.91e-16*exp(-33654/(Tg*11600))
% Rv579:CO2 + CO2v7 => CO + O + CO2v7    k= 3.91e-16*exp(-31148/(Tg*11600))
% Rv580:CO2 + CO2v8 => CO + O + CO2v8    k= 3.91e-16*exp(-28643/(Tg*11600))
% Rv581:CO2 + CO2v5 => CO + O + CO2      k= 3.91e-16*exp(-36160/(Tg*11600))
% Rv582:CO2 + CO2v6 => CO + O + CO2      k= 3.91e-16*exp(-33654/(Tg*11600))
% Rv583:CO2 + CO2v7 => CO + O + CO2      k= 3.91e-16*exp(-31148/(Tg*11600))
% Rv584:CO2 + CO2v8 => CO + O + CO2      k= 3.91e-16*exp(-28643/(Tg*11600))
% Rv585:CO2v5 + O => CO + O2             k= 2.8e-17*exp(-18206/(Tg*11600))
% Rv586:CO2v6 + O => CO + O2             k= 2.8e-17*exp(-16640/(Tg*11600))
% Rv587:CO2v7 + O => CO + O2             k= 2.8e-17*exp(-15074/(Tg*11600))
% Rv588:CO2v8 + O => CO + O2             k= 2.8e-17*exp(-13508/(Tg*11600))
% Rv589:CO2v5 + C => 2CO                 k= 1.0e-21
% Rv590:CO2v6 + C => 2CO                 k= 1.0e-21
% Rv591:CO2v7 + C => 2CO                 k= 1.0e-21
% Rv592:CO2v8 + C => 2CO                 k= 1.0e-21
% Rv593:O^- + CO2v5 + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv594:O^- + CO2v6 + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv595:O^- + CO2v7 + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv596:O^- + CO2v8 + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv597:O^- + CO2v5 + CO => CO3^- + CO   k= 1.5e-40
% Rv598:O^- + CO2v6 + CO => CO3^- + CO   k= 1.5e-40
% Rv599:O^- + CO2v7 + CO => CO3^- + CO   k= 1.5e-40
% Rv600:O^- + CO2v8 + CO => CO3^- + CO   k= 1.5e-40
% Rv601:O^- + CO2v5 + O2 => CO3^- + O2   k= 3.1e-40
% Rv602:O^- + CO2v6 + O2 => CO3^- + O2   k= 3.1e-40
% Rv603:O^- + CO2v7 + O2 => CO3^- + O2   k= 3.1e-40
% Rv604:O^- + CO2v8 + O2 => CO3^- + O2   k= 3.1e-40
% Rv605:O^- + CO2 + CO2v5 => CO3^- + CO2v5 k= 9.0e-41
% Rv606:O^- + CO2 + CO2v6 => CO3^- + CO2v6 k= 9.0e-41
% Rv607:O^- + CO2 + CO2v7 => CO3^- + CO2v7 k= 9.0e-41
% Rv608:O^- + CO2 + CO2v8 => CO3^- + CO2v8 k= 9.0e-41
% Rv609:e + O2^+ + CO2v5 => O2 + CO2v5   k= 1.0e-38
% Rv610:e + O2^+ + CO2v6 => O2 + CO2v6   k= 1.0e-38
% Rv611:e + O2^+ + CO2v7 => O2 + CO2v7   k= 1.0e-38
% Rv612:e + O2^+ + CO2v8 => O2 + CO2v8   k= 1.0e-38
% Rv613:CO2v5 + CO^+ => CO2^+ + CO       k= 1.0e-15
% Rv614:CO2v6 + CO^+ => CO2^+ + CO       k= 1.0e-15
% Rv615:CO2v7 + CO^+ => CO2^+ + CO       k= 1.0e-15
% Rv616:CO2v8 + CO^+ => CO2^+ + CO       k= 1.0e-15
% Rv617:O3^- + CO2v5 => CO3^- + O2       k= 5.5e-16
% Rv618:O3^- + CO2v6 => CO3^- + O2       k= 5.5e-16
% Rv619:O3^- + CO2v7 => CO3^- + O2       k= 5.5e-16
% Rv620:O3^- + CO2v8 => CO3^- + O2       k= 5.5e-16
% Rv621:O^+ + CO2v5 => O2^+ + CO         k= 9.4e-16
% Rv622:O^+ + CO2v6 => O2^+ + CO         k= 9.4e-16
% Rv623:O^+ + CO2v7 => O2^+ + CO         k= 9.4e-16
% Rv624:O^+ + CO2v8 => O2^+ + CO         k= 9.4e-16
% Rv625:O^+ + CO2v5 => O + CO2^+         k= 4.5e-16
% Rv626:O^+ + CO2v6 => O + CO2^+         k= 4.5e-16
% Rv627:O^+ + CO2v7 => O + CO2^+         k= 4.5e-16
% Rv628:O^+ + CO2v8 => O + CO2^+         k= 4.5e-16
% Rv629:C^+ + CO2v5 => CO^+ + CO         k= 1.1e-15
% Rv630:C^+ + CO2v6 => CO^+ + CO         k= 1.1e-15
% Rv631:C^+ + CO2v7 => CO^+ + CO         k= 1.1e-15
% Rv632:C^+ + CO2v8 => CO^+ + CO         k= 1.1e-15
% Rv633:O2^+ + CO2v5 + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv634:O2^+ + CO2v6 + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv635:O2^+ + CO2v7 + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv636:O2^+ + CO2v8 + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv637:O2^- + CO2v5 + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv638:O2^- + CO2v6 + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv639:O2^- + CO2v7 + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv640:O2^- + CO2v8 + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv641:O2^- + CO2v5 + O2 => CO4^- + O2  k= 4.7e-41
% Rv642:O2^- + CO2v6 + O2 => CO4^- + O2  k= 4.7e-41
% Rv643:O2^- + CO2v7 + O2 => CO4^- + O2  k= 4.7e-41
% Rv644:O2^- + CO2v8 + O2 => CO4^- + O2  k= 4.7e-41
% Rv645:O2^- + CO2 + CO2v5 => CO4^- + CO2v5   k= 1.0e-41
% Rv646:O2^- + CO2 + CO2v6 => CO4^- + CO2v6   k= 1.0e-41
% Rv647:O2^- + CO2 + CO2v7 => CO4^- + CO2v7   k= 1.0e-41
% Rv648:O2^- + CO2 + CO2v8 => CO4^- + CO2v8   k= 1.0e-41
% Rv649:O + O2 + CO2v5 => O3 + CO2v5     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv650:O + O2 + CO2v6 => O3 + CO2v6     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv651:O + O2 + CO2v7 => O3 + CO2v7     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv652:O + O2 + CO2v8 => O3 + CO2v8     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv653:O^- + CO2v5 => O + CO2v5 + e     k= 4.0e-18
% Rv654:O^- + CO2v6 => O + CO2v6 + e     k= 4.0e-18
% Rv655:O^- + CO2v7 => O + CO2v7 + e     k= 4.0e-18
% Rv656:O^- + CO2v8 => O + CO2v8 + e     k= 4.0e-18
% Rv657:O + O + CO2v5 => O2 + CO2v5     k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv658:O + O + CO2v6 => O2 + CO2v6     k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv659:O + O + CO2v7 => O2 + CO2v7     k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv660:O + O + CO2v8 => O2 + CO2v8     k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv661:e + O2 + CO2v5 => O2^- + CO2v5      k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv662:e + O2 + CO2v6 => O2^- + CO2v6      k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv663:e + O2 + CO2v7 => O2^- + CO2v7      k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv664:e + O2 + CO2v8 => O2^- + CO2v8      k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv665:e + CO2 => e + CO2vk                k= 查表                                         EV=2.5eV
% Rv666:e + CO2vk => CO2^+ + 2e             k= K(1)                                         ER=13.8eV
% Rv667:e + CO2vk => O2^+ + C + 2e          k= ?
% Rv668:e + CO2vk => 2e + O + CO^+          k= 查表                                         ER=17.0eV
% Rv669:e + CO2vk => 2e + CO + O^+          k= 查表                                         ER=16.6eV
% Rv670:e + CO2vk => CO + O + e             k= 查表                                         ER=8.96eV
% Rv671:e + CO2vk => CO + O-                k= 查表
% Rv672:CO2v1 + CO2 => CO2vk + CO2      k= 1.43e-11*exp(-252*(Tg*11600)^(-1/3)+685*(Tg*11600)^(-2/3))
% Rv673:CO2v1 + CO => CO2vk + CO        k= 4.29e-12*exp(-252*(Tg*11600)^(-1/3)+685*(Tg*11600)^(-2/3))
% Rv674:CO2v1 + O2 => CO2vk + O2        k= 5.72e-12*exp(-252*(Tg*11600)^(-1/3)+685*(Tg*11600)^(-2/3))
% Rv675:CO2vk + CO2 => CO2va + CO2      k= 1.071e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv676:CO2vk + CO => CO2va + CO        k= 3.213e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv677:CO2vk + O2 => CO2va + O2        k= 3.213e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv678:CO2vk + CO2 => CO2vb + CO2      k= 2.897e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv679:CO2vk + CO => CO2vb + CO        k= 2.028e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv680:CO2vk + O2 => CO2vb + O2        k= 2.028e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv681:CO2 + CO2vk => CO + O + CO2vk   k= 3.91e-16*exp(-26230/(Tg*11600))
% Rv682:CO2 + CO2vk => CO + O + CO2     k= 3.91e-16*exp(-26230/(Tg*11600))
% Rv683:CO2vk + O => CO + O2            k= 2.8e-17*exp(-12000/(Tg*11600))
% Rv684:CO2vk + C => 2CO                k= 1.0e-21
% Rv685:O^- + CO2vk + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv686:O^- + CO2vk + CO => CO3^- + CO   k= 1.5e-40
% Rv687:O^- + CO2vk + O2 => CO3^- + O2   k= 3.1e-40
% Rv688:O^- + CO2 + CO2vk => CO3^- + CO2vk k= 9.0e-41
% Rv689:e + O2^+ + CO2vk => O2 + CO2vk   k= 1.0e-38
% Rv690:CO2vk + CO^+ => CO2^+ + CO       k= 1.0e-15
% Rv691:O3^- + CO2vk => CO3^- + O2       k= 5.5e-16
% Rv692:O^+ + CO2vk => O2^+ + CO         k= 9.4e-16
% Rv693:O^+ + CO2vk => O + CO2^+         k= 4.5e-16
% Rv694:C^+ + CO2vk => CO^+ + CO         k= 1.1e-15
% Rv695:O2^+ + CO2vk + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv696:O2^- + CO2vk + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv697:O2^- + CO2vk + O2 => CO4^- + O2  k= 4.7e-41
% Rv698:O2^- + CO2 + CO2vk => CO4^- + CO2vk   k= 1.0e-41
% Rv699:O + O2 + CO2vk => O3 + CO2vk     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv700:O^- + CO2vk => O + CO2vk + e     k= 4.0e-18
% Rv701:O + O + CO2vk => O2 + CO2vk      k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv702:e + O2 + CO2vk => O2^- + CO2vk   k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv703:CO2vk + CO2 => CO2va + CO2vb    k= 5.305e-15*exp(-88*(Tg*11600)^(-1/3)+233*(Tg*11600)^(-2/3))
% Rv704:CO2^+ + CO2e1 + CO => C2O4^+ + CO  k= 3.0e-40
% Rv705:CO2^+ + CO2e2 + CO => C2O4^+ + CO  k= 3.0e-40
% Rv706:CO2^+ + CO2va + CO => C2O4^+ + CO  k= 3.0e-40
% Rv707:CO2^+ + CO2vb + CO => C2O4^+ + CO  k= 3.0e-40
% Rv708:CO2^+ + CO2vc + CO => C2O4^+ + CO  k= 3.0e-40
% Rv709:CO2^+ + CO2vd + CO => C2O4^+ + CO  k= 3.0e-40
% Rv710:CO2^+ + CO2v1 + CO => C2O4^+ + CO  k= 3.0e-40
% Rv711:CO2^+ + CO2v2 + CO => C2O4^+ + CO  k= 3.0e-40
% Rv712:CO2^+ + CO2v3 + CO => C2O4^+ + CO  k= 3.0e-40
% Rv713:CO2^+ + CO2v4 + CO => C2O4^+ + CO  k= 3.0e-40
% Rv714:CO2^+ + CO2v1a + CO => C2O4^+ + CO  k= 3.0e-40
% Rv715:CO2^+ + CO2v1b + CO => C2O4^+ + CO  k= 3.0e-40
% Rv716:CO2^+ + CO2v1c + CO => C2O4^+ + CO  k= 3.0e-40
% Rv717:CO2^+ + CO2v5 + CO => C2O4^+ + CO  k= 3.0e-40
% Rv718:CO2^+ + CO2v6 + CO => C2O4^+ + CO  k= 3.0e-40
% Rv719:CO2^+ + CO2v7 + CO => C2O4^+ + CO  k= 3.0e-40
% Rv720:CO2^+ + CO2v8 + CO => C2O4^+ + CO  k= 3.0e-40
% Rv721:CO2^+ + CO2vk + CO => C2O4^+ + CO  k= 3.0e-40
% Rv722:C2O3^+ + CO + CO2e1 => CO2 + C2O2^+ + CO2e1  k= 2.6e-38
% Rv723:C2O3^+ + CO + CO2e2 => CO2 + C2O2^+ + CO2e2  k= 2.6e-38
% Rv724:C2O3^+ + CO + CO2va => CO2 + C2O2^+ + CO2va  k= 2.6e-38
% Rv725:C2O3^+ + CO + CO2vb => CO2 + C2O2^+ + CO2vb  k= 2.6e-38
% Rv726:C2O3^+ + CO + CO2vc => CO2 + C2O2^+ + CO2vc  k= 2.6e-38
% Rv727:C2O3^+ + CO + CO2vd => CO2 + C2O2^+ + CO2vd  k= 2.6e-38
% Rv728:C2O3^+ + CO + CO2v1 => CO2 + C2O2^+ + CO2v1  k= 2.6e-38
% Rv729:C2O3^+ + CO + CO2v2 => CO2 + C2O2^+ + CO2v2  k= 2.6e-38
% Rv730:C2O3^+ + CO + CO2v3 => CO2 + C2O2^+ + CO2v3  k= 2.6e-38
% Rv731:C2O3^+ + CO + CO2v4 => CO2 + C2O2^+ + CO2v4  k= 2.6e-38
% Rv732:C2O3^+ + CO + CO2v1a => CO2 + C2O2^+ + CO2v1a  k= 2.6e-38
% Rv733:C2O3^+ + CO + CO2v1b => CO2 + C2O2^+ + CO2v1b  k= 2.6e-38
% Rv734:C2O3^+ + CO + CO2v1c => CO2 + C2O2^+ + CO2v1c  k= 2.6e-38
% Rv735:C2O3^+ + CO + CO2v5 => CO2 + C2O2^+ + CO2v5  k= 2.6e-38
% Rv736:C2O3^+ + CO + CO2v6 => CO2 + C2O2^+ + CO2v6  k= 2.6e-38
% Rv737:C2O3^+ + CO + CO2v7 => CO2 + C2O2^+ + CO2v7  k= 2.6e-38
% Rv738:C2O3^+ + CO + CO2v8 => CO2 + C2O2^+ + CO2v8  k= 2.6e-38
% Rv739:C2O3^+ + CO + CO2vk => CO2 + C2O2^+ + CO2vk  k= 2.6e-38
% Rv740:C2O4^+ + CO + CO2e1 => CO2 + C2O3^+ + CO2e1  k= 4.2e-38
% Rv741:C2O4^+ + CO + CO2e2 => CO2 + C2O3^+ + CO2e2  k= 4.2e-38
% Rv742:C2O4^+ + CO + CO2va => CO2 + C2O3^+ + CO2va  k= 4.2e-38
% Rv743:C2O4^+ + CO + CO2vb => CO2 + C2O3^+ + CO2vb  k= 4.2e-38
% Rv744:C2O4^+ + CO + CO2vc => CO2 + C2O3^+ + CO2vc  k= 4.2e-38
% Rv745:C2O4^+ + CO + CO2vd => CO2 + C2O3^+ + CO2vd  k= 4.2e-38
% Rv746:C2O4^+ + CO + CO2v1 => CO2 + C2O3^+ + CO2v1  k= 4.2e-38
% Rv747:C2O4^+ + CO + CO2v2 => CO2 + C2O3^+ + CO2v2  k= 4.2e-38
% Rv748:C2O4^+ + CO + CO2v3 => CO2 + C2O3^+ + CO2v3  k= 4.2e-38
% Rv749:C2O4^+ + CO + CO2v4 => CO2 + C2O3^+ + CO2v4  k= 4.2e-38
% Rv750:C2O4^+ + CO + CO2v1a => CO2 + C2O3^+ + CO2v1a  k= 4.2e-38
% Rv751:C2O4^+ + CO + CO2v1b => CO2 + C2O3^+ + CO2v1b  k= 4.2e-38
% Rv752:C2O4^+ + CO + CO2v1c => CO2 + C2O3^+ + CO2v1c  k= 4.2e-38
% Rv753:C2O4^+ + CO + CO2v5 => CO2 + C2O3^+ + CO2v5  k= 4.2e-38
% Rv754:C2O4^+ + CO + CO2v6 => CO2 + C2O3^+ + CO2v6  k= 4.2e-38
% Rv755:C2O4^+ + CO + CO2v7 => CO2 + C2O3^+ + CO2v7  k= 4.2e-38
% Rv756:C2O4^+ + CO + CO2v8 => CO2 + C2O3^+ + CO2v8  k= 4.2e-38
% Rv757:C2O4^+ + CO + CO2vk => CO2 + C2O3^+ + CO2vk  k= 4.2e-38
% Rv758:C2O2^+ + CO2e1 => CO^+ + CO + CO2e1    k= 1.0e-18
% Rv759:C2O2^+ + CO2e2 => CO^+ + CO + CO2e2    k= 1.0e-18
% Rv760:C2O2^+ + CO2va => CO^+ + CO + CO2va    k= 1.0e-18
% Rv761:C2O2^+ + CO2vb => CO^+ + CO + CO2vb    k= 1.0e-18
% Rv762:C2O2^+ + CO2vc => CO^+ + CO + CO2vc    k= 1.0e-18
% Rv763:C2O2^+ + CO2vd => CO^+ + CO + CO2vd    k= 1.0e-18
% Rv764:C2O2^+ + CO2v1 => CO^+ + CO + CO2v1    k= 1.0e-18
% Rv765:C2O2^+ + CO2v2 => CO^+ + CO + CO2v2    k= 1.0e-18
% Rv766:C2O2^+ + CO2v3 => CO^+ + CO + CO2v3    k= 1.0e-18
% Rv767:C2O2^+ + CO2v4 => CO^+ + CO + CO2v4    k= 1.0e-18
% Rv768:C2O2^+ + CO2v1a => CO^+ + CO + CO2v1a    k= 1.0e-18
% Rv769:C2O2^+ + CO2v1b => CO^+ + CO + CO2v1b    k= 1.0e-18
% Rv770:C2O2^+ + CO2v1c => CO^+ + CO + CO2v1c    k= 1.0e-18
% Rv771:C2O2^+ + CO2v5 => CO^+ + CO + CO2v5    k= 1.0e-18
% Rv772:C2O2^+ + CO2v6 => CO^+ + CO + CO2v6    k= 1.0e-18
% Rv773:C2O2^+ + CO2v7 => CO^+ + CO + CO2v7    k= 1.0e-18
% Rv774:C2O2^+ + CO2v8 => CO^+ + CO + CO2v8    k= 1.0e-18
% Rv775:C2O2^+ + CO2vk => CO^+ + CO + CO2vk    k= 1.0e-18
% Rv776:C2O4^+ + CO2e1 => CO2^+ + CO2 + CO2e1  k= 1.0e-20
% Rv777:C2O4^+ + CO2e2 => CO2^+ + CO2 + CO2e2  k= 1.0e-20
% Rv778:C2O4^+ + CO2va => CO2^+ + CO2 + CO2va  k= 1.0e-20
% Rv779:C2O4^+ + CO2vb => CO2^+ + CO2 + CO2vb  k= 1.0e-20
% Rv780:C2O4^+ + CO2vc => CO2^+ + CO2 + CO2vc  k= 1.0e-20
% Rv781:C2O4^+ + CO2vd => CO2^+ + CO2 + CO2vd  k= 1.0e-20
% Rv782:C2O4^+ + CO2v1 => CO2^+ + CO2 + CO2v1  k= 1.0e-20
% Rv783:C2O4^+ + CO2v2 => CO2^+ + CO2 + CO2v2  k= 1.0e-20
% Rv784:C2O4^+ + CO2v3 => CO2^+ + CO2 + CO2v3  k= 1.0e-20
% Rv785:C2O4^+ + CO2v4 => CO2^+ + CO2 + CO2v4  k= 1.0e-20
% Rv786:C2O4^+ + CO2v1a => CO2^+ + CO2 + CO2v1a  k= 1.0e-20
% Rv787:C2O4^+ + CO2v1b => CO2^+ + CO2 + CO2v1b  k= 1.0e-20
% Rv788:C2O4^+ + CO2v1c => CO2^+ + CO2 + CO2v1c  k= 1.0e-20
% Rv789:C2O4^+ + CO2v5 => CO2^+ + CO2 + CO2v5  k= 1.0e-20
% Rv790:C2O4^+ + CO2v6 => CO2^+ + CO2 + CO2v6  k= 1.0e-20
% Rv791:C2O4^+ + CO2v7 => CO2^+ + CO2 + CO2v7  k= 1.0e-20
% Rv792:C2O4^+ + CO2v8 => CO2^+ + CO2 + CO2v8  k= 1.0e-20
% Rv793:C2O4^+ + CO2vk => CO2^+ + CO2 + CO2vk  k= 1.0e-20
% Rv794:N2^+ + CO2e1 => CO2^+ + N2       k= K(125)*0.916
% Rv795:N2^+ + CO2e2 => CO2^+ + N2       k= K(125)*0.093
% Rv796:N2^+ + CO2va => CO2^+ + N2       k= 9.0e-16
% Rv797:N2^+ + CO2vb => CO2^+ + N2       k= 9.0e-16
% Rv798:N2^+ + CO2vc => CO2^+ + N2       k= 9.0e-16
% Rv799:N2^+ + CO2vd => CO2^+ + N2       k= 9.0e-16
% Rv800:N2^+ + CO2v1 => CO2^+ + N2       k= 9.0e-16
% Rv801:N2^+ + CO2v2 => CO2^+ + N2       k= 9.0e-16
% Rv802:N2^+ + CO2v3 => CO2^+ + N2       k= 9.0e-16
% Rv803:N2^+ + CO2v4 => CO2^+ + N2       k= 9.0e-16
% Rv804:N2^+ + CO2v1a => CO2^+ + N2       k= 9.0e-16
% Rv805:N2^+ + CO2v1b => CO2^+ + N2       k= 9.0e-16
% Rv806:N2^+ + CO2v1c => CO2^+ + N2       k= 9.0e-16
% Rv807:N2^+ + CO2v5 => CO2^+ + N2       k= 9.0e-16
% Rv808:N2^+ + CO2v6 => CO2^+ + N2       k= 9.0e-16
% Rv809:N2^+ + CO2v7 => CO2^+ + N2       k= 9.0e-16
% Rv810:N2^+ + CO2v8 => CO2^+ + N2       k= 9.0e-16
% Rv811:N2^+ + CO2vk => CO2^+ + N2       k= 9.0e-16
% Rv812:N + CO2e1 => NO + CO             k= 1.0e-25
% Rv813:N + CO2e2 => NO + CO             k= 1.0e-25
% Rv814:N + CO2va => NO + CO             k= 1.0e-25
% Rv815:N + CO2vb => NO + CO             k= 1.0e-25
% Rv816:N + CO2vc => NO + CO             k= 1.0e-25
% Rv817:N + CO2vd => NO + CO             k= 1.0e-25
% Rv818:N + CO2v1 => NO + CO             k= 1.0e-25
% Rv819:N + CO2v2 => NO + CO             k= 1.0e-25
% Rv820:N + CO2v3 => NO + CO             k= 1.0e-25
% Rv821:N + CO2v4 => NO + CO             k= 1.0e-25
% Rv822:N + CO2v1a => NO + CO             k= 1.0e-25
% Rv823:N + CO2v1b => NO + CO             k= 1.0e-25
% Rv824:N + CO2v1c => NO + CO             k= 1.0e-25
% Rv825:N + CO2v5 => NO + CO             k= 1.0e-25
% Rv826:N + CO2v6 => NO + CO             k= 1.0e-25
% Rv827:N + CO2v7 => NO + CO             k= 1.0e-25
% Rv828:N + CO2v8 => NO + CO             k= 1.0e-25
% Rv829:N + CO2vk => NO + CO             k= 1.0e-25
% Rv830:NO + O + CO2e1 => CO2e1 + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv831:NO + O + CO2e2 => CO2e2 + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv832:NO + O + CO2va => CO2va + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv833:NO + O + CO2vb => CO2vb + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv834:NO + O + CO2vc => CO2vc + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv835:NO + O + CO2vd => CO2vd + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv836:NO + O + CO2v1 => CO2v1 + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv837:NO + O + CO2v2 => CO2v2 + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv838:NO + O + CO2v3 => CO2v3 + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv839:NO + O + CO2v4 => CO2v4 + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv840:NO + O + CO2v1a => CO2v1a + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv841:NO + O + CO2v1b => CO2v1b + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv842:NO + O + CO2v1c => CO2v1c + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv843:NO + O + CO2v5 => CO2v5 + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv844:NO + O + CO2v6 => CO2v6 + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv845:NO + O + CO2v7 => CO2v7 + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv846:NO + O + CO2v8 => CO2v8 + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv847:NO + O + CO2vk => CO2vk + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv848:NO^- + CO2e1 => NO + CO2e1 + e   k= 8.3e-18
% Rv849:NO^- + CO2e2 => NO + CO2e2 + e   k= 8.3e-18
% Rv850:NO^- + CO2va => NO + CO2va + e   k= 8.3e-18
% Rv851:NO^- + CO2vb => NO + CO2vb + e   k= 8.3e-18
% Rv852:NO^- + CO2vc => NO + CO2vc + e   k= 8.3e-18
% Rv853:NO^- + CO2vd => NO + CO2vd + e   k= 8.3e-18
% Rv854:NO^- + CO2v1 => NO + CO2v1 + e   k= 8.3e-18
% Rv855:NO^- + CO2v2 => NO + CO2v2 + e   k= 8.3e-18
% Rv856:NO^- + CO2v3 => NO + CO2v3 + e   k= 8.3e-18
% Rv857:NO^- + CO2v4 => NO + CO2v4 + e   k= 8.3e-18
% Rv858:NO^- + CO2v1a => NO + CO2v1a + e   k= 8.3e-18
% Rv859:NO^- + CO2v1b => NO + CO2v1b + e   k= 8.3e-18
% Rv860:NO^- + CO2v1c => NO + CO2v1c + e   k= 8.3e-18
% Rv861:NO^- + CO2v5 => NO + CO2v5 + e   k= 8.3e-18
% Rv862:NO^- + CO2v6 => NO + CO2v6 + e   k= 8.3e-18
% Rv863:NO^- + CO2v7 => NO + CO2v7 + e   k= 8.3e-18
% Rv864:NO^- + CO2v8 => NO + CO2v8 + e   k= 8.3e-18
% Rv865:NO^- + CO2vk => NO + CO2vk + e   k= 8.3e-18

global species 

species{1}='e';              mass(1)=9.1e-31;
species{2}='CO';             mass(2)=4.65e-26;
species{3}='O_2';            mass(3)=5.3e-26;
species{4}='O';              mass(4)=2.66e-26;
species{5}='O_3';            mass(5)=7.97e-26;
species{6}='O^-';            mass(6)=2.66e-26;
species{7}='O_2^-';          mass(7)=5.3e-26;
species{8}='CO_2^+';         mass(8)=7.31e-26;
species{9}='CO_4^-';         mass(9)=1.261e-25;
species{10}='CO_2';          mass(10)=7.31e-26;
species{11}='C';             mass(11)=1.99e-26;
species{12}='O_2^+';         mass(12)=5.3e-26;
species{13}='CO_3^-';        mass(13)=9.97e-26;
species{14}='CO^+';          mass(14)=4.65e-26;
species{15}='O_3^-';         mass(15)=7.97e-26;
species{16}='O^+';           mass(16)=2.66e-26;
species{17}='C^+';           mass(17)=1.99e-26;
species{18}='CO_4^+';        mass(18)=1.261e-25;
species{19}='C_2O_2^+';      mass(19)=9.3e-26;
species{20}='C_2O_3^+';      mass(20)=1.196e-25;
species{21}='C_2O_4^+';      mass(21)=1.461e-25;
species{22}='N_2';           mass(22)=4.65e-26;
species{23}='NO';            mass(23)=4.99e-26;
species{24}='N_2O';          mass(24)=7.31e-26;
species{25}='NO_2';          mass(25)=7.65e-26;
species{26}='N';             mass(26)=2.325e-26;
species{27}='NO_3';          mass(27)=1.031e-25;
species{28}='N_2^+';         mass(28)=4.65e-26;
species{29}='NO^+';          mass(29)=4.99e-26;
species{30}='N_2O^+';        mass(30)=7.31e-26;
species{31}='NO_2^+';        mass(31)=7.65e-26;
species{32}='NO^-';          mass(32)=4.99e-26;
species{33}='NO_2^-';        mass(33)=7.65e-26;
species{34}='NO_3^-';        mass(34)=1.031e-25;
% if Vibration==1
    species{35}='CO_2e1';    mass(35)=7.31e-26;
    species{36}='CO_2e2';    mass(36)=7.31e-26;
    species{37}='CO_2va';    mass(37)=7.31e-26;
    species{38}='CO_2vb';    mass(38)=7.31e-26;
    species{39}='CO_2vc';    mass(39)=7.31e-26;
    species{40}='CO_2vd';    mass(40)=7.31e-26;
    species{41}='CO_2v1';    mass(41)=7.31e-26;
    species{42}='CO_2v2';    mass(42)=7.31e-26;
    species{43}='CO_2v3';    mass(43)=7.31e-26;
    species{44}='CO_2v4';    mass(44)=7.31e-26;
    species{45}='CO_2v1a';   mass(45)=7.31e-26;
    species{46}='CO_2v1b';   mass(46)=7.31e-26;
    species{47}='CO_2v1c';   mass(47)=7.31e-26;
    species{48}='CO_2v5';    mass(48)=7.31e-26;
    species{49}='CO_2v6';    mass(49)=7.31e-26;
    species{50}='CO_2v7';    mass(50)=7.31e-26;
    species{51}='CO_2v8';    mass(51)=7.31e-26;
    species{52}='CO_2vk';    mass(52)=7.31e-26;
% end

nsp=length(species);
number_of_reaction=178;
number_of_vreaction=865;
% if Vibration==1
%     number_of_vreaction=703;
% else
%     number_of_vreaction=0;
% end
number_of_reactions = number_of_reaction + number_of_vreaction;
K(1:number_of_reaction)=0;
Kv=zeros(1,number_of_vreaction);
r=zeros(1,8);%CO2vn
Fr=zeros(1,8);
auxv=[];

aux{1}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-e_CO2ionization.lut');
aux{2}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-e_CO2dissociation.lut');
aux{3}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-e_CO2attachment.lut');
aux{5}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-e_O2dissociation.lut');
aux{6}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-e_O2attachment.lut');
aux{13}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-e_CO2elastic.lut');
aux{17}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-e_COattachment.lut');
aux{18}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-e_COexcitation.lut');
aux{22}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-e_O2ionization.lut');
aux{34}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-e_O3attachment_O-.lut');
aux{37}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-e_O3attachment_O2-.lut');
aux{44}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-e_CO2ionization_CO+.lut');
aux{45}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-e_COionization.lut');
aux{53}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-e_CO2ionization_O+.lut');
aux{54}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-e_Oionization.lut');
aux{58}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-e_COionization_C+.lut');
aux{59}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-e_Cionization.lut');
aux{94}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-e_N2elastic.lut');
aux{95}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-e_N2ionization.lut');
aux{96}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-e_N2dissociation.lut');
aux{97}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-NOattachment.lut');
aux{98}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-NOionization.lut');
aux{99}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-N2Oattachment.lut');
aux{100}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-N2Oionization.lut');
if Vibration==1
    auxv{1}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2e1.lut');
    auxv{2}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2e2.lut');
    auxv{15}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2va.lut');
    auxv{16}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2vb.lut');
    auxv{17}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2vc.lut');
    auxv{18}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2vd.lut');
    auxv{31}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v1.lut');
    auxv{32}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v2.lut');
    auxv{33}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v3.lut');
    auxv{34}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v4.lut');
    auxv{39}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v1dissociation.lut');
    auxv{40}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v2dissociation.lut');
    auxv{41}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v3dissociation.lut');
    auxv{42}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v4dissociation.lut');
    auxv{43}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v1attachment.lut');
    auxv{44}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v2attachment.lut');
    auxv{45}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v3attachment.lut');
    auxv{46}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v4attachment.lut');
    auxv{51}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v1ionization_CO+.lut');
    auxv{52}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v2ionization_CO+.lut');
    auxv{53}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v3ionization_CO+.lut');
    auxv{54}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v4ionization_CO+.lut');
    auxv{55}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v1ionization_O+.lut');
    auxv{56}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v2ionization_O+.lut');
    auxv{57}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v3ionization_O+.lut');
    auxv{58}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v4ionization_O+.lut');
    auxv{103}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2vadissociation.lut');
    auxv{104}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2vbdissociation.lut');
    auxv{105}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2vcdissociation.lut');
    auxv{106}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2vddissociation.lut');
    auxv{107}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2vaattachment.lut');
    auxv{108}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2vbattachment.lut');
    auxv{109}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2vcattachment.lut');
    auxv{110}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2vdattachment.lut');
    auxv{119}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2vaionization_CO+.lut');
    auxv{120}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2vbionization_CO+.lut');
    auxv{121}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2vcionization_CO+.lut');
    auxv{122}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2vdionization_CO+.lut');
    auxv{123}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2vaionization_O+.lut');
    auxv{124}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2vbionization_O+.lut');
    auxv{125}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2vcionization_O+.lut');
    auxv{126}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2vdionization_O+.lut');
    auxv{334}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v2-3.lut');
    auxv{336}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v1a.lut');
    auxv{337}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v1b.lut');
    auxv{338}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v1c.lut');
    auxv{345}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v1aionization_CO+.lut');
    auxv{346}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v1bionization_CO+.lut');
    auxv{347}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v1cionization_CO+.lut');
    auxv{348}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v1aionization_O+.lut');
    auxv{349}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v1bionization_O+.lut');
    auxv{350}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v1cionization_O+.lut');
    auxv{351}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v1adissociation.lut');
    auxv{352}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v1bdissociation.lut');
    auxv{353}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v1cdissociation.lut');
    auxv{354}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v1aattachment.lut');
    auxv{355}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v1battachment.lut');
    auxv{356}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v1cattachment.lut');
    auxv{491}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v5.lut');
    auxv{492}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v6.lut');
    auxv{493}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v7.lut');
    auxv{494}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v8.lut');
    auxv{496}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v5-6.lut');
    auxv{507}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v5ionization_CO+.lut');
    auxv{508}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v6ionization_CO+.lut');
    auxv{509}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v7ionization_CO+.lut');
    auxv{510}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v8ionization_CO+.lut');
    auxv{511}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v5ionization_O+.lut');
    auxv{512}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v6ionization_O+.lut');
    auxv{513}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v7ionization_O+.lut');
    auxv{514}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v8ionization_O+.lut');
    auxv{515}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v5dissociation.lut');
    auxv{516}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v6dissociation.lut');
    auxv{517}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v7dissociation.lut');
    auxv{518}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v8dissociation.lut');
    auxv{519}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v5attachment.lut');
    auxv{520}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v6attachment.lut');
    auxv{521}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v7attachment.lut');
    auxv{522}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2v8attachment.lut');
    auxv{665}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2vk.lut');
    auxv{668}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2vkionization_CO+.lut');
    auxv{669}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2vkionization_O+.lut');
    auxv{670}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2vkdissociation.lut');
    auxv{671}=load('C:\zthmatlabhome\MATLAB2\reaction-rate-cross_section-CO2vkattachment.lut');
end
    
%% Function without Te
% 基态存在电子得失
K(7)=5.5e-16;
K(8)=3.0e-16;
K(12)=6.9e-16;
K(15)=4.0e-18; 
K(16)=2.3e-16;
K(29)=5.0e-19;
K(39)=3.3e-16;
K(48)=4.6e-40;
K(68)=5.0e-16;
K(71)=1.0e-38;
K(73)=2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4));
K(106)=1.6e-43;
K(107)=1.6e-43;
K(108)=2.2e-16;
K(109)=1.0e-18;
K(110)=2.5e-16;
K(111)=4.0e-16;
K(112)=1.0e-18;
K(113)=1.0e-18;
K(114)=1.0e-18;
K(115)=1.0e-18;
K(116)=5.0e-16;
K(177)=5.0e-19;
K(178)=8.3e-18;

%基态无电子得失
K(4)=9.0e-16;
K(10)=6.0e-13;
K(11)=1.7e-42*(Tg*1.16e4)^(-1.2);
K(14)=3.31e-16;
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
%K(70)=2.3e-41;
K(70)=0;
K(72)=1.0e-41;
K(77)=3.0e-40;
K(78)=1.1e-15;
K(79)=9.0e-16;
K(80)=2.6e-38;
K(81)=4.2e-38;
K(82)=5.0e-18;
K(83)=1.0e-18;
K(84)=5.0e-13;
K(85)=5.0e-13;
K(86)=6.0e-13;
K(87)=5.0e-13;
K(88)=5.0e-13;
K(89)=6.0e-13;
K(90)=1.0e-20;
K(91)=5.0e-13;
K(92)=5.0e-13;
K(93)=6.0e-13;
K(121)=1.2e-16;
K(122)=7.0e-17;
K(123)=4.9e-16;
K(124)=7.0e-17;
K(125)=9.0e-16;
K(126)=1.9e-17;
K(127)=7.0e-16;
K(128)=6.6e-16;
K(129)=1.0e-20;
K(130)=2.9e-16;
K(131)=2.9e-16;
K(132)=1.0e-15;
K(133)=1.0e-15;
K(134)=2.0e-16;
K(135)=2.0e-15;
K(136)=1.0e-17;
K(137)=5.0e-16;
K(138)=1.0e-17;
K(139)=2.6e-18;
K(140)=7.0e-16;
K(141)=2.8e-16;
K(142)=5.0e-16;
K(143)=9.0e-18;
K(144)=1.0e-16;
K(145)=5.0e-17;
K(146)=9.0e-16;
K(147)=1.8e-17;
K(148)=4.0e-18;
K(149)=3.0e-16;
K(150)=2.0e-17;
K(151)=3.0e-40*(Tg*11600)^(-1);
K(152)=6.0e-13;
K(153)=5.0e-13;
K(154)=6.0e-13;
K(155)=6.0e-13;
K(156)=5.1e-13;
K(157)=8.1e-13;
K(158)=5.8e-13;
K(159)=4.9e-13;
K(160)=4.1e-13;
K(161)=1.3e-13;
K(162)=3.0e-13;
K(163)=3.0e-13;
K(164)=4.1e-13;
K(165)=1.3e-13;
K(166)=4.2e-13;
K(167)=1.0e-13;
K(168)=6.0e-13;
K(169)=5.0e-13;
K(170)=6.0e-13;
K(171)=5.0e-13;
K(172)=6.0e-13;
K(173)=5.0e-13;
K(174)=1.0e-25;
K(175)=1.48e-16*exp(-16967/(Tg*1.16e4));
K(176)=1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4));

%振动态
if Vibration==1
    Kv(7)=K(46)*0.916;
    Kv(8)=K(46)*0.093;
    Kv(11)=K(56)*0.916;
    Kv(12)=K(56)*0.093;
    Kv(13)=K(60)*3.189;
    Kv(14)=K(60)*0.735;
    Kv(19)=7.14e-14*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
    Kv(20)=5.00e-14*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
    Kv(21)=5.00e-14*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
    Kv(22)=1.071e-15*exp(-137*(Tg*1.16e4)^(-1/3));
    Kv(23)=3.213e-15*exp(-137*(Tg*1.16e4)^(-1/3));
    Kv(24)=3.213e-15*exp(-137*(Tg*1.16e4)^(-1/3));
    Kv(25)=1.438e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
    Kv(26)=1.007e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
    Kv(27)=1.007e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
    Kv(28)=1.071e-15*exp(-137*(Tg*1.16e4)^(-1/3));
    Kv(29)=3.213e-15*exp(-137*(Tg*1.16e4)^(-1/3));
    Kv(30)=3.213e-15*exp(-137*(Tg*1.16e4)^(-1/3));
    Kv(59)=4.25e-7*exp(-407*(Tg*11600)^(-1/3)+824*(Tg*11600)^(-2/3));
    Kv(60)=Kv(59)*0.3;
    Kv(61)=Kv(59)*0.4;
    Kv(62)=8.57e-7*exp(-404*(Tg*11600)^(-1/3)+1096*(Tg*11600)^(-2/3));
    Kv(63)=Kv(62)*0.3;
    Kv(64)=Kv(62)*0.4;
    Kv(65)=1.43e-11*exp(-252*(Tg*11600)^(-1/3)+685*(Tg*11600)^(-2/3));
    Kv(66)=Kv(65)*0.3;
    Kv(67)=Kv(65)*0.4;
    %VT弛豫
    r(1,1)=10.555*Tg^(-1/2);%1→0,M=CO2
    r(2,1)=9.01*Tg^(-1/2);%1→0,M=CO
    r(3,1)=9.11*Tg^(-1/2);%1→0,M=O2
    r(1:3,2)=r(1:3,1);%2→1同1→0,均0.29
    r(1,3)=10.191*Tg^(-1/2);%3→2，0.28
    r(2,3)=8.7*Tg^(-1/2);
    r(3,3)=8.794*Tg^(-1/2);
    r(1:3,4)=r(1:3,3);
    r(1:3,5)=r(1:3,1);
    r(1,6)=9.827*Tg^(-1/2);
    r(2,6)=8.391*Tg^(-1/2);
    r(3,6)=8.48*Tg^(-1/2);
    r(1:3,7)=r(1:3,6);
    r(1:3,8)=r(1:3,6);
    Fr=0.5.*(3-exp(-2/3.*r)).*exp(-2/3.*r);
    Kv(68)=(2.0074*Kv(59)+2.002*Kv(62)+1.9697*Kv(65))*Fr(1,2)/Fr(1,1);
    Kv(69)=(2.0074*Kv(60)+2.002*Kv(63)+1.9697*Kv(66))*Fr(2,2)/Fr(2,1);
    Kv(70)=(2.0074*Kv(61)+2.002*Kv(64)+1.9697*Kv(67))*Fr(3,2)/Fr(3,1);
    Kv(71)=(3.0224*Kv(59)+3.006*Kv(62)+2.9106*Kv(65))*Fr(1,3)/Fr(1,1);
    Kv(72)=(3.0224*Kv(60)+3.006*Kv(63)+2.9106*Kv(66))*Fr(2,3)/Fr(2,1);
    Kv(73)=(3.0224*Kv(61)+3.006*Kv(64)+2.9106*Kv(67))*Fr(3,3)/Fr(3,1);
    Kv(74)=(4.0451*Kv(59)+4.012*Kv(62)+3.8238*Kv(65))*Fr(1,4)/Fr(1,1);
    Kv(75)=(4.0451*Kv(60)+4.012*Kv(63)+3.8238*Kv(66))*Fr(2,4)/Fr(2,1);
    Kv(76)=(4.0451*Kv(61)+4.012*Kv(64)+3.8238*Kv(67))*Fr(3,4)/Fr(3,1);
    Kv(77)=1.06e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3));
    Kv(78:83)=2.13e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3));
    Kv(84)=2.157e-15*exp(-88*(Tg*11600)^(-1/3)+233*(Tg*11600)^(-2/3));
    Kv(85)=5.305e-15*exp(-88.2*(Tg*11600)^(-1/3)+233*(Tg*11600)^(-2/3));
    Kv(86)=5.305e-15*exp(-88*(Tg*11600)^(-1/3)+233*(Tg*11600)^(-2/3));
    %V-N
    Kv(87)=3.91e-16*exp(-48688/(Tg*11600));
    Kv(88)=3.91e-16*exp(-47852/(Tg*11600));
    Kv(89)=3.91e-16*exp(-47110/(Tg*11600));
    Kv(90)=3.91e-16*exp(-46368/(Tg*11600));
    Kv(91)=3.91e-16*exp(-46739/(Tg*11600));
    Kv(92)=3.91e-16*exp(-44048/(Tg*11600));
    Kv(93)=3.91e-16*exp(-41449/(Tg*11600));
    Kv(94)=3.91e-16*exp(-38851/(Tg*11600));
    Kv(95)=2.8e-17*exp(-26036/(Tg*11600));
    Kv(96)=2.8e-17*exp(-25514/(Tg*11600));
    Kv(97)=2.8e-17*exp(-25050/(Tg*11600));
    Kv(98)=2.8e-17*exp(-24586/(Tg*11600));
    Kv(99)=2.8e-17*exp(-24818/(Tg*11600));
    Kv(100)=2.8e-17*exp(-23136/(Tg*11600));
    Kv(101)=2.8e-17*exp(-21512/(Tg*11600));
    Kv(102)=2.8e-17*exp(-19888/(Tg*11600));
    %离子反应补充
    Kv(127:136)=1.0e-21;
    Kv(137:146)=9.0e-41;
    Kv(147:156)=1.5e-40;
    Kv(157:166)=3.1e-40;
    Kv(167:174)=1.0e-15;
    Kv(175:184)=5.5e-16;
    Kv(185:194)=9.4e-16;
    Kv(195:202)=4.5e-16;
    Kv(203:210)=1.1e-15;
    Kv(211:220)=4.7e-41;
    %
    Kv(225:234)=1.7e-42*(Tg*1.16e4)^(-1.2);
    Kv(235:244)=4.0e-18;
    Kv(245:254)=9.0e-41;
    Kv(255:264)=3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4));
    %Kv(265:274)=2.3e-41;
    Kv(265:274)=0;
    Kv(275:284)=1.0e-38;
    Kv(285:304)=1.0e-41;
    Kv(305)=2.897e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
    Kv(306)=2.028e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
    Kv(307)=2.028e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
    Kv(308)=1.582e-11*exp(-272*(Tg*1.16e4)^(-1/3)+437*(Tg*1.16e4)^(-2/3));
    Kv(309)=4.745e-11*exp(-272*(Tg*1.16e4)^(-1/3)+437*(Tg*1.16e4)^(-2/3));
    Kv(310)=4.745e-11*exp(-272*(Tg*1.16e4)^(-1/3)+437*(Tg*1.16e4)^(-2/3));
    Kv(311)=4.321e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
    Kv(312)=3.025e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
    Kv(313)=3.025e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
    Kv(314)=1.775e-17*exp(-108*(Tg*1.16e4)^(-1/3)+165*(Tg*1.16e4)^(-2/3));
    Kv(315)=2.13e-17*exp(-108*(Tg*1.16e4)^(-1/3)+165*(Tg*1.16e4)^(-2/3));
    Kv(316)=2.13e-17*exp(-108*(Tg*1.16e4)^(-1/3)+165*(Tg*1.16e4)^(-2/3));
    Kv(317)=2.384e-15*exp(-89*(Tg*1.16e4)^(-1/3)+234*(Tg*1.16e4)^(-2/3));
    Kv(318)=6.48e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3));
    Kv(319)=1.442e-14*exp(-88.2*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3));
    Kv(320)=2.628e-14*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3));
    Kv(321)=1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3));
    Kv(322)=2.927e-17*exp(21.4*(Tg*1.16e4)^(-1/3)-53*(Tg*1.16e4)^(-2/3));
    Kv(323)=4.825e-17*exp(20.1*(Tg*1.16e4)^(-1/3)-65*(Tg*1.16e4)^(-2/3));
    Kv(324)=1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3));
    Kv(325)=2.927e-17*exp(21.4*(Tg*1.16e4)^(-1/3)-53*(Tg*1.16e4)^(-2/3));
    Kv(326)=1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3));
    Kv(327:329)=1.0e-38;
    Kv(330:332)=1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3));
    Kv(357)=8.568e-7*exp(-404*(Tg*1.16e4)^(-1/3)+1096*(Tg*1.16e4)^(-2/3));
    Kv(358)=2.57e-7*exp(-404*(Tg*1.16e4)^(-1/3)+1096*(Tg*1.16e4)^(-2/3));
    Kv(359)=3.427e-7*exp(-404*(Tg*1.16e4)^(-1/3)+1096*(Tg*1.16e4)^(-2/3));
    Kv(360)=1.431e-11*exp(-252*(Tg*1.16e4)^(-1/3)+685*(Tg*1.16e4)^(-2/3));
    Kv(361)=4.293e-12*exp(-252*(Tg*1.16e4)^(-1/3)+685*(Tg*1.16e4)^(-2/3));
    Kv(362)=5.724e-12*exp(-252*(Tg*1.16e4)^(-1/3)+685*(Tg*1.16e4)^(-2/3));
    Kv(363)=7.143e-14*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
    Kv(364)=5e-14*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
    Kv(365)=5e-14*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
    Kv(366)=3.218e-15*exp(-137*(Tg*1.16e4)^(-1/3));
    Kv(367)=9.976e-15*exp(-137*(Tg*1.16e4)^(-1/3));
    Kv(368)=9.976e-15*exp(-137*(Tg*1.16e4)^(-1/3));
    Kv(369)=6.447e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
    Kv(370)=4.513e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
    Kv(371)=4.513e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
    Kv(372)=1.071e-15*exp(-137*(Tg*1.16e4)^(-1/3));
    Kv(373)=3.32e-15*exp(-137*(Tg*1.16e4)^(-1/3));
    Kv(374)=3.32e-15*exp(-137*(Tg*1.16e4)^(-1/3));
    Kv(375)=6.447e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
    Kv(376)=4.513e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
    Kv(377)=4.513e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
    Kv(378)=1.071e-15*exp(-137*(Tg*1.16e4)^(-1/3));
    Kv(379)=3.32e-15*exp(-137*(Tg*1.16e4)^(-1/3));
    Kv(380)=3.32e-15*exp(-137*(Tg*1.16e4)^(-1/3));
    Kv(381)=2.897e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
    Kv(382)=2.028e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
    Kv(383)=2.028e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
    Kv(384)=8.568e-7*exp(-406*(Tg*1.16e4)^(-1/3)+829*(Tg*1.16e4)^(-2/3));
    Kv(385)=2.57e-7*exp(-406*(Tg*1.16e4)^(-1/3)+829*(Tg*1.16e4)^(-2/3));
    Kv(386)=3.427e-7*exp(-406*(Tg*1.16e4)^(-1/3)+829*(Tg*1.16e4)^(-2/3));
    Kv(387)=1.725e-6*exp(-404*(Tg*1.16e4)^(-1/3)+1098*(Tg*1.16e4)^(-2/3));
    Kv(388)=5.175e-7*exp(-404*(Tg*1.16e4)^(-1/3)+1098*(Tg*1.16e4)^(-2/3));
    Kv(389)=6.9e-7*exp(-404*(Tg*1.16e4)^(-1/3)+1098*(Tg*1.16e4)^(-2/3));
    Kv(390)=2.882e-11*exp(-253*(Tg*1.16e4)^(-1/3)+683*(Tg*1.16e4)^(-2/3));
    Kv(391)=8.646e-12*exp(-253*(Tg*1.16e4)^(-1/3)+683*(Tg*1.16e4)^(-2/3));
    Kv(392)=1.153e-11*exp(-253*(Tg*1.16e4)^(-1/3)+683*(Tg*1.16e4)^(-2/3));
    Kv(393)=4.299e-11*exp(-241*(Tg*1.16e4)^(-1/3)+637*(Tg*1.16e4)^(-2/3));
    Kv(394)=4.299e-11*exp(-241*(Tg*1.16e4)^(-1/3)+635*(Tg*1.16e4)^(-2/3));
    Kv(395)=1.071e-15*exp(-87.7*(Tg*1.16e4)^(-1/3)+230*(Tg*1.16e4)^(-2/3));
    Kv(396)=2.157e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3));
    Kv(397)=4.344e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3));
    Kv(398)=6.48e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3));
    Kv(399)=9.667e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3));
    Kv(400)=9.667e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3));
    Kv(401)=4.332e-14*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3));
    Kv(402)=5.848e-14*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3));
    Kv(403)=2.157e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3));
    Kv(404)=2.378e-14*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3));
    Kv(405)=5.292e-14*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3));
    Kv(406)=1.066e-13*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3));
    Kv(407)=1.438e-13*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3));
    Kv(408)=5.305e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3));
    Kv(409)=5.305e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3));
    Kv(410:417)=Kv(87:94);
    Kv(418)=3.91e-16*exp(-46284/(Tg*11600));
    Kv(419)=3.91e-16*exp(-45514/(Tg*11600));
    Kv(420)=3.91e-16*exp(-44744/(Tg*11600));
    Kv(421:423)=Kv(418:420);
    Kv(424)=2.8e-17*exp(-24534/(Tg*11600));
    Kv(425)=2.8e-17*exp(-24052/(Tg*11600));
    Kv(426)=2.8e-17*exp(-23571/(Tg*11600));
    Kv(427:429)=1e-21;
    Kv(430:432)=9.0e-41;
    Kv(433:435)=1.5e-40;
    Kv(436:438)=3.1e-40;
    Kv(439:441)=9.0e-41;
    Kv(442:444)=1.0e-15;
    Kv(445:447)=5.5e-16;
    Kv(448:450)=9.4e-16;
    Kv(451:453)=4.5e-16;
    Kv(454:456)=1.1e-15;
    Kv(457:459)=0;
    Kv(460:462)=1.0e-41;
    Kv(463:465)=4.7e-41;
    Kv(466:468)=1.0e-41;
    Kv(469:471)=1.7e-42*(Tg*1.16e4)^(-1.2);
    Kv(472:474)=4.0e-18;
    Kv(475:477)=3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4));
    Kv(478:490)=2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4));
    Kv(523)=(5.075*Kv(59)+5.02*Kv(62)+4.711*Kv(65))*Fr(1,5)/Fr(1,1);
    Kv(524)=(5.075*Kv(60)+5.02*Kv(63)+4.711*Kv(66))*Fr(2,5)/Fr(2,1);
    Kv(525)=(5.075*Kv(61)+5.02*Kv(64)+4.711*Kv(67))*Fr(3,5)/Fr(3,1);
    Kv(526)=(6.114*Kv(59)+6.03*Kv(62)+5.572*Kv(65))*Fr(1,6)/Fr(1,1);
    Kv(527)=(6.114*Kv(60)+6.03*Kv(63)+5.572*Kv(66))*Fr(2,6)/Fr(2,1);
    Kv(528)=(6.114*Kv(61)+6.03*Kv(64)+5.572*Kv(67))*Fr(3,6)/Fr(3,1);
    Kv(529)=(7.16*Kv(59)+7.042*Kv(62)+6.409*Kv(65))*Fr(1,7)/Fr(1,1);
    Kv(530)=(7.16*Kv(60)+7.042*Kv(63)+6.409*Kv(66))*Fr(2,7)/Fr(2,1);
    Kv(531)=(7.16*Kv(61)+7.042*Kv(64)+6.409*Kv(67))*Fr(3,7)/Fr(3,1);
    Kv(532)=(8.214*Kv(59)+8.056*Kv(62)+7.223*Kv(65))*Fr(1,8)/Fr(1,1);
    Kv(533)=(8.214*Kv(60)+8.056*Kv(63)+7.223*Kv(66))*Fr(2,8)/Fr(2,1);
    Kv(534)=(8.214*Kv(61)+8.056*Kv(64)+7.223*Kv(67))*Fr(3,8)/Fr(3,1);
    Kv(535:542)=2.03e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3));
    Kv(543)=7.199e-17*exp(18.4*(Tg*1.16e4)^(-1/3)-76*(Tg*1.16e4)^(-2/3));
    Kv(544)=4.825e-17*exp(20.1*(Tg*1.16e4)^(-1/3)-65*(Tg*1.16e4)^(-2/3));
    Kv(545)=7.199e-17*exp(18.4*(Tg*1.16e4)^(-1/3)-76*(Tg*1.16e4)^(-2/3));
    Kv(546)=2.927e-17*exp(21.4*(Tg*1.16e4)^(-1/3)-53*(Tg*1.16e4)^(-2/3));
    Kv(547)=4.825e-17*exp(20.1*(Tg*1.16e4)^(-1/3)-65*(Tg*1.16e4)^(-2/3));
    Kv(548)=7.199e-17*exp(18.4*(Tg*1.16e4)^(-1/3)-76*(Tg*1.16e4)^(-2/3));
    Kv(549)=1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3));
    Kv(550)=2.927e-17*exp(21.4*(Tg*1.16e4)^(-1/3)-53*(Tg*1.16e4)^(-2/3));
    Kv(551)=4.825e-17*exp(20.1*(Tg*1.16e4)^(-1/3)-65*(Tg*1.16e4)^(-2/3));
    Kv(552)=7.199e-17*exp(18.4*(Tg*1.16e4)^(-1/3)-76*(Tg*1.16e4)^(-2/3));
    Kv(553)=1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3));
    Kv(554)=2.927e-17*exp(21.4*(Tg*1.16e4)^(-1/3)-53*(Tg*1.16e4)^(-2/3));
    Kv(555)=4.825e-17*exp(20.1*(Tg*1.16e4)^(-1/3)-65*(Tg*1.16e4)^(-2/3));
    Kv(556)=1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3));
    Kv(557)=2.927e-17*exp(21.4*(Tg*1.16e4)^(-1/3)-53*(Tg*1.16e4)^(-2/3));
    Kv(558)=1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3));
    Kv(559:576)=1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3));
    Kv(577)=3.91e-16*exp(-36160/(Tg*11600));
    Kv(578)=3.91e-16*exp(-33654/(Tg*11600));
    Kv(579)=3.91e-16*exp(-31148/(Tg*11600));
    Kv(580)=3.91e-16*exp(-28643/(Tg*11600));
    Kv(581:584)=Kv(577:580);
    Kv(585)=2.8e-17*exp(-18206/(Tg*11600));
    Kv(586)=2.8e-17*exp(-16640/(Tg*11600));
    Kv(587)=2.8e-17*exp(-15074/(Tg*11600));
    Kv(588)=2.8e-17*exp(-13508/(Tg*11600));
    Kv(589:592)=1e-21;
    Kv(593:596)=9.0e-41;
    Kv(597:600)=1.5e-40;
    Kv(601:604)=3.1e-40;
    Kv(605:608)=9.0e-41;
    Kv(609:612)=1.0e-38;
    Kv(613:616)=1.0e-15;
    Kv(617:620)=5.5e-16;
    Kv(621:624)=9.4e-16;
    Kv(625:628)=4.5e-16;
    Kv(629:632)=1.1e-15;
    Kv(633:636)=0;
    Kv(637:640)=1.0e-41;
    Kv(641:644)=4.7e-41;
    Kv(645:648)=1.0e-41;
    Kv(649:652)=1.7e-42*(Tg*1.16e4)^(-1.2);
    Kv(653:656)=4.0e-18;
    Kv(657:660)=3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4));
    Kv(661:664)=2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4));
    Kv(672)=1.43e-11*exp(-252*(Tg*11600)^(-1/3)+685*(Tg*11600)^(-2/3));
    Kv(673)=4.29e-12*exp(-252*(Tg*11600)^(-1/3)+685*(Tg*11600)^(-2/3));
    Kv(674)=5.72e-12*exp(-252*(Tg*11600)^(-1/3)+685*(Tg*11600)^(-2/3));
    Kv(675)=1.071e-15*exp(-137*(Tg*1.16e4)^(-1/3));
    Kv(676)=3.213e-15*exp(-137*(Tg*1.16e4)^(-1/3));
    Kv(677)=3.213e-15*exp(-137*(Tg*1.16e4)^(-1/3));
    Kv(678)=2.897e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
    Kv(679)=2.028e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
    Kv(680)=2.028e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
    Kv(681:682)=3.91e-16*exp(-26230/(Tg*11600));
    Kv(683)=2.8e-17*exp(-12000/(Tg*11600));
    Kv(684)=1.0e-21;
    Kv(685)=9.0e-41;
    Kv(686)=1.5e-40;
    Kv(687)=3.1e-40;
    Kv(688)=9.0e-41;
    Kv(689)=1.0e-38;
    Kv(690)=1.0e-15;
    Kv(691)=5.5e-16;
    Kv(692)=9.4e-16;
    Kv(693)=4.5e-16;
    Kv(694)=1.1e-15;
    Kv(695)=0;
    Kv(696)=1.0e-41;
    Kv(697)=4.7e-41;
    Kv(698)=1.0e-41;
    Kv(699)=1.7e-42*(Tg*1.16e4)^(-1.2);
    Kv(700)=4.0e-18;
    Kv(701)=3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4));
    Kv(702)=2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4));
    Kv(703)=5.305e-15*exp(-88*(Tg*11600)^(-1/3)+233*(Tg*11600)^(-2/3));
    Kv(704:721)=3.0e-40;
    Kv(722:739)=2.6e-38;
    Kv(740:757)=4.2e-38;
    Kv(758:775)=1.0e-18;
    Kv(776:793)=1.0e-20;
    Kv(794)=K(125)*0.916;
    Kv(795)=K(125)*0.093;
    Kv(796:811)=9.0e-16;
    Kv(812:829)=1.0e-25;
    Kv(830:847)=1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4));
    Kv(848:865)=8.3e-18;
    
    %K(1,number_of_reaction+1:number_of_reactions)=Kv(1,1:number_of_vreaction);
end

end

