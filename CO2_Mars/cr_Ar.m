function [nsp,species,X,R,lambda]=cr_model(Te,P,n,area_volume,Tg,only_e_balance)
%CR Model
%[nsp,species,X,]=cr_model(Te,P,n,area_volume,Tg,only_e_balance)

for i=1:30
    spec=['Ar^',num2str(i)];
    species{i}=spec;      charge(i)=0;        mass(i)=6.026e-26;
end
species{31}='e';         charge(31)=-1;       mass(31)=9.1e-31;        %kg
species{32}='Ar^+';      charge(32)=1;        mass(32)=6.026e-26;      %Argon ions
nsp=length(species);

if(Te==0)
    return;
end

ng=3.54e13*P*1e9;%m-3

%Rate coefficents
%Units = m^3.s^-1

number_of_reactions=219;
K(1:number_of_reactions)=0;
R(1:number_of_reactions)=0;
X(nsp,number_of_reactions+1)=0;

%First: Reactions involving generation/loss of electrons
%ionisation coefficents

K(219)=lookup('reaction-rate-cross_section-i1-jion-.lut',Te);    % e + Ar^1 --> Ar^+ +2e :electron impact ionisation from Ar^1 state to Ar^+
K(92)=lookup('reaction-rate-cross_section-i2-j100-.lut',Te);    % e + Ar^2 --> Ar^+ +2e :electron impact ionisation from Ar^2 state to Ar^+
K(93)=lookup('reaction-rate-cross_section-i3-j100-.lut',Te);    % e + Ar^3 --> Ar^+ +2e :electron impact ionisation from Ar^3 state to Ar^+
K(94)=lookup('reaction-rate-cross_section-i4-j100-.lut',Te);    % e + Ar^4 --> Ar^+ +2e :electron impact ionisation from Ar^4 state to Ar^+
K(95)=lookup('reaction-rate-cross_section-i5-j100-.lut',Te);    % e + Ar^5 --> Ar^+ +2e :electron impact ionisation from Ar^5 state to Ar^+
K(96)=lookup('reaction-rate-cross_section-i6-j100-.lut',Te);    % e + Ar^6 --> Ar^+ +2e :electron impact ionisation from Ar^6 state to Ar^+
K(97)=lookup('reaction-rate-cross_section-i7-j100-.lut',Te);    % e + Ar^7 --> Ar^+ +2e :electron impact ionisation from Ar^7 state to Ar^+
K(98)=lookup('reaction-rate-cross_section-i8-j100-.lut',Te);    % e + Ar^8 --> Ar^+ +2e :electron impact ionisation from Ar^8 state to Ar^+
K(99)=lookup('reaction-rate-cross_section-i9-j100-.lut',Te);    % e + Ar^9 --> Ar^+ +2e :electron impact ionisation from Ar^9 state to Ar^+
K(100)=lookup('reaction-rate-cross_section-i10-j100-.lut',Te);    % e + Ar^10 --> Ar^+ +2e :electron impact ionisation from Ar^10 state to Ar^+
K(101)=lookup('reaction-rate-cross_section-i11-j100-.lut',Te);    % e + Ar^11 --> Ar^+ +2e :electron impact ionisation from Ar^11 state to Ar^+
K(102)=lookup('reaction-rate-cross_section-i12-j100-.lut',Te);    % e + Ar^12 --> Ar^+ +2e :electron impact ionisation from Ar^12 state to Ar^+

%ionisation reactions from Ar12 and above (all use the same K!)
K(201)=K(102);   %Ar13--->Ar+ 
K(202)=K(102);   %Ar14--->Ar+ 
K(203)=K(102);   %Ar15--->Ar+ 
K(204)=K(102);   %Ar16--->Ar+ 
K(205)=K(102);   %Ar17--->Ar+ 
K(206)=K(102);   %Ar18--->Ar+ 
K(207)=K(102);   %Ar19--->Ar+ 
K(208)=K(102);   %Ar20--->Ar+ 
K(209)=K(102);   %Ar21--->Ar+
K(210)=K(102);   %Ar22--->Ar+
K(211)=K(102);   %Ar23--->Ar+ 
K(212)=K(102);   %Ar24--->Ar+ 
K(213)=K(102);   %Ar25--->Ar+ 
K(214)=K(102);   %Ar26--->Ar+ 
K(215)=K(102);   %Ar27--->Ar+ 
K(216)=K(102);   %Ar28--->Ar+
K(217)=K(102);   %Ar29--->Ar+ 
K(218)=K(102);   %Ar30--->Ar+ 

%ionisation reactions
R(219)=K(219)*n(31)*n(1); % e + Ar^1 --> Ar^+ +2e :electron impact ionisation from Ar^1 state to Ar^+
R(92)=K(92)*n(31)*n(2);    % e + Ar^2 --> Ar^+ +2e :electron impact ionisation from Ar^2 state to Ar^+
R(93)=K(93)*n(31)*n(3);    % e + Ar^3 --> Ar^+ +2e :electron impact ionisation from Ar^3 state to Ar^+
R(94)=K(94)*n(31)*n(4);    % e + Ar^4 --> Ar^+ +2e :electron impact ionisation from Ar^4 state to Ar^+
R(95)=K(95)*n(31)*n(5);    % e + Ar^5 --> Ar^+ +2e :electron impact ionisation from Ar^5 state to Ar^+
R(96)=K(96)*n(31)*n(6);    % e + Ar^6 --> Ar^+ +2e :electron impact ionisation from Ar^6 state to Ar^+
R(97)=K(97)*n(31)*n(7);    % e + Ar^7 --> Ar^+ +2e :electron impact ionisation from Ar^7 state to Ar^+
R(98)=K(98)*n(31)*n(8);    % e + Ar^8 --> Ar^+ +2e :electron impact ionisation from Ar^8 state to Ar^+
R(99)=K(99)*n(31)*n(9);    % e + Ar^9 --> Ar^+ +2e :electron impact ionisation from Ar^9 state to Ar^+
R(100)=K(100)*n(31)*n(10);    % e + Ar^10 --> Ar^+ +2e :electron impact ionisation from Ar^10 state to Ar^+
R(101)=K(101)*n(31)*n(11);    % e + Ar^11 --> Ar^+ +2e :electron impact ionisation from Ar^11 state to Ar^+
R(102)=K(102)*n(31)*n(12);    % e + Ar^12 --> Ar^+ +2e :electron impact ionisation from Ar^12 state to Ar^+

R(201)=K(201)*n(31)*n(13);    % e + Ar^13 --> Ar^+ +2e :electron impact ionisation from Ar^13 state to Ar^+
R(202)=K(202)*n(31)*n(14);    % e + Ar^14 --> Ar^+ +2e :electron impact ionisation from Ar^14 state to Ar^+
R(203)=K(203)*n(31)*n(15);    % e + Ar^15 --> Ar^+ +2e :electron impact ionisation from Ar^15 state to Ar^+
R(204)=K(204)*n(31)*n(16);    % e + Ar^16 --> Ar^+ +2e :electron impact ionisation from Ar^16 state to Ar^+
R(205)=K(205)*n(31)*n(17);    % e + Ar^17 --> Ar^+ +2e :electron impact ionisation from Ar^17 state to Ar^+
R(206)=K(206)*n(31)*n(18);    % e + Ar^18 --> Ar^+ +2e :electron impact ionisation from Ar^18 state to Ar^+
R(207)=K(207)*n(31)*n(19);    % e + Ar^19 --> Ar^+ +2e :electron impact ionisation from Ar^19 state to Ar^+
R(208)=K(208)*n(31)*n(20);    % e + Ar^20 --> Ar^+ +2e :electron impact ionisation from Ar^20 state to Ar^+
R(209)=K(209)*n(31)*n(21);    % e + Ar^21 --> Ar^+ +2e :electron impact ionisation from Ar^21 state to Ar^+
R(210)=K(210)*n(31)*n(22);    % e + Ar^22 --> Ar^+ +2e :electron impact ionisation from Ar^22 state to Ar^+
R(211)=K(211)*n(31)*n(23);    % e + Ar^23 --> Ar^+ +2e :electron impact ionisation from Ar^23 state to Ar^+
R(212)=K(212)*n(31)*n(24);    % e + Ar^24 --> Ar^+ +2e :electron impact ionisation from Ar^24 state to Ar^+
R(213)=K(213)*n(31)*n(25);    % e + Ar^25 --> Ar^+ +2e :electron impact ionisation from Ar^25 state to Ar^+
R(214)=K(214)*n(31)*n(26);    % e + Ar^26 --> Ar^+ +2e :electron impact ionisation from Ar^26 state to Ar^+
R(215)=K(215)*n(31)*n(27);    % e + Ar^27 --> Ar^+ +2e :electron impact ionisation from Ar^27 state to Ar^+
R(216)=K(216)*n(31)*n(28);    % e + Ar^28 --> Ar^+ +2e :electron impact ionisation from Ar^28 state to Ar^+
R(217)=K(217)*n(31)*n(29);    % e + Ar^29 --> Ar^+ +2e :electron impact ionisation from Ar^29 state to Ar^+
R(218)=K(218)*n(31)*n(30);    % e + Ar^30 --> Ar^+ +2e :electron impact ionisation from Ar^30 state to Ar^+

%Second: Rest of reactions
if(~only_e_balance)
    %Ground state coefficents
    K(1)=lookup('reaction-rate-cross_section-i1-j2-.lut',Te);      % e + Ar^1 --> Ar^2 +e :electron impact excitation from gnd state to 4s[3/2]2
    K(2)=lookup('reaction-rate-cross_section-i1-j3-.lut',Te);      % e + Ar^1 --> Ar^3 +e :electron impact excitation from gnd state to 4s[3/2]1
    K(3)=lookup('reaction-rate-cross_section-i1-j4-.lut',Te);      % e + Ar^1 --> Ar^4 +e :electron impact excitation from gnd state to 4s[1/2]0
    K(4)=lookup('reaction-rate-cross_section-i1-j5-.lut',Te);      % e + Ar^1 --> Ar^5 +e :electron impact excitation from gnd state to 4s[1/2]1
    K(5)=lookup('reaction-rate-cross_section-i1-j6-.lut',Te);      % e + Ar^1 --> Ar^6 +e :electron impact excitation from gnd state to 4p[1/2]1
    K(6)=lookup('reaction-rate-cross_section-i1-j7-.lut',Te);      % e + Ar^1 --> Ar^7 +e :electron impact excitation from gnd state to 4p[3/2]1,2 & 4p[5/2]2,3
    K(7)=lookup('reaction-rate-cross_section-i1-j8-.lut',Te);      % e + Ar^1 --> Ar^8 +e :electron impact excitation from gnd state to 4p[1/2]0
    K(8)=lookup('reaction-rate-cross_section-i1-j9-.lut',Te);      % e + Ar^1 --> Ar^9 +e :electron impact excitation from gnd state to 4p[3/1]1,2
    K(9)=lookup('reaction-rate-cross_section-i1-j10-.lut',Te);     % e + Ar^1 --> Ar^10 +e :electron impact excitation from gnd state to 4p[1/2]1
    K(10)=lookup('reaction-rate-cross_section-i1-j11-.lut',Te);    % e + Ar^1 --> Ar^11 +e :electron impact excitation from gnd state to 4p[1/2]0
    K(11)=lookup('reaction-rate-cross_section-i1-j12-.lut',Te);    % e + Ar^1 --> Ar^12 +e :electron impact excitation from gnd state to 3d[1/2]0,1 & 3d[3/2]2
    K(12)=lookup('reaction-rate-cross_section-i1-j13-.lut',Te);    % e + Ar^1 --> Ar^13 +e :electron impact excitation from gnd state to 3d[7/2]3,4
    K(13)=lookup('reaction-rate-cross_section-i1-j14-.lut',Te);    % e + Ar^1 --> Ar^14 +e :electron impact excitation from gnd state to 3d[3/2]1 & 3d[5/2]2,3 & 5s
    K(14)=lookup('reaction-rate-cross_section-i1-j15-.lut',Te);    % e + Ar^1 --> Ar^15 +e :electron impact excitation from gnd state to 3d[3/2]2 & 3d[5/2]2,3
    K(15)=lookup('reaction-rate-cross_section-i1-j16-.lut',Te);    % e + Ar^1 --> Ar^16 +e :electron impact excitation from gnd state to 5s
    K(16)=lookup('reaction-rate-cross_section-i1-j17-.lut',Te);    % e + Ar^1 --> Ar^17 +e :electron impact excitation from gnd state to 3d[3/2]1
    K(17)=lookup('reaction-rate-cross_section-i1-j18-.lut',Te);    % e + Ar^1 --> Ar^18 +e :electron impact excitation from gnd state to 5p
    K(18)=lookup('reaction-rate-cross_section-i1-j19-.lut',Te);    % e + Ar^1 --> Ar^19 +e :electron impact excitation from gnd state to 5p
    K(19)=lookup('reaction-rate-cross_section-i1-j20-.lut',Te);    % e + Ar^1 --> Ar^20 +e :electron impact excitation from gnd state to 4d.6s
    K(20)=lookup('reaction-rate-cross_section-i1-j21-.lut',Te);    % e + Ar^1 --> Ar^21 +e :electron impact excitation from gnd state to 4f
    K(21)=lookup('reaction-rate-cross_section-i1-j22-.lut',Te);    % e + Ar^1 --> Ar^22 +e :electron impact excitation from gnd state to 4d,6s
    K(22)=lookup('reaction-rate-cross_section-i1-j23-.lut',Te);    % e + Ar^1 --> Ar^23 +e :electron impact excitation from gnd state to 6p
    K(23)=lookup('reaction-rate-cross_section-i1-j24-.lut',Te);    % e + Ar^1 --> Ar^24 +e :electron impact excitation from gnd state to 4f
    K(24)=lookup('reaction-rate-cross_section-i1-j25-.lut',Te);    % e + Ar^1 --> Ar^25 +e :electron impact excitation from gnd state to 5d,7s
    K(25)=lookup('reaction-rate-cross_section-i1-j26-.lut',Te);    % e + Ar^1 --> Ar^26 +e :electron impact excitation from gnd state to 6p
    K(26)=lookup('reaction-rate-cross_section-i1-j27-.lut',Te);    % e + Ar^1 --> Ar^27 +e :electron impact excitation from gnd state to 5f,5g
    K(27)=lookup('reaction-rate-cross_section-i1-j28-.lut',Te);    % e + Ar^1 --> Ar^28 +e :electron impact excitation from gnd state to 7p
    K(28)=lookup('reaction-rate-cross_section-i1-j29-.lut',Te);    % e + Ar^1 --> Ar^29 +e :electron impact excitation from gnd state to 5d,7s
    K(29)=lookup('reaction-rate-cross_section-i1-j30-.lut',Te);    % e + Ar^1 --> Ar^30 +e :electron impact excitation from gnd state to 6d,8s

    %4s[3/2]2 coefficents
    K(30)=lookup('reaction-rate-cross_section-i2-j3-.lut',Te);    % e + Ar^2 --> Ar^3 +e :electron impact excitation from Ar^2 (4s[3/2]2) state to 4s[3/2]1
    K(31)=lookup('reaction-rate-cross_section-i2-j4-.lut',Te);    % e + Ar^2 --> Ar^4 +e :electron impact excitation from Ar^2 (4s[3/2]2) state to 4s[1/2]0
    K(32)=lookup('reaction-rate-cross_section-i2-j5-.lut',Te);    % e + Ar^2 --> Ar^5 +e :electron impact excitation from Ar^2 (4s[3/2]2) state to 4s[1/2]1
    K(33)=lookup('reaction-rate-cross_section-i2-j6-.lut',Te);    % e + Ar^2 --> Ar^6 +e :electron impact excitation from Ar^2 (4s[3/2]2) state to 4p[1/2]1
    K(34)=lookup('reaction-rate-cross_section-i2-j7-.lut',Te);    % e + Ar^2 --> Ar^7 +e :electron impact excitation from Ar^2 (4s[3/2]2) state to 4p[3/2]1,2 & 4p[5/2]2,3
    K(35)=lookup('reaction-rate-cross_section-i2-j8-.lut',Te);    % e + Ar^2 --> Ar^8 +e :electron impact excitation from Ar^2 (4s[3/2]2) state to 4p[1/2]0
    K(36)=lookup('reaction-rate-cross_section-i2-j9-.lut',Te);    % e + Ar^2 --> Ar^9 +e :electron impact excitation from Ar^2 (4s[3/2]2) state to 4p[3/1]1,2
    K(37)=lookup('reaction-rate-cross_section-i2-j10-.lut',Te);    % e + Ar^2 --> Ar^10 +e :electron impact excitation from Ar^2 (4s[3/2]2) state to 4p[1/2]1
    K(38)=lookup('reaction-rate-cross_section-i2-j11-.lut',Te);    % e + Ar^2 --> Ar^11 +e :electron impact excitation from Ar^2 (4s[3/2]2) state to 4p[1/2]0
    K(39)=lookup('reaction-rate-cross_section-i2-j12-.lut',Te);    % e + Ar^2 --> Ar^12 +e :electron impact excitation from Ar^2 (4s[3/2]2) state to 3d[1/2]0,1 & 3d[3/2]2
    K(40)=lookup('reaction-rate-cross_section-i2-j13-.lut',Te);    % e + Ar^2 --> Ar^13 +e :electron impact excitation from Ar^2 (4s[3/2]2) state to 3d[7/2]3,4
    K(41)=lookup('reaction-rate-cross_section-i2-j14-.lut',Te);    % e + Ar^2 --> Ar^14 +e :electron impact excitation from Ar^2 (4s[3/2]2) state to 3d[3/2]1 & 3d[5/2]2,3 & 5s
    K(42)=lookup('reaction-rate-cross_section-i2-j15-.lut',Te);    % e + Ar^2 --> Ar^15 +e :electron impact excitation from Ar^2 (4s[3/2]2) state to 3d[3/2]2 & 3d[5/2]2,3
    K(43)=lookup('reaction-rate-cross_section-i2-j16-.lut',Te);    % e + Ar^2 --> Ar^16 +e :electron impact excitation from Ar^2 (4s[3/2]2) state to 5s
    K(44)=lookup('reaction-rate-cross_section-i2-j17-.lut',Te);    % e + Ar^2 --> Ar^17 +e :electron impact excitation from Ar^2 (4s[3/2]2) state to 3d[3/2]1
    K(45)=lookup('reaction-rate-cross_section-i2-j18-.lut',Te);    % e + Ar^2 --> Ar^18 +e :electron impact excitation from Ar^2 (4s[3/2]2) state to 5p
    K(46)=lookup('reaction-rate-cross_section-i2-j19-.lut',Te);    % e + Ar^2 --> Ar^19 +e :electron impact excitation from Ar^2 (4s[3/2]2) state to 5p

    %4s[3/2]1 coefficents
    K(47)=lookup('reaction-rate-cross_section-i3-j4-.lut',Te);    % e + Ar^3 --> Ar^4 +e :electron impact excitation from Ar^3 (4s[3/2]1) state to 4s[1/2]0
    K(48)=lookup('reaction-rate-cross_section-i3-j5-.lut',Te);    % e + Ar^3 --> Ar^5 +e :electron impact excitation from Ar^3 (4s[3/2]1) state to 4s[1/2]1
    K(49)=lookup('reaction-rate-cross_section-i3-j6-.lut',Te);    % e + Ar^3 --> Ar^6 +e :electron impact excitation from Ar^3 (4s[3/2]1) state to 4p[1/2]1
    K(50)=lookup('reaction-rate-cross_section-i3-j7-.lut',Te);    % e + Ar^3 --> Ar^7 +e :electron impact excitation from Ar^3 (4s[3/2]1) state to 4p[3/2]1,2 & 4p[5/2]2,3
    K(51)=lookup('reaction-rate-cross_section-i3-j8-.lut',Te);    % e + Ar^3 --> Ar^8 +e :electron impact excitation from Ar^3 (4s[3/2]1) state to 4p[1/2]0
    K(52)=lookup('reaction-rate-cross_section-i3-j9-.lut',Te);    % e + Ar^3 --> Ar^9 +e :electron impact excitation from Ar^3 (4s[3/2]1) state to 4p[3/1]1,2
    K(53)=lookup('reaction-rate-cross_section-i3-j10-.lut',Te);    % e + Ar^3 --> Ar^10 +e :electron impact excitation from Ar^3 (4s[3/2]1) state to 4p[1/2]1
    K(54)=lookup('reaction-rate-cross_section-i3-j11-.lut',Te);    % e + Ar^3 --> Ar^11 +e :electron impact excitation from Ar^3 (4s[3/2]1) state to 4p[1/2]0
    K(55)=lookup('reaction-rate-cross_section-i3-j12-.lut',Te);    % e + Ar^3 --> Ar^12 +e :electron impact excitation from Ar^3 (4s[3/2]1) state to 3d[1/2]0,1 & 3d[3/2]2
    K(56)=lookup('reaction-rate-cross_section-i3-j13-.lut',Te);    % e + Ar^3 --> Ar^13 +e :electron impact excitation from Ar^3 (4s[3/2]1) state to 3d[7/2]3,4
    K(57)=lookup('reaction-rate-cross_section-i3-j14-.lut',Te);    % e + Ar^3 --> Ar^14 +e :electron impact excitation from Ar^3 (4s[3/2]1) state to 3d[3/2]1 & 3d[5/2]2,3 & 5s
    K(58)=lookup('reaction-rate-cross_section-i3-j15-.lut',Te);    % e + Ar^3 --> Ar^15 +e :electron impact excitation from Ar^3 (4s[3/2]1) state to 3d[3/2]2 & 3d[5/2]2,3
    K(59)=lookup('reaction-rate-cross_section-i3-j16-.lut',Te);    % e + Ar^3 --> Ar^16 +e :electron impact excitation from Ar^3 (4s[3/2]1) state to 5s
    K(60)=lookup('reaction-rate-cross_section-i3-j17-.lut',Te);    % e + Ar^3 --> Ar^17 +e :electron impact excitation from Ar^3 (4s[3/2]1) state to 3d[3/2]1
    K(61)=lookup('reaction-rate-cross_section-i3-j18-.lut',Te);    % e + Ar^3 --> Ar^18 +e :electron impact excitation from Ar^3 (4s[3/2]1) state to 5p
    K(62)=lookup('reaction-rate-cross_section-i3-j19-.lut',Te);    % e + Ar^3 --> Ar^19 +e :electron impact excitation from Ar^3 (4s[3/2]1) state to 5p

    %4s[1/2]0 coefficents
    K(63)=lookup('reaction-rate-cross_section-i4-j5-.lut',Te);    % e + Ar^4 --> Ar^5 +e :electron impact excitation from Ar^4 (4s[1/2]0) state to 4s[1/2]1
    K(64)=lookup('reaction-rate-cross_section-i4-j6-.lut',Te);    % e + Ar^4 --> Ar^6 +e :electron impact excitation from Ar^4 (4s[1/2]0) state to 4p[1/2]1
    K(65)=lookup('reaction-rate-cross_section-i4-j7-.lut',Te);    % e + Ar^4 --> Ar^7 +e :electron impact excitation from Ar^4 (4s[1/2]0) state to 4p[3/2]1,2 & 4p[5/2]2,3
    K(66)=lookup('reaction-rate-cross_section-i4-j8-.lut',Te);    % e + Ar^4 --> Ar^8 +e :electron impact excitation from Ar^4 (4s[1/2]0) state to 4p[1/2]0
    K(67)=lookup('reaction-rate-cross_section-i4-j9-.lut',Te);    % e + Ar^4 --> Ar^9 +e :electron impact excitation from Ar^4 (4s[1/2]0) state to 4p[3/1]1,2
    K(68)=lookup('reaction-rate-cross_section-i4-j10-.lut',Te);    % e + Ar^4 --> Ar^10 +e :electron impact excitation from Ar^4 (4s[1/2]0) state to 4p[1/2]1
    K(69)=lookup('reaction-rate-cross_section-i4-j11-.lut',Te);    % e + Ar^4 --> Ar^11 +e :electron impact excitation from Ar^4 (4s[1/2]0) state to 4p[1/2]0
    K(70)=lookup('reaction-rate-cross_section-i4-j12-.lut',Te);    % e + Ar^4 --> Ar^12 +e :electron impact excitation from Ar^4 (4s[1/2]0) state to 3d[1/2]0,1 & 3d[3/2]2
    K(71)=lookup('reaction-rate-cross_section-i4-j13-.lut',Te);    % e + Ar^4 --> Ar^13 +e :electron impact excitation from Ar^4 (4s[1/2]0) state to 3d[7/2]3,4
    K(72)=lookup('reaction-rate-cross_section-i4-j14-.lut',Te);    % e + Ar^4 --> Ar^14 +e :electron impact excitation from Ar^4 (4s[1/2]0) state to 3d[3/2]1 & 3d[5/2]2,3 & 5s
    K(73)=lookup('reaction-rate-cross_section-i4-j15-.lut',Te);    % e + Ar^4 --> Ar^15 +e :electron impact excitation from Ar^4 (4s[1/2]0) state to 3d[3/2]2 & 3d[5/2]2,3
    K(74)=lookup('reaction-rate-cross_section-i4-j16-.lut',Te);    % e + Ar^4 --> Ar^16 +e :electron impact excitation from Ar^4 (4s[1/2]0) state to 5s
    K(75)=lookup('reaction-rate-cross_section-i4-j17-.lut',Te);    % e + Ar^4 --> Ar^17 +e :electron impact excitation from Ar^4 (4s[1/2]0) state to 3d[3/2]1
    K(76)=lookup('reaction-rate-cross_section-i4-j18-.lut',Te);    % e + Ar^4 --> Ar^18 +e :electron impact excitation from Ar^4 (4s[1/2]0) state to 5p
    K(77)=lookup('reaction-rate-cross_section-i4-j19-.lut',Te);    % e + Ar^4 --> Ar^19 +e :electron impact excitation from Ar^4 (4s[1/2]0) state to 5p

    %4s[1/2]1 coefficents
    K(78)=lookup('reaction-rate-cross_section-i5-j6-.lut',Te);    % e + Ar^5 --> Ar^6 +e :electron impact excitation from Ar^5 (4s[1/2]1) state to 4p[1/2]1
    K(79)=lookup('reaction-rate-cross_section-i5-j7-.lut',Te);    % e + Ar^5 --> Ar^7 +e :electron impact excitation from Ar^5 (4s[1/2]1) state to 4p[3/2]1,2 & 4p[5/2]2,3
    K(80)=lookup('reaction-rate-cross_section-i5-j8-.lut',Te);    % e + Ar^5 --> Ar^8 +e :electron impact excitation from Ar^5 (4s[1/2]1) state to 4p[1/2]0
    K(81)=lookup('reaction-rate-cross_section-i5-j9-.lut',Te);    % e + Ar^5 --> Ar^9 +e :electron impact excitation from Ar^5 (4s[1/2]1) state to 4p[3/1]1,2
    K(82)=lookup('reaction-rate-cross_section-i5-j10-.lut',Te);    % e + Ar^5 --> Ar^10 +e :electron impact excitation from Ar^5 (4s[1/2]1) state to 4p[1/2]1
    K(83)=lookup('reaction-rate-cross_section-i5-j11-.lut',Te);    % e + Ar^5 --> Ar^11 +e :electron impact excitation from Ar^5 (4s[1/2]1) state to 4p[1/2]0
    K(84)=lookup('reaction-rate-cross_section-i5-j12-.lut',Te);    % e + Ar^5 --> Ar^12 +e :electron impact excitation from Ar^5 (4s[1/2]1) state to 3d[1/2]0,1 & 3d[3/2]2
    K(85)=lookup('reaction-rate-cross_section-i5-j13-.lut',Te);    % e + Ar^5 --> Ar^13 +e :electron impact excitation from Ar^5 (4s[1/2]1) state to 3d[7/2]3,4
    K(86)=lookup('reaction-rate-cross_section-i5-j14-.lut',Te);    % e + Ar^5 --> Ar^14 +e :electron impact excitation from Ar^5 (4s[1/2]1) state to 3d[3/2]1 & 3d[5/2]2,3 & 5s
    K(87)=lookup('reaction-rate-cross_section-i5-j15-.lut',Te);    % e + Ar^5 --> Ar^15 +e :electron impact excitation from Ar^5 (4s[1/2]1) state to 3d[3/2]2 & 3d[5/2]2,3
    K(88)=lookup('reaction-rate-cross_section-i5-j16-.lut',Te);    % e + Ar^5 --> Ar^16 +e :electron impact excitation from Ar^5 (4s[1/2]1) state to 5s
    K(89)=lookup('reaction-rate-cross_section-i5-j17-.lut',Te);    % e + Ar^5 --> Ar^17 +e :electron impact excitation from Ar^5 (4s[1/2]1) state to 3d[3/2]1
    K(90)=lookup('reaction-rate-cross_section-i5-j18-.lut',Te);    % e + Ar^5 --> Ar^18 +e :electron impact excitation from Ar^5 (4s[1/2]1) state to 5p
    K(91)=lookup('reaction-rate-cross_section-i5-j19-.lut',Te);    % e + Ar^5 --> Ar^19 +e :electron impact excitation from Ar^5 (4s[1/2]1) state to 5p

    %Einstein coefficients for photon emission
    
    %Transitions back to ground state:
    K(103)=5.6E+07;    %Ar3--->Ar1	+ hv (106.6nm)
    K(104)=5.10E+08;    %Ar5--->Ar1	+ hv (104.8nm)
    K(105)=2.91E+08;    %Ar14--->Ar1 + hv (87.9nm)
    K(106)=1.12E+08;    %Ar16--->Ar1 + hv (86.9nm)
    
    %Transitions to 4s[3/2]2:
    K(107)=1.89E+07;    %Ar6--->Ar2	+ hv (912nm)
    K(108)=1.8E+07;    %Ar7--->Ar2	+ hv (790nm)
    K(109)=2.21E+06;    %Ar9--->Ar2	+ hv (709nm)
    K(110)=6.39E+06;    %Ar10--->Ar2 + hv (696nm)
    K(111)=6.09E+05;    %Ar18--->Ar2 + hv (418nm)
    K(112)=4.55E+05;    %Ar19--->Ar2 + hv (394nm)

    %Transitions to 4s[3/2]1:
    K(113)=5.43E+06;    %Ar6--->Ar3	+ hv (966nm)
    K(114)=5.14E+07;    %Ar7--->Ar3	+ hv (831nm)
    K(115)=4.02E+07;    %Ar8--->Ar3	+ hv (751nm)
    K(116)=8.49E+06;    %Ar9--->Ar3	+ hv (741nm)
    K(117)=1.83E+06;    %Ar10--->Ar3 + hv (727nm)
    K(118)=2.36E+05;    %Ar11--->Ar3 + hv (668nm)
    K(119)=7.46E+05;    %Ar18--->Ar3 + hv (429nm)
    K(120)=1.33E+05;    %Ar19--->Ar3 + hv (404nm)   
    
    %Transitions to 4s[1/2]0:
    K(121)=9.80E+05;    %Ar6--->Ar4	+ hv (1047nm)
    K(122)=1.05E+06;    %Ar7--->Ar4	+ hv (890nm)
    K(123)=1.17E+07;    %Ar10--->Ar4 + hv (772nm)
    
    %Transitions to 4s[1/2]1:
    K(124)=1.90E+05;    %Ar6--->Ar5	+ hv (1149nm)
    K(125)=2.52E+06;    %Ar7--->Ar5	+ hv (963nm)
    K(126)=1.81E+07;    %Ar9--->Ar5 + hv (845nm)
    K(127)=1.53E+07;    %Ar10--->Ar5 + hv (826nm)
    K(128)=4.45E+07;    %Ar11--->Ar5 + hv (750nm)
    K(129)=3.83E+04;    %Ar18--->Ar5 + hv (462nm)
    K(130)=5.68E+05;    %Ar19--->Ar5 + hv (433nm)
    
    %Transitions From Ar12
    K(131)=6.8e6;     %Ar12--->Ar6	+ hv (1269nm)
    K(132)=7.29e5;    %Ar12--->Ar7	+ hv (1614nm)
    
    %Transitions From Ar13
    K(133)=5.54e6;    %Ar13--->Ar7	+ hv (1140nm)
    K(134)=3.67e5;    %Ar13--->Ar6	+ hv (1412nm)
    
    %Transitions From Ar14
    K(135)=4.31e6;    %Ar14--->Ar7	+ hv (1272nm)
    K(136)=2.8e6;     %Ar14--->Ar8	+ hv (1517nm)
    K(137)=1.33e5;    %Ar14--->Ar9	+ hv (1559nm)
    K(138)=3.3e5;     %Ar14--->Ar10	+ hv (1627nm) 
    
    %Transitions From Ar15
    K(139)=8.3e5;    %Ar15--->Ar7	+ hv (1113nm)
    K(140)=5.1e6;    %Ar15--->Ar9	+ hv (1327nm)
    
    %Transitions From Ar16
    K(141)=2.43e6;    %Ar16--->Ar6	+ hv (921nm)
    K(142)=9.86e5;    %Ar16--->Ar7	+ hv (1091nm)
    K(143)=3.8e5;     %Ar16--->Ar8	+ hv (1266nm)
    K(144)=6e6;       %Ar16--->Ar9	+ hv (1295nm)
    K(145)=2.13e6;    %Ar16--->Ar10	+ hv (1341nm)
    K(146)=1.9e6;     %Ar16--->Ar11	+ hv (1606nm)
    
    %Transitions From Ar17
    K(147)=3.9e5;     %Ar17--->Ar7	+ hv (1043nm)
    K(148)=2.36e6;    %Ar17--->Ar8	+ hv (1202nm)
    K(149)=3.4e5;     %Ar17--->Ar10	+ hv (1270nm)
    K(150)=5.2e6;     %Ar17--->Ar11	+ hv (1504nm)
    
    %Transition From Ar18 and Ar19 that are not accounted for previously:
    K(151)=4.85e4;    %Ar18--->Ar4	+ hv (445nm)
    K(152)=2.85e5;    %Ar18--->Ar5	+ hv (462nm)
    K(153)=5.5e5;     %Ar19--->Ar4	+ hv (417nm)
    K(154)=1.31e6;    %Ar19--->Ar5	+ hv (433nm)
    
    %Transition From Ar20
    K(155)=2.24e6;    %Ar20--->Ar6	+ hv (657nm)
    K(156)=6.71e5;    %Ar20--->Ar7	+ hv (739nm)
    K(157)=9.2e5;     %Ar20--->Ar8	+ hv (816nm)
    K(158)=3.21e5;    %Ar20--->Ar9	+ hv (828nm)
    K(159)=6.17e5;    %Ar20--->Ar10	+ hv (846nm)
     
    %Transition From Ar21
    K(160)=2.635e6;    %Ar21--->Ar12 + hv (1213nm)
    K(161)=6.5e6;      %Ar21--->Ar13 + hv (1359nm)
    K(162)=2.11e6;     %Ar21--->Ar14 + hv (1519nm)
    
    %Transition From Ar22
    K(163)=7.46e5;    %Ar22--->Ar6	+ hv (599nm)
    K(164)=1.94e5;    %Ar22--->Ar7	+ hv (666nm)
    K(165)=1.04e5;    %Ar22--->Ar8	+ hv (728nm)
    K(166)=8.75e5;    %Ar22--->Ar9	+ hv (737nm)
    K(167)=8.17e5;    %Ar22--->Ar10 + hv (752nm)
    K(168)=3.59e5;    %Ar22--->Ar11 + hv (828nm)
    
    %Transition From Ar23
    K(169)=1.1e5;    %Ar23--->Ar2	+ hv (356nm)
    K(170)=2.39e5;   %Ar23--->Ar3	+ hv (364nm)
    K(171)=7e4;      %Ar23--->Ar4	+ hv (375nm)
    K(172)=4.03e5;   %Ar23--->Ar5	+ hv (387nm)
    
    %Transition From Ar24
    K(173)=1.58e6;   %Ar24--->Ar12	+ hv (1034nm)
    K(174)=3.46e5;   %Ar24--->Ar13	+ hv (1138nm)
    K(175)=5.3e6;    %Ar24--->Ar16	+ hv (1491nm)
    K(176)=7.7e6;    %Ar24--->Ar17	+ hv (1591nm)
    
    %Transition From Ar25
    K(177)=1.496e6;   %Ar25--->Ar6	+ hv (552nm)
    K(178)=4.28e5;    %Ar25--->Ar7	+ hv (608nm)
    K(179)=1.61e5;    %Ar25--->Ar8	+ hv (659nm)
    K(180)=1.46e5;    %Ar25--->Ar9	+ hv (667nm)
    K(181)=2.42e5;    %Ar25--->Ar10	+ hv (679nm)

    %Transition From Ar26
    K(182)=6.7e4;   %Ar26--->Ar3	+ hv (346nm)
    K(183)=1.2e5;   %Ar26--->Ar4	+ hv (356nm)
    K(184)=2.93e5;  %Ar26--->Ar5	+ hv (367nm)
    
    %Transition From Ar27
    K(185)=5.276e6;  %Ar27--->Ar12	+ hv (931nm)
    K(186)=6.5e6;    %Ar27--->Ar13	+ hv (1015nm)
    K(187)=1.35e5;   %Ar27--->Ar14	+ hv (1102nm)
    K(188)=3.5e6;    %Ar27--->Ar15	+ hv (1257nm)
    K(189)=1.6e5;    %Ar27--->Ar17	+ hv (1360nm)
    
    %Transition From Ar28
    K(190)=4.5e5;  %Ar28--->Ar25	+ hv (9611nm)?
    
    %Transition From Ar29
    K(191)=5.37e5;  %Ar29--->Ar6	+ hv (512nm)
    K(192)=4.08e5;  %Ar29--->Ar7	+ hv (561nm)
    K(193)=9.4e4;   %Ar29--->Ar8	+ hv (604nm)
    K(194)=4.53e5;  %Ar29--->Ar9	+ hv (611nm)
    K(195)=4.6e5;   %Ar29--->Ar10	+ hv (621nm)
 
    %Transition From Ar30
    K(196)=1.202e6; %Ar30--->Ar6	+ hv (508nm)
    K(197)=3.22e5;  %Ar30--->Ar7	+ hv (555nm)
    K(198)=2.13e5;  %Ar30--->Ar8	+ hv (597nm)
    K(199)=7.8e4;   %Ar30--->Ar9	+ hv (604nm)
    K(200)=8.5e4;   %Ar30--->Ar10	+ hv (614nm)
        
    
    %Ground state reactions
    R(1)=K(1)*n(31)*n(1);      % e + Ar^1 --> Ar^2 +e :electron impact excitation from gnd state to 4s[3/2]2
    R(2)=K(2)*n(31)*n(1);      % e + Ar^1 --> Ar^3 +e :electron impact excitation from gnd state to 4s[3/2]1
    R(3)=K(3)*n(31)*n(1);      % e + Ar^1 --> Ar^4 +e :electron impact excitation from gnd state to 4s[1/2]0
    R(4)=K(4)*n(31)*n(1);      % e + Ar^1 --> Ar^5 +e :electron impact excitation from gnd state to 4s[1/2]1
    R(5)=K(5)*n(31)*n(1);      % e + Ar^1 --> Ar^6 +e :electron impact excitation from gnd state to 4p[1/2]1
    R(6)=K(6)*n(31)*n(1);      % e + Ar^1 --> Ar^7 +e :electron impact excitation from gnd state to 4p[3/2]1,2 & 4p[5/2]2,3
    R(7)=K(7)*n(31)*n(1);      % e + Ar^1 --> Ar^8 +e :electron impact excitation from gnd state to 4p[1/2]0
    R(8)=K(8)*n(31)*n(1);      % e + Ar^1 --> Ar^9 +e :electron impact excitation from gnd state to 4p[3/1]1,2
    R(9)=K(9)*n(31)*n(1);      % e + Ar^1 --> Ar^10 +e :electron impact excitation from gnd state to 4p[1/2]1
    R(10)=K(10)*n(31)*n(1);    % e + Ar^1 --> Ar^11 +e :electron impact excitation from gnd state to 4p[1/2]0
    R(11)=K(11)*n(31)*n(1);    % e + Ar^1 --> Ar^12 +e :electron impact excitation from gnd state to 3d[1/2]0,1 & 3d[3/2]2
    R(12)=K(12)*n(31)*n(1);    % e + Ar^1 --> Ar^13 +e :electron impact excitation from gnd state to 3d[7/2]3,4
    R(13)=K(13)*n(31)*n(1);    % e + Ar^1 --> Ar^14 +e :electron impact excitation from gnd state to 3d[3/2]1 & 3d[5/2]2,3 & 5s
    R(14)=K(14)*n(31)*n(1);    % e + Ar^1 --> Ar^15 +e :electron impact excitation from gnd state to 3d[3/2]2 & 3d[5/2]2,3
    R(15)=K(15)*n(31)*n(1);    % e + Ar^1 --> Ar^16 +e :electron impact excitation from gnd state to 5s
    R(16)=K(16)*n(31)*n(1);    % e + Ar^1 --> Ar^17 +e :electron impact excitation from gnd state to 3d[3/2]1
    R(17)=K(17)*n(31)*n(1);    % e + Ar^1 --> Ar^18 +e :electron impact excitation from gnd state to 5p
    R(18)=K(18)*n(31)*n(1);    % e + Ar^1 --> Ar^19 +e :electron impact excitation from gnd state to 5p
    R(19)=K(19)*n(31)*n(1);    % e + Ar^1 --> Ar^20 +e :electron impact excitation from gnd state to 4d.6s
    R(20)=K(20)*n(31)*n(1);    % e + Ar^1 --> Ar^21 +e :electron impact excitation from gnd state to 4f
    R(21)=K(21)*n(31)*n(1);    % e + Ar^1 --> Ar^22 +e :electron impact excitation from gnd state to 4d,6s
    R(22)=K(22)*n(31)*n(1);    % e + Ar^1 --> Ar^23 +e :electron impact excitation from gnd state to 6p
    R(23)=K(23)*n(31)*n(1);    % e + Ar^1 --> Ar^24 +e :electron impact excitation from gnd state to 4f
    R(24)=K(24)*n(31)*n(1);    % e + Ar^1 --> Ar^25 +e :electron impact excitation from gnd state to 5d,7s
    R(25)=K(25)*n(31)*n(1);    % e + Ar^1 --> Ar^26 +e :electron impact excitation from gnd state to 6p
    R(26)=K(26)*n(31)*n(1);    % e + Ar^1 --> Ar^27 +e :electron impact excitation from gnd state to 5f,5g
    R(27)=K(27)*n(31)*n(1);    % e + Ar^1 --> Ar^28 +e :electron impact excitation from gnd state to 7p
    R(28)=K(28)*n(31)*n(1);    % e + Ar^1 --> Ar^29 +e :electron impact excitation from gnd state to 5d,7s
    R(29)=K(29)*n(31)*n(1);    % e + Ar^1 --> Ar^30 +e :electron impact excitation from gnd state to 6d,8s

    %4s[3/2]2 reactions
    R(30)=K(30)*n(31)*n(2);    % e + Ar^2 --> Ar^3 +e :electron impact excitation from Ar^2 (4s[3/2]2) state to 4s[3/2]1
    R(31)=K(31)*n(31)*n(2);    % e + Ar^2 --> Ar^4 +e :electron impact excitation from Ar^2 (4s[3/2]2) state to 4s[1/2]0
    R(32)=K(32)*n(31)*n(2);    % e + Ar^2 --> Ar^5 +e :electron impact excitation from Ar^2 (4s[3/2]2) state to 4s[1/2]1
    R(33)=K(33)*n(31)*n(2);    % e + Ar^2 --> Ar^6 +e :electron impact excitation from Ar^2 (4s[3/2]2) state to 4p[1/2]1
    R(34)=K(34)*n(31)*n(2);    % e + Ar^2 --> Ar^7 +e :electron impact excitation from Ar^2 (4s[3/2]2) state to 4p[3/2]1,2 & 4p[5/2]2,3
    R(35)=K(35)*n(31)*n(2);    % e + Ar^2 --> Ar^8 +e :electron impact excitation from Ar^2 (4s[3/2]2) state to 4p[1/2]0
    R(36)=K(36)*n(31)*n(2);    % e + Ar^2 --> Ar^9 +e :electron impact excitation from Ar^2 (4s[3/2]2) state to 4p[3/1]1,2
    R(37)=K(37)*n(31)*n(2);    % e + Ar^2 --> Ar^10 +e :electron impact excitation from Ar^2 (4s[3/2]2) state to 4p[1/2]1
    R(38)=K(38)*n(31)*n(2);    % e + Ar^2 --> Ar^11 +e :electron impact excitation from Ar^2 (4s[3/2]2) state to 4p[1/2]0
    R(39)=K(39)*n(31)*n(2);    % e + Ar^2 --> Ar^12 +e :electron impact excitation from Ar^2 (4s[3/2]2) state to 3d[1/2]0,1 & 3d[3/2]2
    R(40)=K(40)*n(31)*n(2);    % e + Ar^2 --> Ar^13 +e :electron impact excitation from Ar^2 (4s[3/2]2) state to 3d[7/2]3,4
    R(41)=K(41)*n(31)*n(2);    % e + Ar^2 --> Ar^14 +e :electron impact excitation from Ar^2 (4s[3/2]2) state to 3d[3/2]1 & 3d[5/2]2,3 & 5s
    R(42)=K(42)*n(31)*n(2);    % e + Ar^2 --> Ar^15 +e :electron impact excitation from Ar^2 (4s[3/2]2) state to 3d[3/2]2 & 3d[5/2]2,3
    R(43)=K(43)*n(31)*n(2);    % e + Ar^2 --> Ar^16 +e :electron impact excitation from Ar^2 (4s[3/2]2) state to 5s
    R(44)=K(44)*n(31)*n(2);    % e + Ar^2 --> Ar^17 +e :electron impact excitation from Ar^2 (4s[3/2]2) state to 3d[3/2]1
    R(45)=K(45)*n(31)*n(2);    % e + Ar^2 --> Ar^18 +e :electron impact excitation from Ar^2 (4s[3/2]2) state to 5p
    R(46)=K(46)*n(31)*n(2);    % e + Ar^2 --> Ar^19 +e :electron impact excitation from Ar^2 (4s[3/2]2) state to 5p

    %4s[3/2]1 reactions
    R(47)=K(47)*n(31)*n(3);    % e + Ar^3 --> Ar^4 +e :electron impact excitation from Ar^3 (4s[3/2]1) state to 4s[1/2]0
    R(48)=K(48)*n(31)*n(3);    % e + Ar^3 --> Ar^5 +e :electron impact excitation from Ar^3 (4s[3/2]1) state to 4s[1/2]1
    R(49)=K(49)*n(31)*n(3);    % e + Ar^3 --> Ar^6 +e :electron impact excitation from Ar^3 (4s[3/2]1) state to 4p[1/2]1
    R(50)=K(50)*n(31)*n(3);    % e + Ar^3 --> Ar^7 +e :electron impact excitation from Ar^3 (4s[3/2]1) state to 4p[3/2]1,2 & 4p[5/2]2,3
    R(51)=K(51)*n(31)*n(3);    % e + Ar^3 --> Ar^8 +e :electron impact excitation from Ar^3 (4s[3/2]1) state to 4p[1/2]0
    R(52)=K(52)*n(31)*n(3);    % e + Ar^3 --> Ar^9 +e :electron impact excitation from Ar^3 (4s[3/2]1) state to 4p[3/1]1,2
    R(53)=K(53)*n(31)*n(3);    % e + Ar^3 --> Ar^10 +e :electron impact excitation from Ar^3 (4s[3/2]1) state to 4p[1/2]1
    R(54)=K(54)*n(31)*n(3);    % e + Ar^3 --> Ar^11 +e :electron impact excitation from Ar^3 (4s[3/2]1) state to 4p[1/2]0
    R(55)=K(55)*n(31)*n(3);    % e + Ar^3 --> Ar^12 +e :electron impact excitation from Ar^3 (4s[3/2]1) state to 3d[1/2]0,1 & 3d[3/2]2
    R(56)=K(56)*n(31)*n(3);    % e + Ar^3 --> Ar^13 +e :electron impact excitation from Ar^3 (4s[3/2]1) state to 3d[7/2]3,4
    R(57)=K(57)*n(31)*n(3);    % e + Ar^3 --> Ar^14 +e :electron impact excitation from Ar^3 (4s[3/2]1) state to 3d[3/2]1 & 3d[5/2]2,3 & 5s
    R(58)=K(58)*n(31)*n(3);    % e + Ar^3 --> Ar^15 +e :electron impact excitation from Ar^3 (4s[3/2]1) state to 3d[3/2]2 & 3d[5/2]2,3
    R(59)=K(59)*n(31)*n(3);    % e + Ar^3 --> Ar^16 +e :electron impact excitation from Ar^3 (4s[3/2]1) state to 5s
    R(60)=K(60)*n(31)*n(3);    % e + Ar^3 --> Ar^17 +e :electron impact excitation from Ar^3 (4s[3/2]1) state to 3d[3/2]1
    R(61)=K(61)*n(31)*n(3);    % e + Ar^3 --> Ar^18 +e :electron impact excitation from Ar^3 (4s[3/2]1) state to 5p
    R(62)=K(62)*n(31)*n(3);    % e + Ar^3 --> Ar^19 +e :electron impact excitation from Ar^3 (4s[3/2]1) state to 5p

    %4s[1/2]0 reactions
    R(63)=K(63)*n(31)*n(4);    % e + Ar^4 --> Ar^5 +e :electron impact excitation from Ar^4 (4s[1/2]0) state to 4s[1/2]1
    R(64)=K(64)*n(31)*n(4);    % e + Ar^4 --> Ar^6 +e :electron impact excitation from Ar^4 (4s[1/2]0) state to 4p[1/2]1
    R(65)=K(65)*n(31)*n(4);    % e + Ar^4 --> Ar^7 +e :electron impact excitation from Ar^4 (4s[1/2]0) state to 4p[3/2]1,2 & 4p[5/2]2,3
    R(66)=K(66)*n(31)*n(4);    % e + Ar^4 --> Ar^8 +e :electron impact excitation from Ar^4 (4s[1/2]0) state to 4p[1/2]0
    R(67)=K(67)*n(31)*n(4);    % e + Ar^4 --> Ar^9 +e :electron impact excitation from Ar^4 (4s[1/2]0) state to 4p[3/1]1,2
    R(68)=K(68)*n(31)*n(4);    % e + Ar^4 --> Ar^10 +e :electron impact excitation from Ar^4 (4s[1/2]0) state to 4p[1/2]1
    R(69)=K(69)*n(31)*n(4);    % e + Ar^4 --> Ar^11 +e :electron impact excitation from Ar^4 (4s[1/2]0) state to 4p[1/2]0
    R(70)=K(70)*n(31)*n(4);    % e + Ar^4 --> Ar^12 +e :electron impact excitation from Ar^4 (4s[1/2]0) state to 3d[1/2]0,1 & 3d[3/2]2
    R(71)=K(71)*n(31)*n(4);    % e + Ar^4 --> Ar^13 +e :electron impact excitation from Ar^4 (4s[1/2]0) state to 3d[7/2]3,4
    R(72)=K(72)*n(31)*n(4);    % e + Ar^4 --> Ar^14 +e :electron impact excitation from Ar^4 (4s[1/2]0) state to 3d[3/2]1 & 3d[5/2]2,3 & 5s
    R(73)=K(73)*n(31)*n(4);    % e + Ar^4 --> Ar^15 +e :electron impact excitation from Ar^4 (4s[1/2]0) state to 3d[3/2]2 & 3d[5/2]2,3
    R(74)=K(74)*n(31)*n(4);    % e + Ar^4 --> Ar^16 +e :electron impact excitation from Ar^4 (4s[1/2]0) state to 5s
    R(75)=K(75)*n(31)*n(4);    % e + Ar^4 --> Ar^17 +e :electron impact excitation from Ar^4 (4s[1/2]0) state to 3d[3/2]1
    R(76)=K(76)*n(31)*n(4);    % e + Ar^4 --> Ar^18 +e :electron impact excitation from Ar^4 (4s[1/2]0) state to 5p
    R(77)=K(77)*n(31)*n(4);    % e + Ar^4 --> Ar^19 +e :electron impact excitation from Ar^4 (4s[1/2]0) state to 5p

    %4s[1/2]1 reactions
    R(78)=K(78)*n(31)*n(5);    % e + Ar^5 --> Ar^6 +e :electron impact excitation from Ar^5 (4s[1/2]1) state to 4p[1/2]1
    R(79)=K(79)*n(31)*n(5);    % e + Ar^5 --> Ar^7 +e :electron impact excitation from Ar^5 (4s[1/2]1) state to 4p[3/2]1,2 & 4p[5/2]2,3
    R(80)=K(80)*n(31)*n(5);    % e + Ar^5 --> Ar^8 +e :electron impact excitation from Ar^5 (4s[1/2]1) state to 4p[1/2]0
    R(81)=K(81)*n(31)*n(5);    % e + Ar^5 --> Ar^9 +e :electron impact excitation from Ar^5 (4s[1/2]1) state to 4p[3/1]1,2
    R(82)=K(82)*n(31)*n(5);    % e + Ar^5 --> Ar^10 +e :electron impact excitation from Ar^5 (4s[1/2]1) state to 4p[1/2]1
    R(83)=K(83)*n(31)*n(5);    % e + Ar^5 --> Ar^11 +e :electron impact excitation from Ar^5 (4s[1/2]1) state to 4p[1/2]0
    R(84)=K(84)*n(31)*n(5);    % e + Ar^5 --> Ar^12 +e :electron impact excitation from Ar^5 (4s[1/2]1) state to 3d[1/2]0,1 & 3d[3/2]2
    R(85)=K(85)*n(31)*n(5);    % e + Ar^5 --> Ar^13 +e :electron impact excitation from Ar^5 (4s[1/2]1) state to 3d[7/2]3,4
    R(86)=K(86)*n(31)*n(5);    % e + Ar^5 --> Ar^14 +e :electron impact excitation from Ar^5 (4s[1/2]1) state to 3d[3/2]1 & 3d[5/2]2,3 & 5s
    R(87)=K(87)*n(31)*n(5);    % e + Ar^5 --> Ar^15 +e :electron impact excitation from Ar^5 (4s[1/2]1) state to 3d[3/2]2 & 3d[5/2]2,3
    R(88)=K(88)*n(31)*n(5);    % e + Ar^5 --> Ar^16 +e :electron impact excitation from Ar^5 (4s[1/2]1) state to 5s
    R(89)=K(89)*n(31)*n(5);    % e + Ar^5 --> Ar^17 +e :electron impact excitation from Ar^5 (4s[1/2]1) state to 3d[3/2]1
    R(90)=K(90)*n(31)*n(5);    % e + Ar^5 --> Ar^18 +e :electron impact excitation from Ar^5 (4s[1/2]1) state to 5p
    R(91)=K(91)*n(31)*n(5);    % e + Ar^5 --> Ar^19 +e :electron impact excitation from Ar^5 (4s[1/2]1) state to 5p
    
    %Spontaneous emission reactions
      %Return to ground state
      R(103)=K(103)*n(3);   %Ar3--->Ar1	+ hv    
      R(104)=K(104)*n(5);   %Ar5--->Ar1	+ hv    
      R(105)=K(105)*n(14);  %Ar14--->Ar1 + hv   
      R(106)=K(106)*n(16);  %Ar16--->Ar1 + hv   
      
      %Return to 4s[3/1]2
      R(107)=K(107)*n(6);   %Ar6--->Ar2	+ hv   
      R(108)=K(108)*n(7);   %Ar7--->Ar2	+ hv   
      R(109)=K(109)*n(9);   %Ar9--->Ar2	+ hv   
      R(110)=K(110)*n(10);  %Ar10--->Ar2 + hv  
      R(111)=K(111)*n(18);  %Ar18--->Ar2 + hv  
      R(112)=K(112)*n(19);  %Ar19--->Ar2 + hv  
      
      %Return to 4s[3/1]1
      R(113)=K(113)*n(6);   %Ar6--->Ar3	+ hv   
      R(114)=K(114)*n(7);   %Ar7--->Ar3	+ hv   
      R(115)=K(115)*n(8);   %Ar8--->Ar3	+ hv   
      R(116)=K(116)*n(9);   %Ar9--->Ar3	+ hv   
      R(117)=K(117)*n(10);  %Ar10--->Ar3 + hv  
      R(118)=K(118)*n(11);  %Ar11--->Ar3 + hv  
      R(119)=K(119)*n(18);  %Ar18--->Ar3 + hv  
      R(120)=K(120)*n(19);  %Ar19--->Ar3 + hv  
      
      %Return to 4s[1/2]0
      R(121)=K(121)*n(6);   %Ar6--->Ar4	+ hv
      R(122)=K(122)*n(7);   %Ar7--->Ar4	+ hv
      R(123)=K(123)*n(10);   %Ar10--->Ar4 + hv
      
      %Return to 4s[1/2]1
      R(124)=K(124)*n(6);    %Ar6--->Ar5 + hv
      R(125)=K(125)*n(7);    %Ar7--->Ar5 + hv
      R(126)=K(126)*n(9);    %Ar9--->Ar5 + hv
      R(127)=K(127)*n(10);   %Ar10--->Ar5 + hv
      R(128)=K(128)*n(11);   %Ar11--->Ar5 + hv
      R(129)=K(129)*n(18);   %Ar18--->Ar5 + hv
      R(130)=K(130)*n(19);   %Ar19--->Ar5 + hv
      
      %Transitions from Ar12
      R(131)=K(131)*n(12);   %Ar12--->Ar6 + hv
      R(132)=K(132)*n(12);   %Ar12--->Ar7 + hv
      
      %Transitions from Ar13
      R(133)=K(133)*n(13);   %Ar13--->Ar7 + hv
      R(134)=K(134)*n(13);   %Ar13--->Ar6 + hv
      
      %Transitions from Ar14
      R(135)=K(135)*n(14);   %Ar14--->Ar7 + hv
      R(136)=K(136)*n(14);   %Ar14--->Ar8 + hv
      R(137)=K(137)*n(14);   %Ar14--->Ar9 + hv
      R(138)=K(138)*n(14);   %Ar14--->Ar10 + hv
      
      %Transitions from Ar15
      R(139)=K(139)*n(15);   %Ar15--->Ar7 + hv
      R(140)=K(140)*n(15);   %Ar15--->Ar9 + hv
      
      %Transitions from Ar16
      R(141)=K(141)*n(16);   %Ar16--->Ar6 + hv
      R(142)=K(142)*n(16);   %Ar16--->Ar7 + hv
      R(143)=K(143)*n(16);   %Ar16--->Ar8 + hv
      R(144)=K(144)*n(16);   %Ar16--->Ar9 + hv
      R(145)=K(145)*n(16);   %Ar16--->Ar10 + hv
      R(146)=K(146)*n(16);   %Ar16--->Ar11 + hv
      
      %Transitions from Ar17
      R(147)=K(147)*n(17);   %Ar17--->Ar7 + hv
      R(148)=K(148)*n(17);   %Ar17--->Ar8 + hv
      R(149)=K(149)*n(17);   %Ar17--->Ar10 + hv
      R(150)=K(150)*n(17);   %Ar17--->Ar11 + hv
      
      %Transitions from Ar18 & Ar19
      R(151)=K(151)*n(18);   %Ar18--->Ar4 + hv
      R(152)=K(152)*n(18);   %Ar18--->Ar5 + hv
      R(153)=K(153)*n(19);   %Ar19--->Ar4 + hv
      R(154)=K(154)*n(19);   %Ar19--->Ar5 + hv
      
      %Transitions from Ar20
      R(155)=K(155)*n(20);   %Ar20--->Ar6 + hv
      R(156)=K(156)*n(20);   %Ar20--->Ar7 + hv
      R(157)=K(157)*n(20);   %Ar20--->Ar8 + hv
      R(158)=K(158)*n(20);   %Ar20--->Ar9 + hv
      R(159)=K(159)*n(20);   %Ar20--->Ar10 + hv
      
      %Transitions from Ar21
      R(160)=K(160)*n(21);   %Ar21--->Ar12 + hv
      R(161)=K(161)*n(21);   %Ar21--->Ar13 + hv
      R(162)=K(162)*n(21);   %Ar21--->Ar14 + hv
      
      %Transitions from Ar22
      R(163)=K(163)*n(22);   %Ar22--->Ar6 + hv
      R(164)=K(164)*n(22);   %Ar22--->Ar7 + hv
      R(165)=K(165)*n(22);   %Ar22--->Ar8 + hv
      R(166)=K(166)*n(22);   %Ar22--->Ar9 + hv
      R(167)=K(167)*n(22);   %Ar22--->Ar10 + hv
      R(168)=K(168)*n(22);   %Ar22--->Ar11 + hv
      
      %Transitions from Ar23
      R(169)=K(169)*n(23);   %Ar23--->Ar2 + hv
      R(170)=K(170)*n(23);   %Ar23--->Ar3 + hv
      R(171)=K(171)*n(23);   %Ar23--->Ar4 + hv
      R(172)=K(172)*n(23);   %Ar23--->Ar5 + hv
      
      %Transitions from Ar24
      R(173)=K(173)*n(24);   %Ar24--->Ar2 + hv
      R(174)=K(174)*n(24);   %Ar24--->Ar3 + hv
      R(175)=K(175)*n(24);   %Ar24--->Ar4 + hv
      R(176)=K(176)*n(24);   %Ar24--->Ar5 + hv
      
      %Transitions from Ar25
      R(177)=K(177)*n(25);   %Ar25--->Ar5 + hv
      R(178)=K(178)*n(25);   %Ar25--->Ar7 + hv
      R(179)=K(179)*n(25);   %Ar25--->Ar8 + hv
      R(180)=K(180)*n(25);   %Ar25--->Ar9 + hv
      R(181)=K(181)*n(25);   %Ar25--->Ar10 + hv
      
      %Transitions from Ar26
      R(182)=K(182)*n(26);   %Ar26--->Ar3 + hv
      R(183)=K(183)*n(26);   %Ar26--->Ar4 + hv
      R(184)=K(184)*n(26);   %Ar26--->Ar5 + hv
      
      %Transitions from Ar27
      R(185)=K(185)*n(27);   %Ar27--->Ar12 + hv
      R(186)=K(186)*n(27);   %Ar27--->Ar13 + hv
      R(187)=K(187)*n(27);   %Ar27--->Ar14 + hv
      R(188)=K(188)*n(27);   %Ar27--->Ar15 + hv      
      R(189)=K(189)*n(27);   %Ar27--->Ar17 + hv
      
      %Transitions from Ar28
      R(190)=K(190)*n(28);   %Ar28--->Ar25 + hv
      
      %Transitions from Ar29
      R(191)=K(191)*n(29);   %Ar29--->Ar6 + hv
      R(192)=K(192)*n(29);   %Ar29--->Ar7 + hv
      R(193)=K(193)*n(29);   %Ar29--->Ar8 + hv
      R(194)=K(194)*n(29);   %Ar29--->Ar9 + hv
      R(195)=K(195)*n(29);   %Ar29--->Ar10 + hv
      
      
      %Transitions from Ar30
      R(196)=K(191)*n(30);   %Ar30--->Ar6 + hv
      R(197)=K(192)*n(30);   %Ar30--->Ar7 + hv
      R(198)=K(193)*n(30);   %Ar30--->Ar8 + hv
      R(199)=K(194)*n(30);   %Ar30--->Ar9 + hv
      R(200)=K(195)*n(30);   %Ar30--->Ar10 + hv
      
end
%Wall losses
%Excited states de-excite on the walls
for loop=2:30 
    W(loop)=-area_volume*0.25*n(loop)*sqrt(8*1.6e-19*Tg/pi/mass(loop));%metastables
end
W(1)=sum(W(2:30));
W(32)=-area_volume*n(32)*sqrt(1.6e-19*Te/mass(32));%ions
W(31)=W(32);%electron loss = ions loss


% electrons, species 31
isp=31;
%Producing reactions
X(isp,219)=R(219);    % e + Ar^1 --> Ar^+ +2e :electron impact ionisation
X(isp,92)=R(92);    % e + Ar^2 --> Ar^+ +2e :electron impact ionisation from Ar^2 state to Ar^+
X(isp,93)=R(93);    % e + Ar^3 --> Ar^+ +2e :electron impact ionisation from Ar^3 state to Ar^+
X(isp,94)=R(94);    % e + Ar^4 --> Ar^+ +2e :electron impact ionisation from Ar^4 state to Ar^+
X(isp,95)=R(95);    % e + Ar^5 --> Ar^+ +2e :electron impact ionisation from Ar^5 state to Ar^+
X(isp,96)=R(96);    % e + Ar^6 --> Ar^+ +2e :electron impact ionisation from Ar^6 state to Ar^+
X(isp,97)=R(97);    % e + Ar^7 --> Ar^+ +2e :electron impact ionisation from Ar^7 state to Ar^+
X(isp,98)=R(98);    % e + Ar^8 --> Ar^+ +2e :electron impact ionisation from Ar^8 state to Ar^+
X(isp,99)=R(99);    % e + Ar^9 --> Ar^+ +2e :electron impact ionisation from Ar^9 state to Ar^+
X(isp,100)=R(100);    % e + Ar^10 --> Ar^+ +2e :electron impact ionisation from Ar^10 state to Ar^+
X(isp,101)=R(101);    % e + Ar^11 --> Ar^+ +2e :electron impact ionisation from Ar^11 state to Ar^+
X(isp,102)=R(102);    % e + Ar^12 --> Ar^+ +2e :electron impact ionisation from Ar^12 state to Ar^+
X(isp,201)=R(201);    % e + Ar^13 --> Ar^+ +2e :electron impact ionisation 
X(isp,202)=R(202);    % e + Ar^14 --> Ar^+ +2e :electron impact ionisation 
X(isp,203)=R(203);    % e + Ar^15 --> Ar^+ +2e :electron impact ionisation 
X(isp,204)=R(204);    % e + Ar^16 --> Ar^+ +2e :electron impact ionisation 
X(isp,205)=R(205);    % e + Ar^17 --> Ar^+ +2e :electron impact ionisation 
X(isp,206)=R(206);    % e + Ar^18 --> Ar^+ +2e :electron impact ionisation 
X(isp,207)=R(207);    % e + Ar^19 --> Ar^+ +2e :electron impact ionisation 
X(isp,208)=R(208);    % e + Ar^20 --> Ar^+ +2e :electron impact ionisation 
X(isp,209)=R(209);    % e + Ar^21 --> Ar^+ +2e :electron impact ionisation 
X(isp,210)=R(210);    % e + Ar^22 --> Ar^+ +2e :electron impact ionisation 
X(isp,211)=R(211);    % e + Ar^23 --> Ar^+ +2e :electron impact ionisation 
X(isp,212)=R(212);    % e + Ar^24 --> Ar^+ +2e :electron impact ionisation 
X(isp,213)=R(213);    % e + Ar^25 --> Ar^+ +2e :electron impact ionisation 
X(isp,214)=R(214);    % e + Ar^26 --> Ar^+ +2e :electron impact ionisation 
X(isp,215)=R(215);    % e + Ar^27 --> Ar^+ +2e :electron impact ionisation 
X(isp,216)=R(216);    % e + Ar^28 --> Ar^+ +2e :electron impact ionisation 
X(isp,217)=R(217);    % e + Ar^29 --> Ar^+ +2e :electron impact ionisation 
X(isp,218)=R(218);    % e + Ar^30 --> Ar^+ +2e :electron impact ionisation 

%Destroying reactions

if(~only_e_balance)
% Ar^1: ground state argon, species 1
isp=1;
%Producing reactions
X(isp,103:106)=R(103:106);
%Destroying reactions
X(isp,1:29)=-R(1:29);
X(isp,219)=-R(219);

% Ar^2: 4s[3/2]2, species 2
isp=2;
%Producing reactions
X(isp,1)=R(1);
%spontaneous emission
X(isp,107:112)=R(107:112);
X(isp,169)=R(169);
%Destroying reactions
X(isp,30:46)=-R(30:46);
X(isp,92)=-R(92);

% Ar^3, species 3
isp=3;
%Producing reactions
X(isp,2)=R(2);
X(isp,30)=R(30);
%spontaneous emission
X(isp,113:120)=R(113:120);
X(isp,170)=R(170);
X(isp,182)=R(182);
X(isp,174)=R(174);
%Destroying reactions
X(isp,47:62)=-R(47:62);
X(isp,93)=-R(93);
X(isp,103)=-R(103);

% Ar^4, species 4
isp=4;
%Producing reactions
X(isp,3)=R(3);
X(isp,31)=R(31);
X(isp,47)=R(47);
%spontaneous emission
X(isp,121:123)=R(121:123);
X(isp,151)=R(151);
X(isp,153)=R(153);
X(isp,171)=R(171);
X(isp,183)=R(183);
%Destroying reactions
X(isp,63:77)=-R(63:77);
X(isp,94)=-R(94);

% Ar^5, species 5
isp=5;
%Producing reactions
X(isp,4)=R(4);
X(isp,32)=R(32);
X(isp,48)=R(48);
X(isp,63)=R(63);
%spontaneous emission
X(isp,124:130)=R(124:130);
X(isp,152)=R(152);
X(isp,154)=R(154);
X(isp,172)=R(172);
X(isp,176)=R(176);
X(isp,177)=R(177);
X(isp,184)=R(184);
%Destroying reactions
X(isp,78:91)=-R(78:91);
X(isp,95)=-R(95);
X(isp,104)=-R(104);

% Ar^6, species 6
isp=6;
%Producing reactions
X(isp,5)=R(5);
X(isp,33)=R(33);
X(isp,49)=R(49);
X(isp,64)=R(64);
X(isp,78)=R(78);
X(isp,131)=R(131);
X(isp,134)=R(134);
X(isp,141)=R(141);
X(isp,155)=R(155);
X(isp,163)=R(163);
X(isp,191)=R(191);
X(isp,196)=R(196);
%Destroying reactions
X(isp,96)=-R(96);
X(isp,107)=-R(107);
X(isp,113)=-R(113);
X(isp,121)=-R(121);
X(isp,124)=-R(124);

% Ar^7, species 7
isp=7;
%Producing reactions
X(isp,6)=R(6);
X(isp,34)=R(34);
X(isp,50)=R(50);
X(isp,65)=R(65);
X(isp,79)=R(79);
X(isp,132)=R(132);
X(isp,133)=R(133);
X(isp,135)=R(135);
X(isp,139)=R(139);
X(isp,142)=R(142);
X(isp,147)=R(147);
X(isp,156)=R(156);
X(isp,164)=R(164);
X(isp,178)=R(178);
X(isp,192)=R(192);
X(isp,197)=R(197);
%Destroying reactions
X(isp,97)=-R(97);
X(isp,108)=-R(108);
X(isp,114)=-R(114);
X(isp,122)=-R(122);
X(isp,125)=-R(125);

% Ar^8, species 8
isp=8;
%Producing reactions
X(isp,7)=R(7);
X(isp,35)=R(35);
X(isp,51)=R(51);
X(isp,66)=R(66);
X(isp,80)=R(80);
X(isp,136)=R(136);
X(isp,143)=R(143);
X(isp,148)=R(148);
X(isp,157)=R(157);
X(isp,165)=R(165);
X(isp,179)=R(179);
X(isp,193)=R(193);
X(isp,198)=R(198);
%Destroying reactions
X(isp,98)=-R(98);
X(isp,115)=-R(115);

% Ar^9, species 9
isp=9;
%Producing reactions
X(isp,8)=R(8);
X(isp,36)=R(36);
X(isp,52)=R(52);
X(isp,67)=R(67);
X(isp,81)=R(81);
X(isp,137)=R(137);
X(isp,140)=R(140);
X(isp,144)=R(144);
X(isp,158)=R(158);
X(isp,166)=R(166);
X(isp,180)=R(180);
X(isp,194)=R(194);
X(isp,199)=R(199);
%Destroying reactions
X(isp,99)=-R(99);
X(isp,109)=-R(109);
X(isp,116)=-R(116);
X(isp,126)=-R(126);


% Ar^10, species 10
isp=10;
%Producing reactions
X(isp,9)=R(9);
X(isp,37)=R(37);
X(isp,53)=R(53);
X(isp,68)=R(68);
X(isp,82)=R(82);
X(isp,138)=R(138);
X(isp,145)=R(145);
X(isp,149)=R(149);
X(isp,159)=R(159);
X(isp,167)=R(167);
X(isp,181)=R(181);
X(isp,195)=R(195);
X(isp,200)=R(200);
%Destroying reactions
X(isp,100)=-R(100);
X(isp,110)=-R(110);
X(isp,117)=-R(117);
X(isp,123)=-R(123);
X(isp,127)=-R(127);

% Ar^11, species 11
isp=11;
%Producing reactions
X(isp,10)=R(10);
X(isp,38)=R(38);
X(isp,54)=R(54);
X(isp,69)=R(69);
X(isp,83)=R(83);
X(isp,146)=R(146);
X(isp,150)=R(150);
X(isp,168)=R(168);
%Destroying reactions
X(isp,101)=-R(101);
X(isp,118)=-R(118);
X(isp,120)=-R(120);
X(isp,128)=-R(128);

% Ar^12, species 12
isp=12;
%Producing reactions
X(isp,11)=R(11);
X(isp,39)=R(39);
X(isp,55)=R(55);
X(isp,70)=R(70);
X(isp,84)=R(84);
X(isp,160)=R(160);
X(isp,173)=R(173);
X(isp,185)=R(185);
%Destroying reactions
X(isp,102)=-R(102);
X(isp,131)=-R(131);
X(isp,132)=-R(132);

% Ar^13, species 13
isp=13;
%Producing reactions
X(isp,12)=R(12);
X(isp,40)=R(40);
X(isp,56)=R(56);
X(isp,71)=R(71);
X(isp,85)=R(85);
X(isp,161)=R(161);
X(isp,174)=R(174);
X(isp,186)=R(186);
%Destroying reactions
X(isp,133)=-R(133);
X(isp,134)=-R(134);
X(isp,201)=-R(201);

% Ar^14, species 14
isp=14;
%Producing reactions
X(isp,13)=R(13);
X(isp,41)=R(41);
X(isp,57)=R(57);
X(isp,72)=R(72);
X(isp,86)=R(86);
X(isp,162)=R(162);
X(isp,187)=R(187);
%Destroying reactions
X(isp,105)=-R(105);
X(isp,135)=-R(135);
X(isp,136)=-R(136);
X(isp,137)=-R(137);
X(isp,138)=-R(138);
X(isp,202)=-R(202);


% Ar^15, species 15
isp=15;
%Producing reactions
X(isp,14)=R(14);
X(isp,42)=R(42);
X(isp,58)=R(58);
X(isp,73)=R(73);
X(isp,87)=R(87);
X(isp,188)=R(188);
%Destroying reactions
X(isp,139:140)=-R(139:140);
X(isp,203)=-R(203);

% Ar^16, species 16
isp=16;
%Producing reactions
X(isp,15)=R(15);
X(isp,43)=R(43);
X(isp,59)=R(59);
X(isp,74)=R(74);
X(isp,88)=R(88);
X(isp,175)=R(175);
%Destroying reactions
X(isp,106)=-R(106);
X(isp,141:146)=-R(141:146);
X(isp,204)=-R(204);

% Ar^17, species 17
isp=17;
%Producing reactions
X(isp,16)=R(16);
X(isp,44)=R(44);
X(isp,60)=R(60);
X(isp,75)=R(75);
X(isp,89)=R(89);
X(isp,176)=R(176);
X(isp,189)=R(189);
%Destroying reactions
X(isp,147:150)=-R(147:150);
X(isp,205)=-R(205);


% Ar^18, species 18
isp=18;
%Producing reactions
X(isp,17)=R(17);
X(isp,45)=R(45);
X(isp,61)=R(61);
X(isp,76)=R(76);
X(isp,90)=R(90);
%Destroying reactions
X(isp,111)=-R(111);
X(isp,119)=-R(119);
X(isp,129)=-R(129);
X(isp,151)=-R(151);
X(isp,152)=-R(152);
X(isp,206)=-R(206);

% Ar^19, species 19
isp=19;
%Producing reactions
X(isp,18)=R(18);
X(isp,46)=R(46);
X(isp,62)=R(62);
X(isp,77)=R(77);
X(isp,91)=R(91);
%Destroying reactions
X(isp,112)=-R(112);
X(isp,120)=-R(120);
X(isp,130)=-R(130);
X(isp,153)=-R(153);
X(isp,154)=-R(154);
X(isp,207)=-R(207);

% Ar^20, species 20
isp=20;
%Producing reactions
X(isp,19)=R(19);
%Destroying reactions
X(isp,155:159)=-R(155:159);
X(isp,208)=-R(208);

% Ar^21, species 21
isp=21;
%Producing reactions
X(isp,20)=R(20);
%Destroying reactions
X(isp,160:162)=-R(160:162);
X(isp,209)=-R(209);


% Ar^22, species 22
isp=22;
%Producing reactions
X(isp,21)=R(21);
%Destroying reactions
X(isp,163:168)=-R(163:168);
X(isp,210)=-R(210);

% Ar^23, species 23
isp=23;
%Producing reactions
X(isp,22)=R(22);
%Destroying reactions
X(isp,169:172)=-R(169:172);
X(isp,211)=-R(211);

% Ar^24, species 24
isp=24;
%Producing reactions
X(isp,23)=R(23);
%Destroying reactions
X(isp,173:176)=-R(173:176);
X(isp,212)=-R(212);

% Ar^25, species 25
isp=25;
%Producing reactions
X(isp,24)=R(24);
X(isp,190)=R(190);
%Destroying reactions
X(isp,177:181)=-R(177:181);
X(isp,213)=-R(213);

% Ar^26, species 26
isp=26;
%Producing reactions
X(isp,25)=R(25);
%Destroying reactions
X(isp,182:184)=-R(182:184);
X(isp,214)=-R(214);

% Ar^27, species 27
isp=27;
%Producing reactions
X(isp,26)=R(26);
%Destroying reactions
X(isp,185:189)=-R(185:189);
X(isp,215)=-R(215);

% Ar^28, species 28
isp=28;
%Producing reactions
X(isp,27)=R(27);
%Destroying reactions
X(isp,190)=-R(190);
X(isp,216)=-R(216);

% Ar^29, species 29
isp=29;
%Producing reactions
X(isp,28)=R(28);
%Destroying reactions
X(isp,217)=-R(217);
X(isp,191:195)=-R(191:195);


% Ar^30, species 30
isp=30;
%Producing reactions
X(isp,29)=R(29);
%Destroying reactions
X(isp,218)=-R(218);
X(isp,196:200)=-R(196:200);

% Argon ions, species 32
isp=32;
%Producing reactions
X(isp,219)=R(219);  % e + Ar^1 --> Ar^+ +2e :electron impact ionisation from Ar^1 state to Ar^+
X(isp,92)=R(92);    % e + Ar^2 --> Ar^+ +2e :electron impact ionisation from Ar^2 state to Ar^+
X(isp,93)=R(93);    % e + Ar^3 --> Ar^+ +2e :electron impact ionisation from Ar^3 state to Ar^+
X(isp,94)=R(94);    % e + Ar^4 --> Ar^+ +2e :electron impact ionisation from Ar^4 state to Ar^+
X(isp,95)=R(95);    % e + Ar^5 --> Ar^+ +2e :electron impact ionisation from Ar^5 state to Ar^+
X(isp,96)=R(96);    % e + Ar^6 --> Ar^+ +2e :electron impact ionisation from Ar^6 state to Ar^+
X(isp,97)=R(97);    % e + Ar^7 --> Ar^+ +2e :electron impact ionisation from Ar^7 state to Ar^+
X(isp,98)=R(98);    % e + Ar^8 --> Ar^+ +2e :electron impact ionisation from Ar^8 state to Ar^+
X(isp,99)=R(99);    % e + Ar^9 --> Ar^+ +2e :electron impact ionisation from Ar^9 state to Ar^+
X(isp,100)=R(100);    % e + Ar^10 --> Ar^+ +2e :electron impact ionisation from Ar^10 state to Ar^+
X(isp,101)=R(101);    % e + Ar^11 --> Ar^+ +2e :electron impact ionisation from Ar^11 state to Ar^+
X(isp,102)=R(102);    % e + Ar^12 --> Ar^+ +2e :electron impact ionisation from Ar^12 state to Ar^+
X(isp,201)=R(201);    % e + Ar^13 --> Ar^+ +2e :electron impact ionisation 
X(isp,202)=R(202);    % e + Ar^14 --> Ar^+ +2e :electron impact ionisation 
X(isp,203)=R(203);    % e + Ar^15 --> Ar^+ +2e :electron impact ionisation 
X(isp,204)=R(204);    % e + Ar^16 --> Ar^+ +2e :electron impact ionisation 
X(isp,205)=R(205);    % e + Ar^17 --> Ar^+ +2e :electron impact ionisation 
X(isp,206)=R(206);    % e + Ar^18 --> Ar^+ +2e :electron impact ionisation 
X(isp,207)=R(207);    % e + Ar^19 --> Ar^+ +2e :electron impact ionisation 
X(isp,208)=R(208);    % e + Ar^20 --> Ar^+ +2e :electron impact ionisation 
X(isp,209)=R(209);    % e + Ar^21 --> Ar^+ +2e :electron impact ionisation 
X(isp,210)=R(210);    % e + Ar^22 --> Ar^+ +2e :electron impact ionisation 
X(isp,211)=R(211);    % e + Ar^23 --> Ar^+ +2e :electron impact ionisation 
X(isp,212)=R(212);    % e + Ar^24 --> Ar^+ +2e :electron impact ionisation 
X(isp,213)=R(213);    % e + Ar^25 --> Ar^+ +2e :electron impact ionisation 
X(isp,214)=R(214);    % e + Ar^26 --> Ar^+ +2e :electron impact ionisation 
X(isp,215)=R(215);    % e + Ar^27 --> Ar^+ +2e :electron impact ionisation 
X(isp,216)=R(216);    % e + Ar^28 --> Ar^+ +2e :electron impact ionisation 
X(isp,217)=R(217);    % e + Ar^29 --> Ar^+ +2e :electron impact ionisation 
X(isp,218)=R(218);    % e + Ar^30 --> Ar^+ +2e :electron impact ionisation     
%Destroying reactions
end;

%Wall losses
X(:,(length(K)+1))=W';

%Map reactions to wavelengths
lambda=zeros(1,length(R));
lambda(103)=	106.6;
lambda(104)=	104.8;
lambda(105)=	87.9;
lambda(106)=	86.9;
lambda(107)=	912;
lambda(108)=	790;
lambda(109)=	709;
lambda(110)=	696;
lambda(111)=	418;
lambda(112)=	394;
lambda(113)=	966;
lambda(114)=	831;
lambda(115)=	751;
lambda(116)=	741;
lambda(117)=	727;
lambda(118)=	668;
lambda(119)=	429;
lambda(120)=	404;
lambda(121)=	1047;
lambda(122)=	890;
lambda(123)=	772;
lambda(124)=	1149;
lambda(125)=	963;
lambda(126)=	845;
lambda(127)=	826;
lambda(128)=	750;
lambda(129)=	462;
lambda(130)=	433;
lambda(131)=	1269;
lambda(132)=	1614;
lambda(133)=	1140;
lambda(134)=	1412;
lambda(135)=	1272;
lambda(136)=	1517;
lambda(137)=	1559;
lambda(138)=	1627;
lambda(139)=	1113;
lambda(140)=	1327;
lambda(141)=	921;
lambda(142)=	1091;
lambda(143)=	1266;
lambda(144)=	1295;
lambda(145)=	1341;
lambda(146)=	1606;
lambda(147)=	1043;
lambda(148)=	1202;
lambda(149)=	1270;
lambda(150)=	1504;
lambda(151)=	445;
lambda(152)=	462;
lambda(153)=	417;
lambda(154)=	433;
lambda(155)=	657;
lambda(156)=	739;
lambda(157)=	816;
lambda(158)=	828;
lambda(159)=	846;
lambda(160)=	1213;
lambda(161)=	1359;
lambda(162)=	1519;
lambda(163)=	599;
lambda(164)=	666;
lambda(165)=	728;
lambda(166)=	737;
lambda(167)=	752;
lambda(168)=	828;
lambda(169)=	356;
lambda(170)=	364;
lambda(171)=	375;
lambda(172)=	387;
lambda(173)=	1034;
lambda(174)=	1138;
lambda(175)=	1491;
lambda(176)=	1591;
lambda(177)=	552;
lambda(178)=	608;
lambda(179)=	659;
lambda(180)=	667;
lambda(181)=	679;
lambda(182)=	346;
lambda(183)=	356;
lambda(184)=	367;
lambda(185)=	931;
lambda(186)=	1015;
lambda(187)=	1102;
lambda(188)=	1257;
lambda(189)=	1360;
lambda(190)=	9611;
lambda(191)=	512;
lambda(192)=	561;
lambda(193)=	604;
lambda(194)=	611;
lambda(195)=	621;
lambda(196)=	508;
lambda(197)=	555;
lambda(198)=	597;
lambda(199)=	604;
lambda(200)=	614;
