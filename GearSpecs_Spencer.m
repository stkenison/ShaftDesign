clc; clear all; close all

%Material Properties, 1045 HR Steel
s_yield = 45000; s_ult = 82000; 

%constants for analysis
H = 9.5; %input horsepower
w_p = 3000; %input shaft rotation speed
w_g = 400; %shaft rotation speeds
phi_n = 20; %normal pressure angle
v_helix = 30; %helix angle
R_gear = 0.9; %reliability of gears 
Lifetime = 1000; %lifetime of gears
k_gear = 1; %gear constants, full-depth teeth 
R_bearing = 0.975; %combined reliability of gears
af_bearing = 1.4; % bearing constants
SF_bend = 1.2; %required AGMA safety factor against bending
SH_wear = 1.2; %%required AGMA safety factor against wear

%gear calculations
m = w_g/w_p %module, eqn. 18-2
phi_t = atand((tand(phi_n))/(cosd(v_helix))) %transverse pressure angle, eqn. 13-19
N_p = 2*k_gear*cosd(v_helix)/((1+2*m)*(sind(phi_t))^2)*(m+sqrt(m^2+(1+2*m)*(sind(phi_t))^2)) %minimum number of teeth on pinion, 13-22
N_p = 6 %actual number of teeth on pinion
N_g = (N_p^2*(sind(phi_t))^2-4*k_gear^2*(cosd(v_helix))^2)/(4*k_gear*cosd(v_helix)-2*N_p*(sind(phi_t))^2) %minimum number of teeth on gear, eqn. 13-23
N_g = N_p/m %actual number of teeth on gear

%torque calculations
T_p = H/w_p*550/(2*pi)*60 %torque on pinion shaft in lbf*ft %eqn. 18-1
T_g = H/w_g*550/(2*pi)*60 %torque on gear shaft in lbf*ft %eqn. 18-1

%gear diameter calculations
syms P
P_min = vpasolve(15==N_p/P+N_g/P+2/P+2*1.5+2*1.5,P) %Minimum diametral pitch, eqn.18-3
P_n = 6 %actual diametral pitch, taken from table 13-2 
d_p = N_p/P_n %actual diameter of pinion, eqn. 13-1
d_g = N_g/P_n %actual diameter of gear, eqn. 13-1

%velocities & loads
V = pi*d_g*w_g/12 %pitch line velocity, 13-34
w_t = 33000*H/V %tangential force, 13-35
w_r = w_t*tand(phi_t) %radial force, 13-40
w_a = w_t*tand(v_helix) %axial force, 13-40

%face width is typically 3-5 times circular pitch
p_n = pi/P_n %circular pitch, 13-4
F = p_n*5 %recommended face width assumption using general standard of x3-5, %13-4
F = 2 %actual face width

%gear wear & bending
m_g = N_g/N_p %gear ratio
a = 1/P_n %addendum, Table 13-4
r_p = d_p/2 %pinion radius, eqn. 13-6
r_g = d_g/2 %pinion gear, eqn. 13-6
rb_p = r_p*cosd(phi_t) %base radius, eqn. 13-6
rb_g = r_g*cosd(phi_t) %gear radius, eqn. 13-6
Z = sqrt((r_p+a)^2-rb_p^2)+sqrt((r_g+a)^2-rb_g^2)-(r_p+r_g)*sind(phi_t) %stress-cycle factor for pitting resistance, eqn. 14-25
p_N = p_n*cosd(phi_n) %normal circular pitch, eqn. 14-24

%trying to determine if m_n EQN 14-21 should be used. Decided it shouldn't
P_t = P_n*cosd(v_helix) %tangential circular pitch, eqn. 13-18
p_t = pi/P_t %tangential diametral pitch, eqn. 13-4
p_x = p_t/tand(v_helix) %axial diametral pitch, 13-7
m_F = F/p_x %face contact ratio, eqn 14-22
m_N = p_N/(0.95*Z) %load sharing ratio, eqn 14-21

I = cosd(phi_t)*sind(phi_t)/(2*m_N)*(m_g)/(m_g+1) %geometry factor of pitting resistance, 14-23

%dynamic factor
Q_v = 5 %quality number, 3-7 is most commercial gears and 8-12 is precision, section eqn 14-7
B = 0.25*(12-Q_v)^(2/3) %constant for quality number, eqn 14-28
A = 50+56*(1-B) % constant for quality number, eqn 14-28
K_v = ((A+sqrt(V))/A)^B %dynamic factor, eqn 14-27,

%Compressive Stress
C_mc = 1 %eqn 14-31, uncrowned teeth
C_pf = F/(10*d_p)-0.0375+0.0125*F %pinion proportional factor, eqn 14-32
s1 = 0 %gear distance from centerline
s = 3.75 %distance between center of bearings
C_pm = 1 %pinion proportional factor, eqn 14-31
A = .127 %Table 14-9, commercial enclosed unit
B = .0158 %Table 14-9, commercial enclosed unit
C = -.930*10^-4 %Table 14-9, commercial enclosed unit
C_ma = A+B*F+C*F^2 %Mesh alignment factor, eqn 14-34
C_e = 1 %Mesh alignment correction factor, eqn. 14-35, not adjusted at assembly or lapped
K_m = 1+C_mc*(C_pf*C_pm+C_ma*C_e) %load distribution factor, eqn 14-30
C_p = 2300 %elastic coefficient, Table 14-8 for steel gear and pinion
K_o = 1.5 %overload factors, Figure 14-17
K_s = 1 %size factor, given
Cf = 1 %surface condition factor, given
sigma_comp = C_p*sqrt(w_t*K_o*K_v*K_s*K_m/(F*d_p)*Cf/I) %contact stress, 14-16, contact stress

%bending stress
dshaft_p = .75; rshaft_p = dshaft_p/2; %diameter of shaft from iterations 
dshaft_g = 1.1; rshaft_g = dshaft_g/2; %diameter of from iterations 
P_d  = P_n/cosd(phi_t) %transverse diametral pitch, eqn 13-6
K_B_p = 3; K_B_g = 1; %rim thickness factor, figure 14-16
J_p = .46*.97 %bending strength, Figure 14-7 & 14-8
J_g = .52*.93 %bending strength, Figure 14-7 & 14-8
sigma_p = w_t*K_o*K_v*K_s*(P_d/F)*(K_m*K_B_p/J_p) %bending stress, AGMA eqn 14-15
sigma_g = w_t*K_o*K_v*K_s*(P_d/F)*(K_m*K_B_g/J_g) %bending stress, AGMA eqn 14-15

%lifetime
L_p = Lifetime*60*w_p %dimensionless lifetime, eqn 14-15
L_g = Lifetime*60*w_g %dimensionless lifetime, eqn 14-15
Z_N_p = 2.466*L_p^-.056 %stress-cycle factor for pitting resistance, Figure 14-15
Z_N_g = 2.466*L_g^-.056 %stress-cycle factor for pitting resistance, Figure 14-15

%hardness ratio factor
C_H_p = 1 %eqn 14-37
HB_P = 325 %selected
HB_G = 280 %selected
HardnessRatio = (HB_P/HB_G) %Table 14-5
Aprime = 0 %Table 14-5, piecewise(AAA<=1.2,0,1.2<AAA&&AAA<1.7,8.98*10^-3*(H_BP/H_BG)-8.29*10^-3,AAA>1.7,.00698)
C_H_g = 1+Aprime*(m_g-1) %eqn 14-37
K_t = 1 %temperature factor, Section 14-15
K_r = .85; %reliability factor, Table 14-10, based on R=0.9

%required hardness values
S_c_p = SH_wear*sigma_comp/Z_N_p %AGMA surface endurance strength, eqn 14-18
S_c_g = SH_wear*sigma_comp/Z_N_g %AGMA surface endurance strength, eqn 14-18
S_c_MATERIAL_p = 349*HB_P+34300 %material surface endurance strength, figure 14-5, grade 2 through hardened
S_c_MATERIAL_g = 349*HB_G+34300 %material surface endurance strength, figure 14-5, grade 2 through hardened

%wear strength
WEAR_STRENGTH_p = S_c_MATERIAL_p*Z_N_p*C_H_p/(K_t*K_r) %safety factor, pitting eqn 14-42
WEAR_STRENGTH_g = S_c_MATERIAL_g*Z_N_p*C_H_g/(K_t*K_r) %safety factor, pitting eqn 14-42

%Gear Bending Stength, St
St_p = 108.6*HB_P+15890 %gear bending strength, eqn 14-42
St_g = 108.6*HB_G+15890 %gear bending strength, eqn 14-42

%gear contanct endurance strength
sigma_call_p = S_c_p*Z_N_p/(SH_wear*K_t*K_r) %AGMA allowable contact stress, eqn 14-18
sigma_call_g = S_c_g*Z_N_g*C_H_g/(SH_wear*K_t*K_r) %AGMA allowable contact stress, eqn 14-18

%stress cycle factor
Y_N_p = (1.3558*L_p^-.0178+1.6831*L_p^-.0323)/2; %Stress-cycle factor for bending strength, figure 14-14
Y_N_g = (1.3558*L_g^-.0178+1.6831*L_g^-.0323)/2; %Stress-cycle factor for bending strength, figure 14-14

%Gear Bending Endurance Strength
sigma_all_p = St_p*Y_N_p/(SF_bend*K_t*K_r) %Allowable bending stress, AGMA, eqn 14-17
sigma_all_g = St_g*Y_N_p/(SF_bend*K_t*K_r) %Allowable bending stress, AGMA, eqn 14-17


%actual calculated safety factor
SH_p = S_c_MATERIAL_p*Z_N_p*C_H_p/(K_t*K_r)/sigma_comp %safety factor, pitting eqn 14-42
SH_g = S_c_MATERIAL_g*Z_N_g*C_H_g/(K_t*K_r)/sigma_comp %safety factor, pitting eqn 14-42

%Bending factor of safety
SF_p = St_p*Y_N_p/(K_t*K_r)/sigma_p %safety factor, bending eqn 14-41
SF_g = St_g*Y_N_g/(K_t*K_r)/sigma_g %safety factor, bending eqn 14-41
