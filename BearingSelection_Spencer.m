clc
clear all
close all

%reaction forces
Ra = 206.53682812427652465532769127602;
Rb = 325.68995382055760119396397627685;
Rc = 338.26405565207855332197237745012;
Rd = 372.33719255023731734689949456434;

%given values
Lifetime = 1000; w_input = 3000; w_output = 400;
R = 0.975; rating_life = 10^6;af = 1.4;

%lifetimes
Ld_input = Lifetime*60*w_input
Ld_output = Lifetime*60*w_output

F_r_a = Ra*((Ld_input/rating_life)/(0.02+(4.459-0.02)*log(1/R)^(1/1.483)))^(3/10) %straight roller bearing
F_r_b = Rb*((Ld_input/rating_life)/(0.02+(4.459-0.02)*log(1/R)^(1/1.483)))^(1/3) %radial bearing
F_r_c = Ra*((Ld_output/rating_life)/(0.02+(4.459-0.02)*log(1/R)^(1/1.483)))^(1/3) %radial bearing
F_r_d = Rb*((Ld_output/rating_life)/(0.02+(4.459-0.02)*log(1/R)^(1/1.483)))^(3/10) %straight roller bearing

C10_a = vpa(af*F_r_a) %straight roller bearing
C10_b = vpa(af*F_r_b) %radial bearing
C10_c = vpa(af*F_r_c) %radial bearing
C10_d = vpa(af*F_r_d) %straight roller bearing