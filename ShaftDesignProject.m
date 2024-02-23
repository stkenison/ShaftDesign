clc;clear all;close all;

%material properties for 1030 HR steel
s_ult = 68000; %ultimate strength
s_yield = 37500; %yield strength
if s_ult>200000, se_prime = 100000; %set se_prime value
else, se_prime = s_ult/2; end

%input forces, diameter OR desired safety factor for all cases
problems_to_solve = 30;
M = [[400 400 400 300 300 200 300 200 200 200] [100 100 100 100 100] [100 100 100 100 100] [300 300 300 400 300 200 200 300 300 100]];
T = [[300 300 200 300 100 100 300 100 200 100] [100 200 300 300 200] [200 200 200 200 100] [100 100 300 100 100 300 100 200 300 200]];
Diam_or_SF = [[1 .75 .75 1 1 .75 .75 .75 1 1] [1 1 1.5 1.5 1] [1.5 1 1.5 1.5 1] [1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5]]; %Diameter or safety factor desired

%define concentration types for all cases
%1 - sharp fillet, 2 - wide fillet, 3 - keyway, 4 - retaining ring groove
Concentration_Type = [[4 4 4 1 3 4 4 1 3 2] [2 4 1 3 1] [2 4 1 3 4] [3 1 3 1 2 4 2 2 2 3]];

%define problem type for all cases
%1 - SF Goodman, 2 - SF first cycle yield, 3 - SF first cycle conservative,
%4 - Solve for diameter
Problem_Type = [[1 1 1 1 1 1 1 1 1 1] [2 2 2 2 2] [3 3 3 3 3] [4 4 4 4 4 4 4 4 4 4]];

answer = zeros(problems_to_solve,1); %array to hold answers
for i = 1:problems_to_solve
    s_e = CalculateSE(se_prime,s_ult,Diam_or_SF(i)); %calculate se value
    [kf kfs] = CalculateKfKfs(Diam_or_SF(i),Concentration_Type(i),s_ult); %calculate fatigue stress concentration factors
    [sa sm] = CalcEquivCombStress(Diam_or_SF(i),M(i),T(i),kf,kfs); %calculate combined equivalent stresses
    if Problem_Type(i)==4
        answer(i) = DiamSolver(M(i),T(i),Diam_or_SF(i),Concentration_Type(i),s_yield,s_ult,se_prime); %solve for minimum diameter
    else
        answer(i) = SFsolver(M(i),T(i),Diam_or_SF(i),sa,sm,kf,kfs,s_e,s_yield,s_ult,Problem_Type(i)); %solve for safety factor
    end
    
end

[linspace(1,problems_to_solve,problems_to_solve)' answer]

%fuction to calculate endurance limit
function s_e = CalculateSE(se_prime,s_ult,D)
    k_a = (2*(s_ult/1000)^(-.217)); %surface stress concentration factor
    k_b = (D/0.3)^(-.107); %size stress concentration factor
    k_c = 1; %load stress concentration factor for torsion
    k_d = 1; %temperature stress concentation factor
    k_e = 1; %reliability stress concentration factor for 50% reliability
    s_e = se_prime*k_a*k_b*k_c*k_d*k_e;
end

%fuction to calculate fatigue stress concentration factors
function [kf kfs] = CalculateKfKfs(D,Type,s_ult)
    switch Type
        case 1 %sharp shoulder fillet
            kt = 2.7;kts = 2.2; r = D*.02;
        case 2 %well-rounded shoulder fillet
            kt = 1.7;kts = 1.5; r = D*.1;
        case 3 %keyseat
            kt = 2.14;kts = 3.0;r = D*.02;
        case 4 %retaining ring groove
            kt = 5;kts = 3;r = 0.01;
    end
    
    %convert stress concentrations to fatigue stress concentrations
    %using eqn 6-32 through 6-36
    a_bending = (0.246-3.08*10^-3*(s_ult/1000)+1.51*10^-5*(s_ult/1000)^2-2.67*10^-8*(s_ult/1000)^3)^2;
    a_torsion = (0.190-2.51*10^-3*(s_ult/1000)+1.35*10^-5*(s_ult/1000)^2-2.67*10^-8*(s_ult/1000)^3)^2;
    kf = 1+(kt-1)/(1+sqrt(a_bending/r));
    kfs = 1+(kts-1)/(1+sqrt(a_torsion/r));
end

%fuction to calculate equivalent combined stresses using eqn 7-11 and 7-12
function [sa sm] = CalcEquivCombStress(D,M,T,kf,kfs)
    sa = sqrt(((32*kf*M)/(pi*D^3))^2);
    sm = sqrt(3*((16*kfs*T)/(pi*D^3))^2);
end

%function to solve for safety factor using DEGoodman Method
function n = DEGoodman(sa,sm,s_e,s_ult)
    syms n
    n = vpasolve(sa/s_e+sm/s_ult==1/n);
end

%function to solve for yielding safety factor using von mises stress 
function n = VonMisesFirstCycleYield(M,T,D,kf,kfs,s_yield)
    syms n
    n = vpasolve(s_yield/n==sqrt((M*kf*D/2/(pi*D^4/64))^2+3*(T*kfs*D/2/(pi*D^4/32))^2));
end

%function to solve for yielding safety factor using conservative estimate
function n = ConservativeFirstCycleYield(sa, sm, s_yield)
    syms n
    n = vpasolve(sa+sm==s_yield/n);
end

%Call solution function for safety factor depending on problem type
function solution = SFsolver(M,T,D,sa,sm,kf,kfs,s_e,s_yield,s_ult,problem_type)
    switch problem_type
        case 1, solution = DEGoodman(sa,sm,s_e,s_ult);
        case 2, solution = VonMisesFirstCycleYield(M,T,D,kf,kfs,s_yield);
        case 3, solution = ConservativeFirstCycleYield(sa, sm, s_yield);
    end
end

function diam = DiamSolver(M,T,SF,Type,s_yield,s_ult,se_prime)
    syms d real positive
    k_a = (2*(s_ult/1000)^(-.217)); %surface stress concentration factor
    k_b = (d/0.3)^(-.107); %size stress concentration factor
    k_c = 1; %load stress concentration factor for torsion
    k_d = 1; %temperature stress concentation factor
    k_e = 1; %reliability stress concentration factor for 50% reliability
    s_e = se_prime*k_a*k_b*k_c*k_d*k_e;

    %set kt, kts, and r values based on stress concentration type
    switch Type
        case 1 %sharp shoulder fillet
            kt = 2.7;kts = 2.2; r = d*.02;
        case 2 %well-rounded shoulder fillet
            kt = 1.7;kts = 1.5; r = d*.1;
        case 3 %keyseat
            kt = 2.14;kts = 3.0;r = d*.02;
        case 4 %retaining ring groove
            kt = 5;kts = 3;r = 0.01;
    end
    
    %convert stress concentrations to fatigue stress concentrations
    %using eqn 6-32 through 6-36
    a_bending = (0.246-3.08*10^-3*(s_ult/1000)+1.51*10^-5*(s_ult/1000)^2-2.67*10^-8*(s_ult/1000)^3)^2;
    a_torsion = (0.190-2.51*10^-3*(s_ult/1000)+1.35*10^-5*(s_ult/1000)^2-2.67*10^-8*(s_ult/1000)^3)^2;
    kf = 1+(kt-1)/(1+sqrt(a_bending/r));
    kfs = 1+(kts-1)/(1+sqrt(a_torsion/r));
    
    %calculate combined stresses
    sa = sqrt(((32*kf*M)/(pi*d^3))^2);
    sm = sqrt(3*((16*kfs*T)/(pi*d^3))^2);
    
    %compute conservative and infinite life diameters
    diam_cons = vpasolve(sa+sm==s_yield/SF,d,[0 Inf]);
    diam_inf_life = vpasolve(sa/s_e+sm/s_ult==1/SF,d,[0 Inf]);
    
    %output maximum required diameter between two methods as solution
    diam = max([diam_cons diam_inf_life]);
end


