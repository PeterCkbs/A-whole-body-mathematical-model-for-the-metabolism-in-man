function [Q,V,Km,Vm,p,I,gamma,mu] = LoadParmModel()
%Script, that loads all the parameters in the model.
%It requires no input
%The parameters are divided into 8 different structures for easier use
%Q: Flowrates in main model
%V: Volumes of compartments in main model
%Km: All limiting velocities in the model
%Vm: All maximum velocities in the model
%p: Parameters in GI, as well as the stoichiometric matrix
%I: Parameters in the insulin sub model
%gamma: Parameters in the gamma sub model
%mu: Parameters associated with hormonal control


%Example of parameter loadout
    %Flowrate of heart for metabolites
    %Q.H
    %Volume of heart
    %V.H
    %Km of GLC->G6P in heart
    %Km.H_GLC_G6P
    %Vmax of GLC->G6P in liver
    %Vm.L_GLC_G6P
    %Hormonal control of glucose uptake in muscle
    %mu.MP_GLC_G6P
    
    
    %Vector of the metabolites in the model
    %[GLC, G6P, GLY, GA3P, PYR, ACoA, OXA, CIT, LAC, AA, FFA, TGL, GLR, KET, PRO, TGL_AP, INS, GLU]
    %Stoichiometric Matrix
    S = zeros(31,18);
    %GLC -> G6P
    S(1,1) = -1; S(1,2) = 1;
    %G6P -> GLC
    S(2,2) = -1; S(2,1) = 1;
    %G6P -> GA3P
    S(3,2) = -1; S(3,4) = 2; 
    %GA3P -> G6P 
    S(4,4) = -2; S(4,2) = 1;
    %G6P -> GLY
    S(5,2) = -1; S(5,3) = 1;
    %GLY -> G6P
    S(6,3) = -1; S(6,2) = 1;
    %GA3P -> PYR
    S(7,4) = -1; S(7,5) = 1; 
    %PYR -> GA3P 
    S(8,5) = -1; S(8,4) = 1; 
    %PYR -> LAC
    S(9,5) = -1; S(9,9) = 1;
    %LAC -> PYR
    S(10,9) = -1; S(10,5) = 1; 
    %PYR -> AA
    S(11,5) = -1; S(11,10) = 1;
    %AA -> PYR
    S(12,10) = -1; S(12,5) = 1;
    %PYR -> ACoA 
    S(13,5) = -1; S(13,6) = 1;
    %ACoA + OXA -> CIT
    S(14,6) = -1; S(14,7) = -1; S(14,8) = 1;
    %CIT -> OXA
    S(15,8) = -1; S(15,7) = 1;
    %OXA -> PYR
    S(16,7) = -1; S(16,5) = 1;
    %PYR -> OXA
    S(17,5) = -1; S(17,7) = 1;
    %GA3P -> GLR
    S(18,4) = -1; S(18,13) = 1;
    %GLR -> GA3P
    S(19,13) = -1; S(19,4) = 1;
    %3FFA + GLR -> TGL
    S(20,11) = -3; S(20,13) = -1; S(20,12) = 1;
    %TGL -> 3FFA + GLR
    S(21,12) = -1; S(21,11) = 3; S(21,13) = 1;
    %FFA -> ACoA
    S(22,11) = -1; S(22,6) = 7;
    %ACoA -> FFA
    S(23,6) = -7; S(23,11) = 1;
    %KET -> ACoA
    S(24,14) = -1; S(24,6) = 1;
    %ACoA -> KET
    S(25,6) = -1; S(25,14) = 1;
    %AA -> PRO
    S(26,10) = -1; S(26,15) = 1;
    %PRO -> AA
    S(27,15) = -1; S(27,10) = 1;
    %3FFA + GLR -> TGL_AP
    S(28,11) = -3; S(28,13) = -1; S(28,16) = 1;
    %TGL_AP -> 3FFA + GLR
    S(29,16) = -1; S(29,11) = 3; S(29,13) = 1;
    %Insulin -> Formation and Breakdown
    S(30,17) = 1;
    %Glucagon -> Formation and Breakdown
    S(31,18) = 1;
    %Convert to structure
    p.S = S;
    
    %Distribution matrix consisting of 1's in the diagonal for the 
    %circulating metabolites
    %[GLC,LAC,AA,FFA,TGL,GLR,KET] are the circulating metabolites.
    m = zeros(18,1); m([1,9:14,17,18]) = 1; p.M = diag(m);
    
    %T-vectors for all the tissues
    %Brain - 11 reactions
    T_B = zeros(11,31); T_B(1,1) = 1; T_B(2,3) = 1; T_B(3,7) = 1;
    T_B(4,9) = 1; T_B(5,10) = 1; T_B(6,13) = 1; T_B(7,14) = 1; 
    T_B(8,15) = 1; T_B(9,16) = 1; T_B(10,17) = 1; T_B(11,24) = 1;
    p.T_B = T_B;
    %Heart - 13 reactions
    T_H = zeros(13,31); T_H(1,1) = 1; T_H(2,3) = 1; T_H(3,7) = 1;
    T_H(4,9) = 1; T_H(5,10) = 1; T_H(6,13) = 1; T_H(7,14) = 1; 
    T_H(8,15) = 1; T_H(9,16) = 1; T_H(10,17) = 1; T_H(11,21) = 1;
    T_H(12,22) = 1; T_H(13,24) = 1;
    p.T_H = T_H;
    %Gut - 12 reactions
    T_G = zeros(12,31); T_G(1,1) = 1; T_G(2,3) = 1; T_G(3,7) = 1;
    T_G(4,9) = 1; T_G(5,10) = 1; T_G(6,12) = 1; T_G(7,13) = 1; 
    T_G(8,14) = 1; T_G(9,15) = 1; T_G(10,16) = 1; T_G(11,17) = 1;
    T_G(12,22) = 1;
    p.T_G = T_G;
    %Liver - 21 reactions + insulin and glucagon production/clearance
    T_L = zeros(23,31); T_L(1,1) = 1; T_L(2,2) = 1; T_L(3,3) = 1; 
    T_L(4,4) = 1; T_L(5,5) = 1; T_L(6,6) = 1; T_L(7,7) = 1; T_L(8,8) = 1;
    T_L(9,9) = 1; T_L(10,10) = 1; T_L(11,12) = 1; T_L(12,13) =1;
    T_L(13,14) = 1; T_L(14,15) = 1; T_L(15,16) = 1; T_L(16,17) = 1;
    T_L(17,19) = 1; T_L(18,20) = 1; T_L(19,22) = 1; T_L(20,23) = 1;
    T_L(21,25) = 1; T_L(22,30) = 1; T_L(23,31) = 1;
    p.T_L = T_L;
    %Kidney - 16 reactions + insulin clearance
    T_K = zeros(17,31); T_K(1,1) = 1; T_K(2,2) = 1; T_K(3,3) = 1;
    T_K(4,4) = 1; T_K(5,7) = 1; T_K(6,8) = 1; T_K(7,9) = 1; T_K(8,10) = 1;
    T_K(9,12) = 1; T_K(10,13) = 1; T_K(11,14) = 1; T_K(12,15) = 1;
    T_K(13,16) = 1; T_K(14,17) = 1; T_K(15,19) = 1; T_K(16,22) = 1;
    T_K(17,30) = 1;
    p.T_K = T_K;
    %Muscle - 18 reactions + insulin clearance
    T_MP = zeros(19,31); T_MP(1,1) = 1; T_MP(2,3) = 1; T_MP(3,5) = 1;
    T_MP(4,6) = 1; T_MP(5,7) = 1; T_MP(6,9) = 1; T_MP(7,10) = 1; T_MP(8,11) = 1;
    T_MP(9,13) = 1; T_MP(10,14) = 1; T_MP(11,15) = 1; T_MP(12,16) = 1;
    T_MP(13,17) = 1; T_MP(14,21) = 1; T_MP(15,22) = 1; T_MP(16,24) = 1;
    T_MP(17,26) = 1; T_MP(18,27) = 1; T_MP(19,30) = 1;
    p.T_MP = T_MP;
    %Adipose - 16 reactions + insulin clearance
    T_AP = zeros(17,31); T_AP(1,1) = 1; T_AP(2,3) = 1; T_AP(3,7) = 1;
    T_AP(4,9) = 1; T_AP(5,10) = 1; T_AP(6,13) = 1; T_AP(7,14) = 1; 
    T_AP(8,15) = 1; T_AP(9,16) = 1; T_AP(10,17) = 1; T_AP(11,18) = 1;
    T_AP(12,21) = 1; T_AP(13,22) = 1; T_AP(14,23) = 1; T_AP(15,28) = 1;
    T_AP(16,29) = 1; T_AP(17,30) = 1;
    p.T_AP = T_AP;
    
    % The Q's in [L/min]
    Q.B = 0.59; Q.H = 4.37; Q.A = 0.25; Q.L = 1.26;
    Q.G = 1.01; Q.K = 1.01; Q.MP = 1.2231; Q.AP = 0.2869;
    % The V's in [L]
    V.B = 0.8; V.H = 1.38; V.L = 2.51; V.G = 1.12;
    V.K = 0.66; V.AP = 1.945; V.MP = 5.835;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% KINETICS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Brain kinetics %%%%%%
    %Vm: [mM/min], Km: [mmol/L]
    %Reaction 1
    Vm.B_GLC_G6P = 0.3963/V.B;
    Km.B_GLC_G6P = 1;
    %Reaction 3
    Vm.B_G6P_GA3P = 3.963e-1/V.B;
    Km.B_G6P_GA3P = 4.2e-3;
    %Reaction 7
    Vm.B_GA3P_PYR = 6.27/V.B;
    Km.B_GA3P_PYR = 1.86;
    %Reaction 9 
    Vm.B_PYR_LAC = 0.5571/V.B;
    Km.B_PYR_LAC = 0.065;
    %Reaction 10
    Vm.B_LAC_PYR = 1.4453/V.B;
    Km.B_LAC_PYR = 4.5068;
    %Reaction 13
    Vm.B_PYR_ACoA = 0.6953*2/V.B;
    Km.B_PYR_ACoA = 0.187;
    %Reaction 14
    Vm.B_ACoA_OXA_CIT = 0.6953 * 1.2/V.B;
    Km.B_ACoA_OXA_CIT = 4.5e-6;
    %Reaction 15
    Vm.B_CIT_OXA = 39.34/V.B;
    Km.B_CIT_OXA = 2.51;
    %Reaction 16
    Vm.B_OXA_PYR = 1e-2/V.B;
    Km.B_OXA_PYR = 0.003;
    %Reaction 17
    Vm.B_PYR_OXA = 1e-2/V.B;
    Km.B_PYR_OXA = 0.187;
    %Reaction 24
    Vm.B_KET_ACoA = 1/V.B;
    Km.B_KET_ACoA = 5;
    
    
    %%%%%% Heart kinetics %%%%%%
    %Vm: [mM/min], Km: [mmol/L]
    %Reaction 1
    Vm.H_GLC_G6P = 0.1243/V.H;
    Km.H_GLC_G6P = 3;
    %Reaction 3
    Vm.H_G6P_GA3P = 3.28e-1/V.H;
    Km.H_G6P_GA3P = 4.2e-3;
    %Reaction 7
    Vm.H_GA3P_PYR = 6.27/V.H;
    Km.H_GA3P_PYR = 1.86;
    %Reaction 9 
    Vm.H_PYR_LAC = 0.6116/V.H;
    Km.H_PYR_LAC = 0.1219;
    %Reaction 10
    Vm.H_LAC_PYR = 0.7928/V.H;
    Km.H_LAC_PYR = 3.3714;
    %Reaction 13
    Vm.H_PYR_ACoA = 0.0454 * 2/V.H;
    Km.H_PYR_ACoA = 0.187;
    %Reaction 14
    Vm.H_ACoA_OXA_CIT = 0.2562 * 1.2/V.H;
    Km.H_ACoA_OXA_CIT = 4.5e-6;
    %Reaction 15
    Vm.H_CIT_OXA = 39.34/V.H;
    Km.H_CIT_OXA = 2.51;
    %Reaction 16
    Vm.H_OXA_PYR = 1e-2/V.H;
    Km.H_OXA_PYR = 0.003;
    %Reaction 17
    Vm.H_PYR_OXA = 1e-2/V.H;
    Km.H_PYR_OXA = 0.187;
    %Reaction 21
    Vm.H_TGL_GLR_FFA = 0.0032/V.H;
    Km.H_TGL_GLR_FFA = 11.21;
    %Reaction 22
    Vm.H_FFA_ACoA = 0.02645 * 2/V.H;
    Km.H_FFA_ACoA = 0.45;
    %Reaction 24
    Vm.H_KET_ACoA = 0.02562/V.H;
    Km.H_KET_ACoA = 0.5;
    
    
    %%%%%% Gut kinetics %%%%%%
    %Vm: [mM/min], Km: [mmol/L]
    %Reaction 1
    Vm.G_GLC_G6P = 0.2564/V.G;
    Km.G_GLC_G6P = 17;
    %Reaction 3
    Vm.G_G6P_GA3P = 3.28e-1/V.G;
    Km.G_G6P_GA3P = 4.2e-3;
    %Reaction 7
    Vm.G_GA3P_PYR = 6.27/V.G;
    Km.G_GA3P_PYR = 1.86;
    %Reaction 9 
    Vm.G_PYR_LAC = 0.6222/V.G;
    Km.G_PYR_LAC = 0.0391;
    %Reaction 10
    Vm.G_LAC_PYR = 0.9748/V.G;
    Km.G_LAC_PYR = 1.1989;
    %Reaction 12 
    Vm.G_AA_PYR = 0.3766 * 3/V.G;
    Km.G_AA_PYR = 2.85 * 3;
    %Reaction 13
    Vm.G_PYR_ACoA = 0.4476 * 3/V.G;
    Km.G_PYR_ACoA = 0.187 * 3;
    %Reaction 14
    Vm.G_ACoA_OXA_CIT = 0.5968 * 1.5/V.G;
    Km.G_ACoA_OXA_CIT = 4.5e-6;
    %Reaction 15
    Vm.G_CIT_OXA = 39.34/V.G;
    Km.G_CIT_OXA = 2.51;
    %Reaction 16
    Vm.G_OXA_PYR = 1e-2/V.G;
    Km.G_OXA_PYR = 0.003;
    %Reaction 17 
    Vm.G_PYR_OXA = 1e-2/V.G;
    Km.G_PYR_OXA = 0.187;
    %Reaction 22
    Vm.G_FFA_ACoA = 0.0213 * 3/V.G;
    Km.G_FFA_ACoA = 0.45 * 3;
    
    
    %%%%%% Liver kinetics %%%%%%
    %Vm: [mM/min], Km: [mmol/L]
    %Reaction 1
    Vm.L_GLC_G6P = 0.4274/V.L;
    Km.L_GLC_G6P = 17;
    %Reaction 2
    Vm.L_G6P_GLC = 0.622/V.L;
    Km.L_G6P_GLC = 1e-5;
    %Reaction 3
    Vm.L_G6P_GA3P = 3.28e-1/V.L;
    Km.L_G6P_GA3P = 0.2575;
    %Reaction 4
    Vm.L_GA3P_G6P = 0.2331/V.L;
    Km.L_GA3P_G6P = 1e-5;
    %Reaction 5
    Vm.L_G6P_GLY = 0.6816/V.L;
    Km.L_G6P_GLY = 0.2575;
    %Reaction 6
    Vm.L_GLY_G6P = 0.389 * 2/V.L;
    Km.L_GLY_G6P = 221.3;
    %Reaction 7
    Vm.L_GA3P_PYR = 6.27/V.L;
    Km.L_GA3P_PYR = 1.86;
    %Reaction 8
    Vm.L_PYR_GA3P = 0.3857 * 2/V.L;
    Km.L_PYR_GA3P = 0.187;
    %Reaction 9 
    Vm.L_PYR_LAC = 0.8949/V.L;
    Km.L_PYR_LAC = 0.2414;
    %Reaction 10
    Vm.L_LAC_PYR = 1.4099/V.L;
    Km.L_LAC_PYR = 1.0518;
    %Reaction 12
    Vm.L_AA_PYR = (0.1856 + 0.343) * 2/V.L;
    Km.L_AA_PYR = 2.85;
    %Reaction 13
    Vm.L_PYR_ACoA = 0.343 * 2/V.L;
    Km.L_PYR_ACoA = 0.187;
    %Reaction 14
    Vm.L_ACoA_OXA_CIT = 0.9881 * 1.1/V.L;
    Km.L_ACoA_OXA_CIT = 4.5e-6;
    %Reaction 15
    Vm.L_CIT_OXA = 39.34/V.L;
    Km.L_CIT_OXA = 2.51;
    %Reaction 16
    Vm.L_OXA_PYR = 1e-2/V.L;
    Km.L_OXA_PYR = 0.003;
    %Reaction 17
    Vm.L_PYR_OXA = 1e-2/V.L;
    Km.L_PYR_OXA = 0.187;
    %Reaction 19
    Vm.L_GLR_GA3P = 0.0806 * 2/V.L;
    Km.L_GLR_GA3P = 0.05;
    %Reaction 20
    Vm.L_GLR_FFA_TGL = 0.0213 * 2/V.L;
    Km.L_GLR_FFA_TGL = 0.0178;
    %Reaction 22
    Vm.L_FFA_ACoA = 0.3/V.L;
    Km.L_FFA_ACoA = 1.385;
    %Reaction 23
    Vm.L_ACoA_FFA = 0.3/V.L;
    Km.L_ACoA_FFA = 0.5;
    %Reaction 25
    Vm.L_ACoA_KET = 1/V.L;
    Km.L_ACoA_KET = 2;
    
    %%%%%% Kidney kinetics %%%%%%
    %Vm: [mM/min], Km: [mmol/L]
    %Reaction 1
    Vm.K_GLC_G6P = 0.2564/V.K;
    Km.K_GLC_G6P = 17;
    %Reaction 2
    Vm.K_G6P_GLC = 0.1554/V.K;
    Km.K_G6P_GLC = 1e-5;
    %Reaction 3
    Vm.K_G6P_GA3P = 3.28e-1/V.K; 
    Km.K_G6P_GA3P = 4.2e-3; 
    %Reaction 4
    Vm.K_GA3P_G6P = 0.1554/V.K; 
    Km.K_GA3P_G6P = 1e-5; 
    %Reaction 7
    Vm.K_GA3P_PYR = 6.27/V.K; 
    Km.K_GA3P_PYR = 1.86; 
    %Reaction 8
    Vm.K_PYR_GA3P = 0.2741*2/V.K; 
    Km.K_PYR_GA3P = 0.187; 
    %Reaction 9 
    Vm.K_PYR_LAC = 0.2454/V.K; 
    Km.K_PYR_LAC = 0.3045; 
    %Reaction 10
    Vm.K_LAC_PYR = 0.4144/V.K; 
    Km.K_LAC_PYR = 0.3239; 
    %Reaction 12
    Vm.K_AA_PYR = 0.082 * 2/V.K; 
    Km.K_AA_PYR = 2.85; 
    %Reaction 13
    Vm.K_PYR_ACoA = 0.1166 * 2/V.K;
    Km.K_PYR_ACoA = 0.187;
    %Reaction 14
    Vm.K_ACoA_OXA_CIT = 0.3660 * 1.2/V.K;
    Km.K_ACoA_OXA_CIT = 4.5e-6;
    %Reaction 15
    Vm.K_CIT_OXA = 39.34/V.K;
    Km.K_CIT_OXA = 2.51;
    %Reaction 16
    Vm.K_OXA_PYR = 1e-2/V.K;
    Km.K_OXA_PYR = 0.003;
    %Reaction 17
    Vm.K_PYR_OXA = 1e-2/V.K;
    Km.K_PYR_OXA = 0.187;
    %Reaction 19
    Vm.K_GLR_GA3P = 0.0368 * 2/V.K;
    Km.K_GLR_GA3P = 0.05;
    %Reaction 22
    Vm.K_FFA_ACoA = 0.0356 * 2/V.K;
    Km.K_FFA_ACoA = 0.45;
    
    
    %%%%%% Muscle periphery kinetics %%%%%%
    %Vm: [mM/min], Km: [mmol/L]
    %Reaction 1
    Vm.MP_GLC_G6P = 0.2331/V.MP;
    Km.MP_GLC_G6P = 5;
    %Reaction 3
    Vm.MP_G6P_GA3P = 3.28e-1/V.MP;
    Km.MP_G6P_GA3P = 4.2e-3;
    %Reaction 5
    Vm.MP_G6P_GLY = 4.987e-1/V.MP;
    Km.MP_G6P_GLY = 8.93e-2/V.MP;
    %Reaction 6
    Vm.MP_GLY_G6P = 4.437e-1/V.MP;
    Km.MP_GLY_G6P = 285.6 * 2;
    %Reaction 7
    Vm.MP_GA3P_PYR = 6.27/V.MP;
    Km.MP_GA3P_PYR = 1.86;
    %Reaction 9
    Vm.MP_PYR_LAC = 0.821/V.MP;
    Km.MP_PYR_LAC = 0.5029;
    %Reaction 10
    Vm.MP_LAC_PYR = 0.6852/V.MP;
    Km.MP_LAC_PYR = 0.3821;
    %Reaction 11
    Vm.MP_PYR_AA = 7.761e-2/V.MP;
    Km.MP_PYR_AA = 7.43e-4;
    %Reaction 13
    Vm.MP_PYR_ACoA = 0.09892 * 2/V.MP;
    Km.MP_PYR_ACoA = 0.187;
    %Reaction 14
    Vm.MP_ACoA_OXA_CIT = 0.6587 * 1.5/V.MP;
    Km.MP_ACoA_OXA_CIT = 4.5e-6;
    %Reaction 15
    Vm.MP_CIT_OXA = 39.34/V.MP;
    Km.MP_CIT_OXA = 2.51;
    %Reaction 16
    Vm.MP_OXA_PYR = 1e-2/V.MP;
    Km.MP_OXA_PYR = 0.003;
    %Reaction 17
    Vm.MP_PYR_OXA = 1e-2/V.MP;
    Km.MP_PYR_OXA = 0.187;
    %Reaction 21
    Vm.MP_TGL_GLR_FFA = 0.0032/V.MP;
    Km.MP_TGL_GLR_FFA = 11.21;
    %Reaction 22
    Vm.MP_FFA_ACoA = 0.08004 * 2/V.MP;
    Km.MP_FFA_ACoA = 0.45;
    %Reaction 25
    Vm.MP_KET_ACoA = 0.06587/V.MP;
    Km.MP_KET_ACoA = 0.5;
    %Reaction 26
    Vm.MP_AA_PRO = 1.92/V.MP;
    Km.MP_AA_PRO = 2.85*4;
    %Reaction 27
    Vm.MP_PRO_AA = 2.4/V.MP;
    Km.MP_PRO_AA = 11541;
    
    
    %%%%%% Adipose kinetics %%%%%%
    %Vm: [mM/min], Km: [mmol/L]
    %Reaction 1
    Vm.AP_GLC_G6P = 0.0777/V.AP;
    Km.AP_GLC_G6P = 5;
    %Reaction 3
    Vm.AP_G6P_GA3P = 3.28e-1/V.AP;
    Km.AP_G6P_GA3P = 4.2e-3;
    %Reaction 7
    Vm.AP_GA3P_PYR = 6.27/V.AP;
    Km.AP_GA3P_PYR = 1.86;
    %Reaction 9
    Vm.AP_PYR_LAC = 0.5166/V.AP;
    Km.AP_PYR_LAC = 0.0022;
    %Reaction 10
    Vm.AP_LAC_PYR = 0.9326/V.AP;
    Km.AP_LAC_PYR = 1.9798;
    %Reaction 13
    Vm.AP_PYR_ACoA = 0.00777 * 11/V.AP;
    Km.AP_PYR_ACoA = 0.187 * 10;
    %Reaction 14
    Vm.AP_ACoA_OXA_CIT = 0.098474 * 1.2/V.AP;
    Km.AP_ACoA_OXA_CIT = 4.5e-6;
    %Reaction 15
    Vm.AP_CIT_OXA = 39.34/V.AP;
    Km.AP_CIT_OXA = 2.51;
    %Reaction 16
    Vm.AP_OXA_PYR = 1e-2/V.AP;
    Km.AP_OXA_PYR = 0.003;
    %Reaction 17 
    Vm.AP_PYR_OXA = 1e-2/V.AP;
    Km.AP_PYR_OXA = 0.187;
    %Reaction 18
    Vm.AP_GA3P_GLR = 0.0155 * 2/V.AP;
    Km.AP_GA3P_GLR = 0.016;
    %Reaction 21
    Vm.AP_TGL_GLR_FFA = 1.7897/V.AP;
    Km.AP_TGL_GLR_FFA = 63.525;
    %Reaction 22
    Vm.AP_FFA_ACoA = 0.01296 * 2/V.AP;
    Km.AP_FFA_ACoA = 0.45;
    %Reaction 23
    Vm.AP_ACoA_FFA = 0.01296 * 2/V.AP;
    Km.AP_ACoA_FFA = 1;
    %Reaction 28
    Vm.AP_GLR_FFA_TGL_AP = 0.3578/V.AP;
    Km.AP_GLR_FFA_TGL_AP = 0.1368;
    %Reaction 29, lipolysis
    Vm.AP_TGL_AP_GLR_FFA = 0.3158 * 2/V.AP;
    Km.AP_TGL_AP_GLR_FFA = 5978;
    
    
    
    %%%% Parameters in the SIMO-Model %%%%
    p.GI_k_js = 0.028237;
    p.GI_k_gl = 0.0180942;
    p.GI_k_gj = 0.0329673; 
    p.GI_k_rj = 0.0344046;
    p.GI_k_lr = 0.0513802;
    
    
    %%%% Parameters in the Insulin submodel %%%%
    I.F_LIC = 0.4; %[mU/min]
    I.F_KIC = 0.3; %[mU/min]
    I.F_PIC = 0.15; %[mU/min]
    %Insulin release model
    I.beta_pir1 = 3.27; I.beta_pir3 = 5.93; I.beta_pir4 = 3.02;
    I.beta_pir5 = 1.11;
    I.beta_pir2 = 132; %[mg/dL]
    I.M_1 = 0.00747; I.M_2 = 0.0958; I.alpha = 0.0482; %[1/min]
    I.beta = 0.931; %[1/min]
    I.Q_0 = 6.33; %[U]
    I.K = 0.00794;
    I.y = 0.575;
    
    %%%% Parameters in the Glucagon submodel %%%%
    gamma.r_MgammaC = 0.91; %[L/min]
    
    %%%% Hormonal Control Parameters %%%%
    %%% LIVER %%%
    %Insulin activation
    mu.L_GLC_G6P = 4.2422; 
    mu.L_PYR_ACoA = 3; 
    mu.L_ACoA_FFA = 3;
    %Insulin-to-Glucagon stimulation:
    mu.L_G6P_GLY = 2; 
    mu.L_G6P_GA3P = 3;
    %Glucagon-to-Insulin stimulation
    mu.L_GLY_G6P = 0.5; 
    mu.L_GA3P_G6P = 0.5; 
    %%% MUSCLE %%%
    %Insulin activation
    mu.MP_GLC_G6P = 2.1557; 
    mu.MP_G6P_GLY = 3; 
    mu.MP_G6P_GA3P = 0.1; 
    mu.MP_PYR_ACoA = 0.5; 
    mu.MP_AA_PRO = 4;
    %Insulin inhibition
    mu.MP_PRO_AA = 0.25; 
    mu.MP_GLY_G6P = 0.2;
    %%% ADIPOSE %%%
    %Insulin activation
    mu.AP_GLC_G6P = 0.6352; 
    mu.AP_TGL_GLR_FFA =  4;
    mu.AP_GLR_FFA_TGL_AP =  2;
    %Glucagon-to-insulin stimulation
    mu.AP_TGL_AP_GLR_FFA =  0.1; 
    
end
    