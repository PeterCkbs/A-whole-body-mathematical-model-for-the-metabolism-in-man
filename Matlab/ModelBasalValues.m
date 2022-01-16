function [x0,I,gamma] = ModelBasalValues(G0,I0,Gamma0,I,gamma,Q)
    %Script that loads the intial values for the differential equations.
    %It also saves some of the values as basal values, that are used as
    %parameters in the model

    %Consists of 141 differential equations
    x0 = ones(141,1);
    
    gamma.G_B_H = G0; %Heart glucose at t=0
    gamma.I_B_H = I0; %Heart insulin at t=0
    
    %%% Insulin submodel initial values %%%
    X_B = (18*G0)^I.beta_pir1 / (I.beta_pir2^I.beta_pir1 +...
        I.beta_pir3 * (18*G0)^I.beta_pir4);
    P_inf = X_B^I.beta_pir5;
    Y_B = P_inf;
    H = I.K;
    %%% Differential equations affecting insulin release %%%
    %P
    x0(139) = P_inf;
    %II
    x0(140) = X_B;
    %QQ
    x0(141) = (H * I.Q_0 + I.y * P_inf)/(H + I.M_1 * Y_B);

    % Circulating insulin initial values %
    %Heart
    I_H = I0;
    %Brain
    I_B = I_H;
    %Gut
    I_G = I_H;
    %Kidney
    I_K = I_H * (1 - I.F_KIC);
    %Muscle
    I_MP = I_H * (1 - I.F_PIC);
    %Adipose
    I_AP = I_H * (1 - I.F_PIC);
    %Liver
    I_L = 1 / Q.L * (Q.H * I_H - Q.B * I_B - Q.K * I_K...
        - Q.MP * I_MP - Q.AP * I_AP);

    % Basal values kept for later use %
    I.S_B = I.M_1 * Y_B * x0(141);
    I.r_B_PIR = Q.L / (1 - I.F_LIC) * I_L - Q.G * I_G - Q.A * I_H;
    
    
    %%% Glucagon submodel initial and basal values %%%
    gamma.r_B_PgammaR = gamma.r_MgammaC * Gamma0;
    gamma.r_PgammaC = gamma.r_MgammaC * Gamma0;
    gamma.M_G_PgammaR = 2.93 - 2.10 * tanh(4.18 * (G0 / gamma.G_B_H - 0.61));
    gamma.M_I_PgammaR = 1.31 - 0.61 * tanh(1.06 * (I_H / gamma.I_B_H - 0.47));
    gamma.r_PgammaR = gamma.M_G_PgammaR * gamma.M_I_PgammaR * gamma.r_B_PgammaR;
    % Circulating glucagon basal values %
    gamma_H = Gamma0;
    gamma_B = gamma_H;
    gamma_G = gamma_H;
    gamma_K = gamma_H;
    gamma_MP = gamma_H;
    gamma_AP = gamma_H;
    gamma_L = (gamma.r_PgammaR - gamma.r_PgammaC + Q.A * gamma_H +...
        Q.G * gamma_G) / Q.L;
    
    
    
    
    %The vector:
    %[GLC,G6P,GLY,GA3P,PYR,ACoA,OXA,CIT,LAC,AA,FFA,TGL,GLR,KET,PRO,TGL_AP,INS,GLU] 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Metabolites for the main model %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %All metabolites are in [mmol/L]
    %%% Brain initial values %%%
    x0(1:18) = [4.44, 0.0187, 0, 0.214, 0.135, 0.005, 0.002, 0.0376,...
        1.24, 2.91, 0.453, 0.7, 0.093, 0.0055, 0, 0, I_B, gamma_B]; 
    %%% Heart initial values %%%
    x0(19:36) = [G0, 0.0013, 0, 0.047, 0.132, 0.006, 0.002, 0.0143,...
        1.13, 2.91, 0.453, 0.7, 0.093, 0.0073, 0, 0, I_H, gamma_H];
    %%% Gut initial values %%%
    x0(37:54) = [4.93, 0.0009, 0, 0.035, 0.195, 0.002, 0.003, 0.03,...
        1.17, 2.65, 0.452, 0.7, 0.093, 0.0073, 0, 0, I_G, gamma_G];
    %%% Liver initial values %%%
    x0(55:72) = [5.36, 0.0212, 221.3, 0.016, 0.22, 0.005, 0.004, 0.055,...
        0.97, 2.32, 0.446, 0.715, 0.03, 0.0092, 0, 0, I_L, gamma_L];
    %%% Kidney initial values %%%
    x0(73:90) = [5.1, 0.001, 0, 0.047, 0.239, 0.008, 0.004, 0.024,...
        0.94, 2.83, 0.451, 0.7, 0.055, 0.0073, 0, 0, I_K, gamma_K];
    %%% Muscle initial values %%%
    x0(91:108) = [4.9 0.0084, 285.63, 0.139, 2.96, 0.0004, 0.04, 0.048,...
        1.28, 3.59, 0.449, 0.7, 0.093, 0.0066, 11541, 0, I_MP, gamma_MP];
    %%% Adipose initial values
    x0(109:126) = [4.85, 0.0006, 0, 0.018, 0.012, 0.0683, 0.0002, 0.0055,...
        1.34, 2.91, 0.518, 0.64, 0.5, 0.0073, 0, 5978, I_AP, gamma_AP];
    %%% Basal values kept for later use %%%
    
    
    %%% SIMO Model initial values %%%
    x0(127:138) = 0; %Nothing at fasting steady state.    
    
    %%% Basal values for hormonal control %%%
    %Tissue specific hormonal control
    gamma.L_SS = gamma_L;
    gamma.AP_SS = gamma_AP;
    I.L_SS = I_L;
    I.MP_SS = I_MP;
    I.AP_SS = I_AP;
    
    
    
end