    function [dzdt] = Model(t,x0,Q,V,Km,Vm,p,I,gamma,mu)
    %The complete model containing 141 differential equations.
    
    %%% Unpack states %%%
    
    %Equations for the circulating metabolites used in metabolic reactions
    %inside tissues.
    %The indexes follow the metabolite vector:
    %[GLC, G6P, GLY, GA3P, PYR, ACoA, OXA, CIT, LAC, AA, FFA, TGL, GLR, KET, PRO, TGL_AP, INS, GLU]
    C_B = x0(1:18);
    C_H = x0(19:36);
    C_G = x0(37:54);
    C_L = x0(55:72);
    C_K = x0(73:90);
    C_MP = x0(91:108);
    C_AP = x0(109:126);
    
    %SIMO-Model equations
    S = x0(127:129);
    J = x0(130:132);
    R = x0(133:135);
    L = x0(136:138);
    
    %Insulin equations:
    P = x0(139);
    II = x0(140);
    QQ = x0(141);
    
    %Calculating the production rate vectors as well as functions used in
    %glucagon and insulin submodel
    [p,I,gamma] = Modelparameterfunctions(Km,Vm,C_B,C_H,C_G,C_L,C_K,...
    C_MP,C_AP,p,J,L,I,QQ,II,gamma,mu,Q,V);
    
    %SIMO-vectors used for macronutrient uptake
    %GUT; consisting of zeroes except for index of GLC and AA.
    SIMO_GI = zeros(18,1); SIMO_GI(1) = p.GI_GLC;
    SIMO_GI(10) = p.GI_AA;
    %Muscle: Consisting of zeroes except for index of TGL
    SIMO_MP = zeros(18,1); SIMO_MP(12) = p.GI_TGL * 0.5;
    %Adipose: Consisting of zeroes except for index of TGL
    SIMO_AP = zeros(18,1); SIMO_AP(12) = p.GI_TGL * 0.5;
    
    %The Differential equations of the main model:
    %Brain
    dC_Bdt = (Q.B * p.M * (C_H - C_B)) / V.B + p.R_B;
    %Heart
    dC_Hdt = (p.M * (Q.B * C_B + Q.L * C_L + Q.K * C_K + Q.MP * C_MP...
        + Q.AP * C_AP - Q.H * C_H)) / V.H + p.R_H;
    %Gut
    dC_Gdt = (Q.G * p.M * (C_H - C_G) + SIMO_GI) / V.G + p.R_G;
    %Liver
    dC_Ldt = (p.M * (Q.A * C_H + Q.G * C_G - Q.L * C_L)) / V.L + p.R_L;
    %Kidney
    dC_Kdt = (Q.K * p.M * (C_H - C_K)) / V.K + p.R_K;
    %Muscle periphery
    dC_MPdt = (Q.MP * p.M * (C_H - C_MP) + SIMO_MP) / V.MP + p.R_MP;
    %Adipose periphery
    dC_APdt = (Q.AP * p.M * (C_H - C_AP) + SIMO_AP) / V.AP + p.R_AP;
    %Return
    dCdt = [dC_Bdt; dC_Hdt; dC_Gdt; dC_Ldt; dC_Kdt; dC_MPdt; dC_APdt];
    
    %SIMO-Model:
    %Stomach
    dSdt = - p.GI_k_js .* S;
    %Jejenum
    dJdt = p.GI_k_js .* S - p.GI_k_gj .* J - p.GI_k_rj .* J;
    %Delay
    dRdt = p.GI_k_rj .* J - p.GI_k_lr .* R;
    %Ileum
    dLdt = p.GI_k_lr .* R - p.GI_k_gl .* L;
    
    %Insulin release equations
    dPdt = I.alpha * (I.P_inf - P);
    dIIdt = I.beta * (I.X - II);
    dQQdt = I.K * (I.Q_0 - QQ) + I.y * P - I.S;
    

    %Return vector
    dzdt = [dCdt;...
        %SIMO
        dSdt; dJdt; dRdt; dLdt;...
        %Insulin:
        dPdt; dIIdt; dQQdt];
    
    end