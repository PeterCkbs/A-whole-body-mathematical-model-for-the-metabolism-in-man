%% Driver script for the model used in the scientific journal
clear, clc

%Loading parameters
[Q,V,Km,Vm,p,I,gamma,mu] = LoadParmModel();

%Set basal values for glucose, insulin and glucagon
G0 = 5; %[mmol/L]
I0 = 15.1765; %[mU/L]
Gamma0 = 100; %[ng/L]

%Defining the rest of the basal values
[x0,I,gamma] = ModelBasalValues(G0,I0,Gamma0,I,gamma,Q);

%Food from t=0, [g], [GLC,AA,TGL]
Food = [0.6,0.24,0.16]'*100;
%Converting it to [mmol]
Conversion = Food .* [1/180*1000,1/89.1*1000,1/860*1000]';
x0(127:129) = Conversion;

%Time [min]
tspan = 0:1:60*72;

%Running the model
options = odeset("NonNegative",1:138);
[T,X] = ode15s(@Model,tspan,x0,options,Q,V,Km,Vm,p,I,gamma,mu);


%The following is only used for plots, it is the indexes of the metabolites
%It is needed to run the plot code for the figures in the article
GLC = 18*[0:6] + 1; G6P = 18*[0:6] + 2; GA3P = 18*[0:6] + 4;
CIT = 18*[0:6] + 8; LAC = 18*[0:6] + 9; AA = 18*[0:6] + 10;
FFA = 18*[0:6] + 11; TGL = 18*[0:6] + 12; GLR = 18*[0:6] + 13; 
KET = 18*[0:6] + 14; ACoA = 18*[0:6] + 6; GLY = 18*[0:6] + 3;
OXA = 18*[0:6] + 7; PYR = 18*[0:6] + 5;
PRO = 105; %Muscle tissue only
TGL_AP = 124; %Adipose tissue only
SIMO_GLC = [127,130,133,136];
SIMO_AA = [128,131,134,137];
SIMO_TGL = [129,132,135,138];
INS = 18*[0:6] + 17; GLU = 18*[0:6] + 18;

% From here, experimentation with plots is possible
% Plots from the article is included in the following sections

%% Generating Figure 7 from article

% First running model
[Q,V,Km,Vm,p,I,gamma,mu] = LoadParmModel();
G0 = 5; I0 = 15.1765; Gamma0 = 100;
[x0,I,gamma] = ModelBasalValues(G0,I0,Gamma0,I,gamma,Q);
Food = [0.6,0.24,0.16]'*100;
Conversion = Food .* [1/180*1000,1/89.1*1000,1/860*1000]';
x0(127:129) = Conversion;
tspan = 0:1:60*72;
options = odeset("NonNegative",1:138);
[T,X] = ode15s(@Model,tspan,x0,options,Q,V,Km,Vm,p,I,gamma,mu);

% Then generating Figure 7
figure()
axis tight
set(gca,'LooseInset',get(gca,'TightInset'));
sgtitle("Metabolite concentrations in Liver")
subplot(2,3,1);
plot(T(1:60*5)/60-5,repelem(X(1,GLC(4)),60*5),"LineWidth",1.5,"Color",[0.4940 0.1840 0.5560])
line(T/60,X(:,GLC(4)),"LineWidth",1.5,"Color",[0.4940 0.1840 0.5560]); grid on;
xline(0,"--",'LabelOrientation','horizontal');
xline(3,"--",'LabelOrientation','horizontal');
xline(18,"--",'LabelOrientation','horizontal');
xline(48,"--",'LabelOrientation','horizontal');
ylabel("GLC concentration [mmol/L]"); xlabel("Time [h]");
ylim([min(X(:,GLC(4)))*0.85,max(X(:,GLC(4)))*1.15])
xlim([-5,72])
legend("L","Location","best");

subplot(2,3,2);
plot(T(1:60*5)/60-5,repelem(X(1,AA(4)),60*5),"LineWidth",1.5,"Color",[0.4940 0.1840 0.5560])
line(T/60,X(:,AA(4)),"LineWidth",1.5,"Color",[0.4940 0.1840 0.5560]); grid on;
ylabel("AA concentration [mmol/L]"); xlabel("Time [h]");
ylim([min(X(:,AA(4)))*0.85,max(X(:,AA(4)))*1.15])
xline(0,"--",'LabelOrientation','horizontal');
xline(3,"--",'LabelOrientation','horizontal');
xline(18,"--",'LabelOrientation','horizontal');
xline(48,"--",'LabelOrientation','horizontal');
xlim([-5,72])

subplot(2,3,3);
plot(T(1:60*5)/60-5,repelem(X(1,TGL(4)),60*5),"LineWidth",1.5,"Color",[0.4940 0.1840 0.5560])
line(T/60,X(:,TGL(4)),"LineWidth",1.5,"Color",[0.4940 0.1840 0.5560]); grid on;
ylabel("TGL concentration [mmol/L]"); xlabel("Time [h]");
ylim([min(X(:,TGL(4)))*0.85,max(X(:,TGL(4)))*1.15])
xline(0,"--",'LabelOrientation','horizontal');
xline(3,"--",'LabelOrientation','horizontal');
xline(18,"--",'LabelOrientation','horizontal');
xline(48,"--",'LabelOrientation','horizontal');
xlim([-5,72])

subplot(2,3,4);
plot(T(1:60*5)/60-5,repelem(X(1,GLY(4)),60*5),"LineWidth",1.5,"Color",[0.4940 0.1840 0.5560])
line(T/60,X(:,GLY(4)),"LineWidth",1.5,'Color',[0.4940 0.1840 0.5560]); grid on;
ylabel("GLY concentration [mmol/L]"); xlabel("Time [h]");
ylim([min(X(:,GLY(4)))*0.85,max(X(:,GLY(4)))*1.15])
xline(0,"--",'LabelOrientation','horizontal');
xline(3,"--",'LabelOrientation','horizontal');
xline(18,"--",'LabelOrientation','horizontal');
xline(48,"--",'LabelOrientation','horizontal');
xlim([-5,72])

subplot(2,3,5);
plot(T(1:60*5)/60-5,repelem(X(1,FFA(4)),60*5),"LineWidth",1.5,"Color",[0.4940 0.1840 0.5560])
line(T/60,X(:,FFA(4)),"LineWidth",1.5,"Color",[0.4940 0.1840 0.5560]); grid on;
ylabel("FFA concentration [mmol/L]"); xlabel("Time [h]");
ylim([min(X(:,FFA(4)))*0.85,max(X(:,FFA(4)))*1.15])
xline(0,"--",'LabelOrientation','horizontal');
xline(3,"--",'LabelOrientation','horizontal');
xline(18,"--",'LabelOrientation','horizontal');
xline(48,"--",'LabelOrientation','horizontal');
xlim([-5,72])

subplot(2,3,6);
plot(T(1:60*5)/60-5,repelem(X(1,GLR(4)),60*5),"LineWidth",1.5,"Color",[0.4940 0.1840 0.5560])
line(T/60,X(:,GLR(4)),"LineWidth",1.5,"Color",[0.4940 0.1840 0.5560]); grid on;
ylabel("GLR concentration [mmol/L]"); xlabel("Time [h]");
ylim([min(X(:,GLR(4)))*0.85,max(X(:,GLR(4)))*1.15])
xline(0,"--",'LabelOrientation','horizontal');
xline(3,"--",'LabelOrientation','horizontal');
xline(18,"--",'LabelOrientation','horizontal');
xline(48,"--",'LabelOrientation','horizontal');
xlim([-5,72])

%% Generating Figure 8a and b from article

% First running model
[Q,V,Km,Vm,p,I,gamma,mu] = LoadParmModel();
G0 = 5; I0 = 15.1765; Gamma0 = 100;
[x0,I,gamma] = ModelBasalValues(G0,I0,Gamma0,I,gamma,Q);
Food = [0.6,0.24,0.16]'*100;
Conversion = Food .* [1/180*1000,1/89.1*1000,1/860*1000]';
x0(127:129) = Conversion;
tspan = 0:1:60*72;
options = odeset("NonNegative",1:138);
[T,X] = ode15s(@Model,tspan,x0,options,Q,V,Km,Vm,p,I,gamma,mu);

%Then generating fluxes using a loop
Flux_B = zeros(length(tspan),18); Flux_H = zeros(length(tspan),18);
Flux_G = zeros(length(tspan),18); Flux_L = zeros(length(tspan),18);
Flux_K = zeros(length(tspan),18); Flux_MP = zeros(length(tspan),18);
Flux_AP = zeros(length(tspan),18);

for i=1:length(tspan)
    C_B = X(i,1:18); 
    C_H = X(i,19:36); 
    C_G = X(i,37:54); 
    C_L = X(i,55:72);
    C_K = X(i,73:90); 
    C_MP = X(i,91:108); 
    C_AP = X(i,109:126); 
    
    J = X(i,130:132);
    L = X(i,136:138); 
    
    II = X(i,140); 
    QQ = X(i,141); 
    Gamma = X(i,135);
    
    [p,I,gamma] = Modelparameterfunctions(Km,Vm,C_B,C_H,C_G,C_L,C_K,...
    C_MP,C_AP,p,J,L,I,QQ,II,gamma,mu,Q,V);
    Flux_B(i,:) = p.R_B*V.B; Flux_H(i,:) = p.R_H*V.B;
    Flux_G(i,:) = p.R_G*V.G; Flux_L(i,:) = p.R_L*V.L;
    Flux_K(i,:) = p.R_K*V.K; Flux_MP(i,:) = p.R_MP*V.MP;
    Flux_AP(i,:) = p.R_AP*V.AP;
end

%Then define vectors of non-zero fluxes

%Brain fluxes:
B_ind = [1,2,4,5,6,7,8,9,14]; %are non-zero
%Heart fluxes:
H_ind = [1,2,4,5,6,7,8,9,11,12,13,14]; %are non-zero
%Gut fluxes:
G_ind = [1,2,4,5,6,7,8,9,10,11]; %are non-zero
%Liver fluxes:
L_ind = [1,2,3,4,5,6,7,8,9,10,11,12,13,14]; %are non-zero
%Kidney fluxes:
K_ind = [1,2,4,5,6,7,8,9,10,11,13]; %are non-zero
%Muscle fluxes:
MP_ind = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]; %are non-zero
%Adipose fluxes:
AP_ind = [1,2,4,5,6,7,8,9,11,12,13,16]; %are non-zero

%Total fluxes:
Flux_total = [Flux_B,Flux_H,Flux_G,Flux_L,Flux_K,Flux_MP,Flux_AP];

% Figure 8b from article
figure();
semilogx(tspan/60,Flux_total(:,GLC),"LineWidth",1.5); grid on; title("Glucose fluxes");
legend(["B","H","G","L","K","MP","AP"],"Location","southeast");
ylabel("[mmol/min]"); xlabel("Time [h]");

%Summin up the fluxes, so it can be used in a barplot
sum_B = sum(Flux_B); sum_H = sum(Flux_H); sum_G = sum(Flux_G); sum_L = sum(Flux_L);
sum_K = sum(Flux_K); sum_MP = sum(Flux_MP); sum_AP = sum(Flux_AP);
Flux_mat = [sum_B;sum_H;sum_G;sum_L;sum_K;sum_MP;sum_AP]';
%And the sum of the fluxes
Fluxsum = sum_B + sum_H + sum_G + sum_L + sum_K + sum_MP + sum_AP;

% Figure 8a from article
figure(); 
title("Sum of Glucose fluxes"); hold on; grid on
bar(1,Flux_mat(1,1)/1000*180, 'FaceColor', [0 0.4470 0.7410])
    bar(2,Flux_mat(1,2)/1000*180, 'FaceColor',[0.8500 0.3250 0.0980])
    bar(3,Flux_mat(1,3)/1000*180, 'FaceColor',[0.9290 0.6940 0.1250])
    bar(4,Flux_mat(1,4)/1000*180, 'FaceColor',[0.4940 0.1840 0.5560])
    bar(5,Flux_mat(1,5)/1000*180, 'FaceColor',[0.4660 0.6740 0.1880])
    bar(6,Flux_mat(1,6)/1000*180, 'FaceColor',[0.3010 0.7450 0.9330])
    bar(7,Flux_mat(1,7)/1000*180, 'FaceColor',[0.6350 0.0780 0.1840])
    bar(8,Fluxsum(1)/1000*180, "b")
    set(gca,"Xtick",1:8,'xticklabel',["B","H","G","L","K","MP","AP","Sum"])
ylabel("g");
hold off;


%% Generating Figure 9 from article

% First running model, but now with multiple meals
options = odeset("NonNegative",1:138);
[Q,V,Km,Vm,p,I,gamma,mu] = LoadParmModel();
G0 = 5; I0 = 15.1765; Gamma0 = 100;
[x0,I,gamma] = ModelBasalValues(G0,I0,Gamma0,I,gamma,Q);
%Defining the meal-sizes
Food = [0.6,0.24,0.16]'*100;
Conversion = Food .* [1/180*1000,1/89.1*1000,1/860*1000]';

%Defining what time points to eat
timings = [1,6,11,16,25,30,35,40,49,54,59,64,72]*60;

%Loop for eating multiple times a day
X = x0;
for i=1:length(timings)
    tspan = 0:1:timings(1);
    if i ~= 1
        x0(127:129) = x0(127:129) + Conversion;
        tspan = timings(i-1):1:timings(i);
    end

    [T,X] = ode15s(@Model,tspan,x0,options,Q,V,Km,Vm,p,I,gamma,mu);
    x0 = X(end,:)';
    if i==1
        X1 = X;
        T1 = length(T)-1;
    else
        X1 = [X1;X(2:end,:)];
        T1 = T1 + length(T)-1;
    end
end
X = X1;
T = 0:1:T1;

% Then generating Figure 9
figure()
axis tight
set(gca,'LooseInset',get(gca,'TightInset'));
sgtitle("Metabolite concentrations in Liver")
subplot(2,3,1);
plot(T(1:60*5)/60-5,repelem(X(1,GLC(4)),60*5),"LineWidth",1.5,"Color",[0.4940 0.1840 0.5560])
line(T/60,X(:,GLC(4)),"LineWidth",1.5,"Color",[0.4940 0.1840 0.5560]); grid on;
xline(0,"--",'LabelOrientation','horizontal');
xline(3,"--",'LabelOrientation','horizontal');
xline(18,"--",'LabelOrientation','horizontal');
xline(48,"--",'LabelOrientation','horizontal');
ylabel("GLC concentration [mmol/L]"); xlabel("Time [h]");
ylim([min(X(:,GLC(4)))*0.85,max(X(:,GLC(4)))*1.15])
xlim([-5,72])
legend("L","Location","best");

subplot(2,3,2);
plot(T(1:60*5)/60-5,repelem(X(1,AA(4)),60*5),"LineWidth",1.5,"Color",[0.4940 0.1840 0.5560])
line(T/60,X(:,AA(4)),"LineWidth",1.5,"Color",[0.4940 0.1840 0.5560]); grid on;
ylabel("AA concentration [mmol/L]"); xlabel("Time [h]");
ylim([min(X(:,AA(4)))*0.85,max(X(:,AA(4)))*1.15])
xline(0,"--",'LabelOrientation','horizontal');
xline(3,"--",'LabelOrientation','horizontal');
xline(18,"--",'LabelOrientation','horizontal');
xline(48,"--",'LabelOrientation','horizontal');
xlim([-5,72])

subplot(2,3,3);
plot(T(1:60*5)/60-5,repelem(X(1,TGL(4)),60*5),"LineWidth",1.5,"Color",[0.4940 0.1840 0.5560])
line(T/60,X(:,TGL(4)),"LineWidth",1.5,"Color",[0.4940 0.1840 0.5560]); grid on;
ylabel("TGL concentration [mmol/L]"); xlabel("Time [h]");
ylim([min(X(:,TGL(4)))*0.85,max(X(:,TGL(4)))*1.15])
xline(0,"--",'LabelOrientation','horizontal');
xline(3,"--",'LabelOrientation','horizontal');
xline(18,"--",'LabelOrientation','horizontal');
xline(48,"--",'LabelOrientation','horizontal');
xlim([-5,72])

subplot(2,3,4);
plot(T(1:60*5)/60-5,repelem(X(1,GLY(4)),60*5),"LineWidth",1.5,"Color",[0.4940 0.1840 0.5560])
line(T/60,X(:,GLY(4)),"LineWidth",1.5,'Color',[0.4940 0.1840 0.5560]); grid on;
ylabel("GLY concentration [mmol/L]"); xlabel("Time [h]");
ylim([min(X(:,GLY(4)))*0.85,max(X(:,GLY(4)))*1.15])
xline(0,"--",'LabelOrientation','horizontal');
xline(3,"--",'LabelOrientation','horizontal');
xline(18,"--",'LabelOrientation','horizontal');
xline(48,"--",'LabelOrientation','horizontal');
xlim([-5,72])

subplot(2,3,5);
plot(T(1:60*5)/60-5,repelem(X(1,FFA(4)),60*5),"LineWidth",1.5,"Color",[0.4940 0.1840 0.5560])
line(T/60,X(:,FFA(4)),"LineWidth",1.5,"Color",[0.4940 0.1840 0.5560]); grid on;
ylabel("FFA concentration [mmol/L]"); xlabel("Time [h]");
ylim([min(X(:,FFA(4)))*0.85,max(X(:,FFA(4)))*1.15])
xline(0,"--",'LabelOrientation','horizontal');
xline(3,"--",'LabelOrientation','horizontal');
xline(18,"--",'LabelOrientation','horizontal');
xline(48,"--",'LabelOrientation','horizontal');
xlim([-5,72])

subplot(2,3,6);
plot(T(1:60*5)/60-5,repelem(X(1,GLR(4)),60*5),"LineWidth",1.5,"Color",[0.4940 0.1840 0.5560])
line(T/60,X(:,GLR(4)),"LineWidth",1.5,"Color",[0.4940 0.1840 0.5560]); grid on;
ylabel("GLR concentration [mmol/L]"); xlabel("Time [h]");
ylim([min(X(:,GLR(4)))*0.85,max(X(:,GLR(4)))*1.15])
xline(0,"--",'LabelOrientation','horizontal');
xline(3,"--",'LabelOrientation','horizontal');
xline(18,"--",'LabelOrientation','horizontal');
xline(48,"--",'LabelOrientation','horizontal');
xlim([-5,72])

%% Generating Figure 10 from article

% First running model, but now only eating every other day
options = odeset("NonNegative",1:138);
[Q,V,Km,Vm,p,I,gamma,mu] = LoadParmModel();
G0 = 5; I0 = 15.1765; Gamma0 = 100;
[x0,I,gamma] = ModelBasalValues(G0,I0,Gamma0,I,gamma,Q);
%Defining the meal-sizes
Food = [0.6,0.24,0.16]'*100;
Conversion = Food .* [1/180*1000,1/89.1*1000,1/860*1000]';

%Defining what time points to eat (regular eating every other day)
days = 13;
timings = [];
for i=0:2:days
    t = [1+(24*i),6+(24*i),11+(24*i),16+(24*i)]*60;
    timings = [timings, t];
end

%Loop for eating multiple times a day
X = x0;
for i=1:length(timings)
    tspan = 0:1:timings(1);
    if i ~= 1
        x0(127:129) = x0(127:129) + Conversion;
        tspan = timings(i-1):1:timings(i);
    end

    [T,X] = ode15s(@Model,tspan,x0,options,Q,V,Km,Vm,p,I,gamma,mu);
    x0 = X(end,:)';
    if i==1
        X1 = X;
        T1 = length(T)-1;
    else
        X1 = [X1;X(2:end,:)];
        T1 = T1 + length(T)-1;
    end
end
X = X1;
T = 0:1:T1;

% Generating figure 10
figure()
sgtitle("Metabolite concentrations")
subplot(2,3,1);
plot(T/(60*24),X(:,TGL_AP),"LineWidth",2,"Color",[0.6350 0.0780 0.1840]); grid on;
ylabel("TGL_{AP} concentration [mmol/L]"); xlabel("Time [days]");
ylim([min(X(:,TGL_AP))*0.85,max(X(:,TGL_AP))*1.15])
legend("TGL_{AP}","Location","northeast");
xlim([0,13]);

subplot(2,3,2);
plot(T/(60*24),X(:,PRO),"LineWidth",2,"Color",[0.3010 0.7450 0.9330]); grid on;
ylabel("PRO concentration [mmol/L]"); xlabel("Time [days]");
ylim([min(X(:,PRO))*0.85,max(X(:,PRO))*1.15])
legend("PRO","Location","northeast");
xlim([0,13]);

subplot(2,3,3);
p = plot(T/(60*24),X(:,GLY(4)),T/(60*24),X(:,GLY(6)));
p(1).LineWidth = 1.5;
p(2).LineWidth = 1.5;
p(1).Color = [0.4940 0.1840 0.5560];
p(2).Color = [0.3010 0.7450 0.9330];
grid on;
ylabel("GLY concentration [mmol/L]"); xlabel("Time [days]");
legend("GLY_{L}","GLY_{MP}","Location","northeast");
xlim([0,13]);

subplot(2,3,4);
plot(T/(60*24),X(:,TGL(4)),"LineWidth",1.5,"Color",[0.4940 0.1840 0.5560]); grid on;
ylabel("TGL concentration [mmol/L]"); xlabel("Time [days]");
ylim([min(X(:,TGL(4)))*0.85,max(X(:,TGL(4)))*1.15])
legend("TGL_L","Location","northeast");
xlim([0,13]);

subplot(2,3,5);
plot(T/(60*24),X(:,FFA(4)),"LineWidth",1.5,"Color",[0.4940 0.1840 0.5560]); grid on;
ylabel("FFA concentration [mmol/L]"); xlabel("Time [days]");
ylim([min(X(:,FFA(4)))*0.85,max(X(:,FFA(4)))*1.15])
legend("FFA_L","Location","northeast");
xlim([0,13]);

subplot(2,3,6);
plot(T/(60*24),X(:,GLR(4)),"LineWidth",1.5,"Color",[0.4940 0.1840 0.5560]); grid on;
ylabel("GLR concentration [mmol/L]"); xlabel("Time [days]");
ylim([min(X(:,GLR(4)))*0.85,max(X(:,GLR(4)))*1.15])
legend("GLR_L","Location","northeast");
xlim([0,13]);
