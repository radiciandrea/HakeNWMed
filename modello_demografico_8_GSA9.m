% Demographic model of Radici et al., 2021
% Simulation of a single - cell protection

if exist('Q') == 0
    clear 
    close
    clc
    Q = [3.3957 0.0004 6.6525 0.0001 164.5678];
    modello = 1; %0 = lineare, 1 = ricker, 2 = BH
    kernel = 3;  %1 = pl; 2 = exp; 3 = gau
end


if exist('c_p') == 0
    c_p = 500;
end
%% IMPOSTAZIONI DA SETTARE

T = 50; % periodo della simulazione
Tmax = 100;

addpath data
load n_celle_regioni_3
load celle_M_4_3_o
load Mc_4_3_o_perc
load R10_EFF9_var

%% INIZIALIZZAZIONE PARAMETRI

limiti_regioni = n_celle_regioni_3;
celle_o = celle_M_4_3_o(limiti_regioni(4,1):end);
celle_tot = celle_M_4_3_o;

M =  Mc_4_3_o_perc;

classi_eta = 0:6;
T = min(T,Tmax);
reg = {'Lazio', 'Toscana', 'Liguria'};

% mortalità istantanea naturale (y^-1)
M_9 = [1.2 0.62	0.44 0.37 0.33 0.31 0.29];

% maturità (frazione)
m_9 = [0 0.25 0.9 1 1 1 1];

% peso medio (kg)
p_9 = [0.008 0.166 0.578 1.2 1.949 2.745 3.529 ];

% sr relazione
% Ricker per la GSA 9 (liguria, toscana, lazio)
alfa_r =  56.0398269589289;
beta_r = 0.000740805842994942;

% lineare
alfa_l = 30.8682878331365;
beta_l = 0;

% BH
alfa_bh = 61.9782853313202;
beta_bh =  0.00125892356243384;

% contributo esterno
C_10 = zeros(T,1) + R_10_v(1:T);

load persistenza_n
load persistenza_sg

%% coefficenti per la pesca E protezione delle celle

coefficienti_fisheries

temp_otb_012 = sum(P_otb_012.*A_f);
temp_otb_1224 = sum(P_otb_1224.*A_f);
temp_ssc_012 = sum(P_ssc_012.*A_f);
temp_ssc_1224 = sum(P_ssc_1224.*A_f);

% si protegge la cella indicata

P_otb_012(id_f == c_p) = 0;
P_otb_1224(id_f == c_p) = 0;
P_ssc_012(id_f == c_p) = 0;
P_ssc_1224(id_f == c_p) = 0;

% si redistribuiscele lo sforzo su tutte le altre celle
P_otb_012 = P_otb_012*temp_otb_012/(sum(P_otb_012.*A_f));
P_otb_1224 = P_otb_1224*temp_otb_1224/(sum(P_otb_1224.*A_f));
P_ssc_012 = P_ssc_012*temp_ssc_012/(sum(P_ssc_012.*A_f));
P_ssc_1224 = P_ssc_1224*temp_ssc_1224/(sum(P_ssc_1224.*A_f));

eff_gns = GNS_9_v(1:T);
eff_otb = OTB_9_v(1:T);

%% INIZIALIZZAZIONE STATO INIZIALE
% questi coefficienti considerano già la cella intera.

%load mare_o
load area_celle_o

area_cella_deg = 0.023438;
area_cella = area_celle_o; %.*mare_o; %km^2

N_9 = [88907 9663.6 1293 297.02 136.8 139.8 33.381]*1000; % 2014
%N_9 = [84450 11107 2902 360.77 94.19 2.944 1.842]*1000; % 2006

% Ogni riga è una diversa cella, ogni colonna una classe d'età
N = nan(length(celle_o), length(classi_eta), T);

for i = 1: length(N_9) 
    N(:,i,1)= m_9(i)*repmat(N_9(i), (limiti_regioni(end) - limiti_regioni(3,2)), 1).*...
    persistenza_sg(limiti_regioni(4,1):end)./...
    repmat(sum(persistenza_sg(limiti_regioni(4,1):limiti_regioni(end))),...
    (limiti_regioni(end) - limiti_regioni(3,2)), 1)+...
    (1-m_9(i))*repmat(N_9(i), (limiti_regioni(end) - limiti_regioni(3,2)), 1).*...
    persistenza_n(limiti_regioni(4,1):end)./...
    repmat(sum(persistenza_n(limiti_regioni(4,1):limiti_regioni(end))),...
    (limiti_regioni(end) - limiti_regioni(3,2)), 1);
end

SSB = nan(length(celle_o), T-1);
L = nan(length(celle_o), T-1);
catches_n = nan(length(celle_o), length(classi_eta), T-1); % quantità assolute
deaths_n = nan(length(celle_o), length(classi_eta), T-1);

landings_gns = nan(T-1,1);
landings_otb = nan(T-1,1);
%% CICLO

for t = 1:T-1
   
    % calcolo da densità di SSB
    
    SSB0_cl = N(:,:,t)./(repmat(persistenza_sg(limiti_regioni(4,1):limiti_regioni(end)).*area_cella(limiti_regioni(4,1):limiti_regioni(end)),1 , length(classi_eta)));
    SSB0_cl(isnan(SSB0_cl)) = 0;
    SSB0_cl(SSB0_cl == Inf) = 0;
    SSB0_cl(SSB0_cl == -Inf) = 0;
    
    SSB0 = sum(repmat(m_9, (limiti_regioni(end) - limiti_regioni(3,2)), 1 ).*...
        repmat(p_9, (limiti_regioni(end) - limiti_regioni(3,2)), 1 ).*...
        SSB0_cl, 2);
    
    % si applica la relazione S/R (
    R0 = (modello == 0)*repmat(alfa_r, (limiti_regioni(end) - limiti_regioni(3,2)), 1).*...
        SSB0 + ...
        (modello == 1)*repmat(alfa_r, (limiti_regioni(end) - limiti_regioni(3,2)), 1).*...
        SSB0.*exp(-SSB0.*...
        repmat(beta_r, (limiti_regioni(end) - limiti_regioni(3,2)), 1)) + ...
        (modello == 2)*repmat(alfa_bh, (limiti_regioni(end) - limiti_regioni(3,2)), 1).*...
        SSB0./(1 + SSB0.*...
        repmat(beta_bh, (limiti_regioni(end) - limiti_regioni(3,2)), 1));
    
    R_9 = (R0).*(area_cella(limiti_regioni(4,1):limiti_regioni(end)).*persistenza_n(limiti_regioni(4,1):end));
    
    % si distribuiscono i recruits secondo la matrice di connettività.
    % c'è un fattore correttivo. Si incolonnano anche i recruits della GSA
    % 10:
    R_10 = C_10(t)*persistenza_sg(1:limiti_regioni(3,2))/sum(persistenza_sg(1:limiti_regioni(3,2)));
    
    R_9_10 = [R_10; R_9];
    
    R_9_10_C = (R_9_10') * M * (sum(R_9_10) / sum(sum((R_9_10') * M)));
    
    
    % Si calcolano i Recruits effettivi, e non la densità. In questo caso,
    % prima o dopo aver fatto
    
    R = R_9_10_C(limiti_regioni(4,1):end);
    
    % Calcolo tutte le classi dell'anno successivo: recruits
    N(:,1,t+1) = R;
    
    % spostamento degli adulti
    for i = 2:7
        N(:,i,t) = N(:,i,t)' * D;
    end
    
    n_temp = N(cella_temp,:,t).*... 
        repmat(A_f, 1, length(classi_eta))./...
        repmat(area_pesca(cella_temp), 1, length(classi_eta));

    
    % considero 2 volte otb 
    F_temp = reshape(sum(...
        repmat([cat_gns_625; cat_gns_625; cat_gns_820; cat_gns_820; cat_otb_33; cat_otb_40], 1,1,N_max).*...
        repmat([eff_gns(t); eff_gns(t); eff_gns(t); eff_gns(t); eff_otb(t); eff_otb(t)], 1, length(classi_eta),N_max)...
        .*repmat(reshape([P_ssc_012 P_ssc_1224 P_ssc_012 P_ssc_1224 P_otb_1224 P_otb_1224]', 6,1, N_max)...
        , 1, length(classi_eta),1)), length(classi_eta), N_max)';
    
    F_gns = reshape(sum(...
        repmat([cat_gns_625; cat_gns_625; cat_gns_820; cat_gns_820; zeros(size(cat_otb_33)); zeros(size(cat_otb_40))], 1,1,N_max).*...
        repmat([eff_gns(t); eff_gns(t); eff_gns(t); eff_gns(t); eff_otb(t); eff_otb(t)], 1, length(classi_eta),N_max)...
        .*repmat(reshape([P_ssc_012 P_ssc_1224 P_ssc_012 P_ssc_1224 P_otb_1224 P_otb_1224]', 6,1, N_max)...
        , 1, length(classi_eta),1)), length(classi_eta), N_max)';
    
    
    F_otb = reshape(sum(...
        repmat([zeros(size(cat_gns_625)); zeros(size(cat_gns_625)); zeros(size(cat_gns_820)); zeros(size(cat_gns_820)); cat_otb_33; cat_otb_40], 1,1,N_max).*...
        repmat([eff_gns(t); eff_gns(t); eff_gns(t); eff_gns(t); eff_otb(t); eff_otb(t)], 1, length(classi_eta),N_max)...
        .*repmat(reshape([P_ssc_012 P_ssc_1224 P_ssc_012 P_ssc_1224 P_otb_1224 P_otb_1224]', 6,1, N_max)...
        , 1, length(classi_eta),1)), length(classi_eta), N_max)';
    
    % Eq Baranov
    catches_temp = F_temp./(F_temp + repmat(M_9, N_max,1)).*...
        (1-exp(-(F_temp + repmat(M_9, N_max,1)))).*n_temp;  

    deaths_temp = (1-exp(-(F_temp + repmat(M_9, N_max,1)))).*n_temp; 
    
    catches_temp_gns = F_gns./(F_temp + repmat(M_9, N_max,1)).*...
        (1-exp(-(F_temp + repmat(M_9, N_max,1)))).*n_temp.*repmat(p_9, N_max, 1);
        
    catches_temp_otb = F_otb./(F_temp + repmat(M_9, N_max,1)).*...
        (1-exp(-(F_temp + repmat(M_9, N_max,1)))).*n_temp.*repmat(p_9, N_max, 1);   
    
    for i = 1:length(celle_o)
        catches_n(i, :, t) = reshape(sum(catches_temp(id_f==celle_o(i),:),1),1, length(classi_eta),1);
        deaths_n(i, :, t) = reshape(sum(deaths_temp(id_f==celle_o(i),:),1),1, length(classi_eta),1);
        L(i, t) = sum(sum(catches_temp_gns(id_f==celle_o(i),2:end))) + ...
            sum(sum(catches_temp_otb(id_f==celle_o(i),2:end)));
    end
    
    %(tonnellate)
    landings_gns(t) = sum(sum(catches_temp_gns(:,2:end)))*0.001; % solo il peso è considerato, e solo individui > 1 anno
    landings_otb(t) = sum(sum(catches_temp_otb(:,2:end)))*0.001;
    
    % classi 2-5
    for i = 2: length(classi_eta)-1
        N(:,i,t+1) = N(:,i-1,t) - deaths_n(:,i-1,t);
    end
    
    % classe 6+
    N(:,length(classi_eta),t+1) = N(:,length(classi_eta),t) + N(:,length(classi_eta)-1,t)...
        - deaths_n(:,length(classi_eta),t) - deaths_n(:,length(classi_eta)-1,t);
        
    % Calcolo gli indicatori (kg)
    SSB(:,t) =sum(repmat(m_9, (limiti_regioni(end) - limiti_regioni(3,2)), 1 )...
        .*repmat(p_9, (limiti_regioni(end) - limiti_regioni(3,2)), 1 )...
        .*N(:,:,t), 2);
    
end

for i = 2:7
    N(:,i,end) = N(:,i,end)' * D;
end

NT =reshape(sum(N), length(classi_eta),T);

% PLOTTING (uno dei due)

landings = landings_gns + landings_otb;

N_tot = reshape(sum(N, 1), 7, T);
SSB_tot = reshape(sum(SSB, 1), T-1, 1);
L_tot = landings;

plot_modello_demografico_8_GSA_9
