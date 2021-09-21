%% PLOTTING

T_sim = 2014:2014+T-1;
T_sim_l = 2014:2014+T-2;

% normale
figure;
for i = 1:length(classi_eta)
    subplot(2,4,i)
    hold on
    plot(T_sim, N_tot(i,:));
    title(strcat('classe (', num2str(classi_eta(i)), ')'));
end

suptitle('GSA 9');


% landings and ssb
figure;
subplot(2,1,1)
hold on
plot(T_sim_l, L_tot);
title('Landings');

subplot(2,1, 2)
hold on
plot(T_sim_l, SSB_tot);
title('SSB');

%Landings by type
figure;
subplot(2,1,1)
hold on
plot(T_sim_l, landings_gns);
title('Landings by GNS');

subplot(2,1, 2)
hold on
plot(T_sim_l, landings_otb);
title('Landings by OTB');