%% Numerical Analysis: SDOF System with Coulomb Friction
%
% This set of codes provide the analysis of results for the Single 
% Degree-of-Freedom Dynamical System subjected to Coulomb Friction Force 
% that is presented in the literature:
%
% L. Marino and A. Cicirello (2020). Experimental investigation of a single-
% degree-of-freedom system with Coulomb friction. Nonlinear Dynamics, 99(3), 
% 1781-1799. doi: 10.1007/s11071-019-05443-2
%
% This set-up is applicable only for the Base motion (with fixed wall) case.
%-------------------------------------------------------------------------%
%% Load the data from Part 1:
load('experimental_SDOF_System_Coulomb_Friction_Part1_m')

%% Analysis of the model updating:
t_i = [1;2;3;4]; % Time-step values

posterior_mean_semc1 = zeros(length(t_i),size(semc_allsamples1,2));
posterior_mean_semc2 = zeros(length(t_i),size(semc_allsamples2,2));
posterior_std_semc1 = zeros(length(t_i),size(semc_allsamples1,2));
posterior_std_semc2 = zeros(length(t_i),size(semc_allsamples2,2));

posterior_bounds_coulomb_semc1 = zeros(length(t_i),2);
posterior_bounds_coulomb_semc2 = zeros(length(t_i),2);
posterior_bounds_omega_semc1 = zeros(length(t_i),2);
posterior_bounds_omega_semc2 = zeros(length(t_i),2);
posterior_bounds_sigma_phi_semc1 = zeros(length(t_i),2);
posterior_bounds_sigma_phi_semc2 = zeros(length(t_i),2);
posterior_bounds_sigma_r_semc1 = zeros(length(t_i),2);
posterior_bounds_sigma_r_semc2 = zeros(length(t_i),2);

for idx = 1:length(t_i)
posterior_mean_semc1(idx,1) = mean(semc_allsamples1(:,1,idx+1));
posterior_mean_semc1(idx,2) = mean(semc_allsamples1(:,2,idx+1));
posterior_mean_semc1(idx,3) = mean(semc_allsamples1(:,3,idx+1));
posterior_mean_semc1(idx,4) = mean(semc_allsamples1(:,4,idx+1));

posterior_mean_semc2(idx,1) = mean(semc_allsamples2(:,1,idx+1));
posterior_mean_semc2(idx,2) = mean(semc_allsamples2(:,2,idx+1));
posterior_mean_semc2(idx,3) = mean(semc_allsamples2(:,3,idx+1));
posterior_mean_semc2(idx,4) = mean(semc_allsamples2(:,4,idx+1));

posterior_std_semc1(idx,1) = std(semc_allsamples1(:,1,idx+1));
posterior_std_semc1(idx,2) = std(semc_allsamples1(:,2,idx+1));
posterior_std_semc1(idx,3) = std(semc_allsamples1(:,3,idx+1));
posterior_std_semc1(idx,4) = std(semc_allsamples1(:,4,idx+1));

posterior_std_semc2(idx,1) = std(semc_allsamples2(:,1,idx+1));
posterior_std_semc2(idx,2) = std(semc_allsamples2(:,2,idx+1));
posterior_std_semc2(idx,3) = std(semc_allsamples2(:,3,idx+1));
posterior_std_semc2(idx,4) = std(semc_allsamples2(:,4,idx+1));

posterior_bounds_coulomb_semc1(idx,:) = prctile(semc_allsamples1(:,1,idx+1), [5, 95]);
posterior_bounds_omega_semc1(idx,:) = prctile(semc_allsamples1(:,2,idx+1), [5, 95]);
posterior_bounds_sigma_phi_semc1(idx,:) = prctile(semc_allsamples1(:,3,idx+1), [5, 95]);
posterior_bounds_sigma_r_semc1(idx,:) = prctile(semc_allsamples1(:,4,idx+1), [5, 95]);

posterior_bounds_coulomb_semc2(idx,:) = prctile(semc_allsamples2(:,1,idx+1), [5, 95]);
posterior_bounds_omega_semc2(idx,:) = prctile(semc_allsamples2(:,2,idx+1), [5, 95]);
posterior_bounds_sigma_phi_semc2(idx,:) = prctile(semc_allsamples2(:,3,idx+1), [5, 95]);
posterior_bounds_sigma_r_semc2(idx,:) = prctile(semc_allsamples2(:,4,idx+1), [5, 95]);
end
posterior_cov_semc1 = (posterior_std_semc1./posterior_mean_semc1).*100;
posterior_cov_semc2 = (posterior_std_semc2./posterior_mean_semc2).*100;

%% Plotting the graphical results:

% To plot the estimates of F(t_i) with the theoretical model (Model 1 vs Model 2):
figure;
% Subplot for Model 1:
subplot(1,2,1)
hold on; grid on; box on;
plot(t_i, coulomb_force_nom, 'k--', 'linewidth', 1)
plot(t_i, coulomb_force, 'x',  'color', '#7E2F8E', 'linewidth', 2)
y_neg_a0 = abs(posterior_mean_semc1(:,1) - posterior_bounds_coulomb_semc1(:,1)); % error in the negative y-direction
y_pos_a0 = abs(posterior_mean_semc1(:,1) - posterior_bounds_coulomb_semc1(:,2)); % error in the positive y-direction
errorbar(t_i, posterior_mean_semc1(:,1), y_neg_a0, y_pos_a0, '-sb','MarkerSize',5,...
    'MarkerEdgeColor','blue','MarkerFaceColor','blue', 'linewidth',1);
legend('Markov model T_1','True values','SEMC estimates','linewidth', 2)
xlim([1 4]); ylim([0.10 1.6]);
xlabel('$t_{s}$ $[mth]$','Interpreter','latex')
ylabel('$F_{\mu}(t_{s})$ $[N]$','Interpreter','latex')
set(gca, 'fontsize', 20)

% Subplot for Model 2:
subplot(1,2,2)
hold on; grid on; box on;
plot((1:0.01:4), force_evolution2((1:0.01:4)), 'k--', 'linewidth', 1)
plot(t_i, coulomb_force, 'x',  'color', '#7E2F8E', 'linewidth', 2)
y_neg_b0 = abs(posterior_mean_semc2(:,1) - posterior_bounds_coulomb_semc2(:,1)); % error in the negative y-direction
y_pos_b0 = abs(posterior_mean_semc2(:,1) - posterior_bounds_coulomb_semc2(:,2)); % error in the positive y-direction
errorbar(t_i, posterior_mean_semc2(:,1), y_neg_b0, y_pos_b0, '-sm','MarkerSize',5,...
    'MarkerEdgeColor','magenta','MarkerFaceColor','magenta', 'linewidth',1);
legend('Markov model T_2','True values','SEMC estimates','linewidth', 2)
xlim([1 4]); ylim([0.10 1.6]);
xlabel('$t_{s}$ $[mth]$','Interpreter','latex')
ylabel('$F_{\mu}(t_{s})$ $[N]$','Interpreter','latex')
set(gca, 'fontsize', 20)

%% To plot the estimates of omega, sigma_v, and sigma_phi with the theoretical model (Model 1 vs Model 2):

% Compute reference value for sigma_phi:
rmse_phi = zeros(length(coulomb_force),1);
for idx = 1:length(coulomb_force)
model_output = blackbox_model(coulomb_force(idx), r_exp(:,idx), driving_force);
phi_mod = model_output.phase_angles;
rmse_phi(idx) = sqrt(mean((phi_exp(:,idx) - phi_mod).^2));   
end
sigma_phi = mean(rmse_phi);

% Compute reference value for sigma_r:
rmse_r = zeros(length(coulomb_force),1);
for idx = 1:length(coulomb_force)
rmse_r(idx) = sqrt(mean((r_exp(:,idx) - r_nom).^2)); 
end
sigma_r = mean(rmse_r);

figure;
subplot(2,2,1)
hold on; grid on; box on;
plot([1, 4], [omega_n, omega_n], 'k--', 'linewidth', 1)
y_neg_c1 = abs(posterior_mean_semc1(:,2) - posterior_bounds_omega_semc1(:,1)); % error in the negative y-direction
y_pos_c1 = abs(posterior_mean_semc1(:,2) - posterior_bounds_omega_semc1(:,2)); % error in the positive y-direction
errorbar(t_i, posterior_mean_semc1(:,2), y_neg_c1, y_pos_c1, '-sb','MarkerSize',5,...
    'MarkerEdgeColor','blue','MarkerFaceColor','blue','linewidth',1);
y_neg_d1 = abs(posterior_mean_semc2(:,2) - posterior_bounds_omega_semc2(:,1)); % error in the negative y-direction
y_pos_d1 = abs(posterior_mean_semc2(:,2) - posterior_bounds_omega_semc2(:,2)); % error in the positive y-direction
errorbar(t_i, posterior_mean_semc2(:,2), y_neg_d1, y_pos_d1, '-sm','MarkerSize',5,...
    'MarkerEdgeColor','magenta','MarkerFaceColor','magenta','linewidth',1);
legend(['True value = ',num2str(omega_n, '%.3f'), ' rad/s'],'SEMC estimates, T_{1}',...
       'SEMC estimates, T_{2}','linewidth', 2)
xlim([1 4]); ylim([0 100]);
xlabel('$t_{s}$ $[mth]$','Interpreter','latex')
ylabel('$\omega_{n}$ $[rad/s]$','Interpreter','latex')
set(gca, 'fontsize', 20)

subplot(2,2,2)
hold on; grid on; box on;
plot([1, 4], [sigma_phi, sigma_phi], 'k--', 'linewidth', 1)
y_neg_c2 = abs(posterior_mean_semc1(:,3) - posterior_bounds_sigma_phi_semc1(:,1)); % error in the negative y-direction
y_pos_c2 = abs(posterior_mean_semc1(:,3) - posterior_bounds_sigma_phi_semc1(:,2)); % error in the positive y-direction
errorbar(t_i, posterior_mean_semc1(:,3), y_neg_c2, y_pos_c2, '-sb','MarkerSize',5,...
    'MarkerEdgeColor','blue','MarkerFaceColor','blue','linewidth',1);
y_neg_d2 = abs(posterior_mean_semc2(:,3) - posterior_bounds_sigma_phi_semc2(:,1)); % error in the negative y-direction
y_pos_d2 = abs(posterior_mean_semc2(:,3) - posterior_bounds_sigma_phi_semc2(:,2)); % error in the positive y-direction
errorbar(t_i, posterior_mean_semc2(:,3), y_neg_d2, y_pos_d2, '-sm','MarkerSize',5,...
    'MarkerEdgeColor','magenta','MarkerFaceColor','magenta','linewidth',1);
legend(['Ref. value = ',num2str(sigma_phi, '%.3f'), '^o'], 'linewidth', 2)
xlim([1 4]); ylim([0 10]);
xlabel('$t_{s}$ $[mth]$','Interpreter','latex')
ylabel('$\sigma_{\phi}$ $[deg]$','Interpreter','latex')
set(gca, 'fontsize', 20)

subplot(2,2,3)
hold on; grid on; box on;
plot([1, 4], [sigma_r, sigma_r], 'k--', 'linewidth', 1)
y_neg_c3 = abs(posterior_mean_semc1(:,4) - posterior_bounds_sigma_r_semc1(:,1)); % error in the negative y-direction
y_pos_c3 = abs(posterior_mean_semc1(:,4) - posterior_bounds_sigma_r_semc1(:,2)); % error in the positive y-direction
errorbar(t_i, posterior_mean_semc1(:,4), y_neg_c3, y_pos_c3, '-sb','MarkerSize',5,...
    'MarkerEdgeColor','blue','MarkerFaceColor','blue','linewidth',1);
y_neg_d3 = abs(posterior_mean_semc2(:,4) - posterior_bounds_sigma_r_semc2(:,1)); % error in the negative y-direction
y_pos_d3 = abs(posterior_mean_semc2(:,4) - posterior_bounds_sigma_r_semc2(:,2)); % error in the positive y-direction
errorbar(t_i, posterior_mean_semc2(:,4), y_neg_d3, y_pos_d3, '-sm','MarkerSize',5,...
    'MarkerEdgeColor','magenta','MarkerFaceColor','magenta','linewidth',1);
legend(['Ref. value = ',num2str(sigma_r, '%.3f')], 'linewidth', 2)
xlim([1 4]); ylim([0 1]);
xlabel('$t_{s}$ $[mth]$','Interpreter','latex')
ylabel('$\sigma_{r}$','Interpreter','latex')
set(gca, 'fontsize', 20)

%% SEMC vs SMC Statistics

dim = size(semc_allsamples1,2); % dimensionality of the problem
target_accept = 0.23 + (0.21./dim);

% Plot the acceptance rate values across iterations:
figure;
hold on; box on; grid on;
plot([1 size(semc_allsamples1,3)],[target_accept target_accept] , 'c','linewidth', 1.5)
plot([1 size(semc_allsamples1,3)],[0.15 0.15] , 'k','linewidth', 1.5)
plot((1:size(semc_allsamples1,3)-1)', SEMC1.acceptance, '--bs', 'MarkerFaceColor','b','linewidth', 1.5)
plot((1:size(semc_allsamples2,3)-1)', SEMC2.acceptance, '--ms', 'MarkerFaceColor','m','linewidth', 1.5)
plot([1 size(semc_allsamples1,3)],[0.5 0.5] , 'k','linewidth', 1.5)
legend('Target acceptance rate', 'Optimum acceptance limits', 'SEMC acceptance rates (T_1)',...
       'SEMC acceptance rates (T_2)', 'location', 'northeast', 'linewidth', 2)
xlabel('$j$','Interpreter','latex'); ylabel('Acceptance rate');
xlim([1 4]); ylim([0 1]); set(gca, 'fontsize', 20)

%% Plot Bar Chart Results for the Log-evidence:

log_evidence1 = SEMC1.log_evidence;
log_evidence2 = SEMC2.log_evidence;

figure;
hold on; box on; grid on;
for i = 1:4
subplot(2,2,i)
x1 = log_evidence1(i+1); x2 = log_evidence2(i+1); 
X = categorical({'T_1', 'T_2'});
Y = [x1; x2];
bar(X,Y)
title(sprintf('j = %1.0f \n', i))
ylim([-180 0])
set(gca,'Fontsize',20)
end
