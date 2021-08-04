%% Load Data:

load('DampedOscillator2D_Part1');

%% Analysis of the model updating for a chosen Markov model:

id = 1; % Takes values between 1 to 6: 1 for Markov model 1, 2 for Markov model 2, and so forth.

k_samples_SEMC = zeros(Nsamples,size(logL,1),runs); k_samples_SMC = zeros(Nsamples,size(logL,1),runs);
c_samples_SEMC = zeros(Nsamples,size(logL,1),runs); c_samples_SMC = zeros(Nsamples,size(logL,1),runs);
for r = 1:runs
for i = 1:size(logL,1)
SEMC_struct = SEMC{r,id}; SMC_struct = SMC{r,id}; 
SEMC_allsamples = SEMC_struct.allsamples; SMC_allsamples = SMC_struct.allsamples; % N x 2 x 7
k_samples_SEMC(:,i,r) = SEMC_allsamples(:,1,i+1); k_samples_SMC(:,i,r) = SMC_allsamples(:,1,i+1);
c_samples_SEMC(:,i,r) = SEMC_allsamples(:,2,i+1); c_samples_SMC(:,i,r) = SMC_allsamples(:,2,i+1);    
end    
end

posterior_mean_semc = zeros(size(logL,1),2,runs); posterior_mean_smc = zeros(size(logL,1),2,runs);
posterior_std_semc = zeros(size(logL,1),2,runs); posterior_std_smc = zeros(size(logL,1),2,runs);
posterior_bounds_k_semc = zeros(size(logL,1),2,runs); posterior_bounds_k_smc = zeros(size(logL,1),2,runs);
posterior_bounds_c_semc = zeros(size(logL,1),2,runs); posterior_bounds_c_smc = zeros(size(logL,1),2,runs);
posterior_cov_semc = zeros(size(logL,1),2,runs); posterior_cov_smc = zeros(size(logL,1),2,runs);
for rdx = 1:runs
for idx = 1:size(logL,1)
posterior_mean_semc(idx,1,rdx) = mean(k_samples_SEMC(:,idx,rdx));
posterior_mean_semc(idx,2,rdx) = mean(c_samples_SEMC(:,idx,rdx));
posterior_mean_smc(idx,1,rdx) = mean(k_samples_SMC(:,idx,rdx));
posterior_mean_smc(idx,2,rdx) = mean(c_samples_SMC(:,idx,rdx));

posterior_std_semc(idx,1,rdx) = std(k_samples_SEMC(:,idx,rdx));
posterior_std_semc(idx,2,rdx) = std(c_samples_SEMC(:,idx,rdx));
posterior_std_smc(idx,1,rdx) = std(k_samples_SMC(:,idx,rdx));
posterior_std_smc(idx,2,rdx) = std(c_samples_SMC(:,idx,rdx));

posterior_bounds_k_semc(idx,:,rdx) = prctile(k_samples_SEMC(:,idx,rdx), [5, 95]);
posterior_bounds_c_semc(idx,:,rdx) = prctile(c_samples_SEMC(:,idx,rdx), [5, 95]);
posterior_bounds_k_smc(idx,:,rdx) = prctile(k_samples_SMC(:,idx,rdx), [5, 95]);
posterior_bounds_c_smc(idx,:,rdx) = prctile(c_samples_SMC(:,idx,rdx), [5, 95]);
end
posterior_cov_semc(:,:,rdx) = (posterior_std_semc(:,:,rdx)./posterior_mean_semc(:,:,rdx)).*100;
posterior_cov_smc(:,:,rdx) = (posterior_std_smc(:,:,rdx)./posterior_mean_smc(:,:,rdx)).*100;
end

% Prediction phase:
SMC_samps_k = zeros(Nsamples,runs); SMC_samps_c = zeros(Nsamples,runs);   % SMC samples at the last iteration
SEMC_samps_k = zeros(Nsamples,runs); SEMC_samps_c = zeros(Nsamples,runs); % SEMC samples at the last iteration
for rdx = 1:runs
SMC_samps_k(:,rdx) = k_samples_SMC(:,end,rdx); SMC_samps_c(:,rdx) = c_samples_SMC(:,end,rdx); 
SEMC_samps_k(:,rdx) = k_samples_SEMC(:,end,rdx); SEMC_samps_c(:,rdx) = c_samples_SEMC(:,end,rdx); 
end

SMC_prediction_mean = zeros(3,2,runs); SEMC_prediction_mean = zeros(3,2,runs);
SMC_prediction_cov = zeros(3,2,runs); SEMC_prediction_cov = zeros(3,2,runs);
prediction_bounds_k_smc = zeros(3,2,runs); prediction_bounds_c_smc = zeros(3,2,runs);
prediction_bounds_k_semc = zeros(3,2,runs); prediction_bounds_c_semc = zeros(3,2,runs);
for rdx = 1:runs
theta_smc = [SMC_samps_k(:,rdx), SMC_samps_c(:,rdx)]; theta_semc = [SEMC_samps_k(:,rdx), SEMC_samps_c(:,rdx)]; 
for i=1:3
if id == 1
dynamic_model = markov_model{1,r};
elseif id == 2
dynamic_model = markov_model{2,r};
elseif id == 3
dynamic_model = markov_model{3,r};
elseif id == 4
dynamic_model = markov_model{4,r};
elseif id == 5
dynamic_model = markov_model{5,r};
elseif id == 6
dynamic_model = markov_model{6,r};
end    
SMC_prediction = dynamic_model([theta_smc(:,1), theta_smc(:,2)]); 
SEMC_prediction = dynamic_model([theta_semc(:,1), theta_semc(:,2)]);

for kd = 1:size(SMC_prediction,1)
    
if SMC_prediction(kd,2) < 0
SMC_prediction(kd,2) = 0;
end

if SEMC_prediction(kd,2) < 0
SEMC_prediction(kd,2) = 0;
end

end

theta_smc = SMC_prediction; theta_semc = SEMC_prediction; 
SMC_prediction_mean(i,:,rdx) = mean(theta_smc); SEMC_prediction_mean(i,:,rdx) = mean(theta_semc);
SMC_prediction_cov(i,:,rdx) = (std(theta_smc)./mean(theta_smc)).*100; SEMC_prediction_cov(i,:,rdx) = (std(theta_semc)./mean(theta_semc)).*100; 
prediction_bounds_k_smc(i,:,rdx) = prctile(theta_smc(:,1), [5, 95]); prediction_bounds_k_semc(i,:,rdx) = prctile(theta_semc(:,1), [5, 95]);
prediction_bounds_c_smc(i,:,rdx) = prctile(theta_smc(:,2), [5, 95]); prediction_bounds_c_semc(i,:,rdx) = prctile(theta_semc(:,2), [5, 95]);
end
end

%% To plot the estimates of k(t_i) and c(t_i) with the Markov model id (SEMC vs SMC):

figure;
for i = 1:runs
subplot(3,4,i)
hold on; grid on; box on;
fill([7 7 9 9],[0 7 7 0],'y','facealpha',0.3, 'handlevisibility', 'off')
stairs(time(:,i,1), k(:,i), 'color', '#C0C0C0', 'linewidth',1, 'handlevisibility', 'off')
stairs([time(end,i,1),9], [k(end,i),k(end,i)], 'color', '#C0C0C0', 'linewidth',1)
plot((1:1:6)', true_val(:,1,i), 'x', 'color', '#7E2F8E', 'linewidth', 2)
plot((7:1:9)', pre_val(:,1,i), 'x', 'color', '#7E2F8E', 'linewidth', 2,'handlevisibility', 'off')
y_neg_a1 = abs(posterior_mean_smc(:,1,i) - posterior_bounds_k_smc(:,1,i)); % error in the negative y-direction
y_pos_a1 = abs(posterior_mean_smc(:,1,i) - posterior_bounds_k_smc(:,2,i)); % error in the positive y-direction
y_neg_a2 = abs(SMC_prediction_mean(:,1,i) - prediction_bounds_k_smc(:,1,i)); % error in the negative y-direction
y_pos_a2 = abs(SMC_prediction_mean(:,1,i) - prediction_bounds_k_smc(:,2,i)); % error in the positive y-direction
errorbar((1:1:9)', [posterior_mean_smc(:,1,i); SMC_prediction_mean(:,1,i)], [y_neg_a1; y_neg_a2], [y_pos_a1; y_pos_a2], '-sb','MarkerSize',5,...
    'MarkerEdgeColor','blue','MarkerFaceColor','blue', 'linewidth',1);
y_neg_a3 = abs(posterior_mean_semc(:,1,i) - posterior_bounds_k_semc(:,1,i)); % error in the negative y-direction
y_pos_a3 = abs(posterior_mean_semc(:,1,i) - posterior_bounds_k_semc(:,2,i)); % error in the positive y-direction
y_neg_a4 = abs(SEMC_prediction_mean(:,1,i) - prediction_bounds_k_semc(:,1,i)); % error in the negative y-direction
y_pos_a4 = abs(SEMC_prediction_mean(:,1,i) - prediction_bounds_k_semc(:,2,i)); % error in the positive y-direction
errorbar((1:1:9)', [posterior_mean_semc(:,1,i);SEMC_prediction_mean(:,1,i)], [y_neg_a3; y_neg_a4], [y_pos_a3; y_pos_a4], '-sr','MarkerSize',5,...
    'MarkerEdgeColor','red','MarkerFaceColor','red','linewidth',1);
xlabel(['$t$ $[mth]$ $|$ Run $\#$', num2str(i, '%d')], 'Interpreter', 'latex'); ylabel('$k(t)$ $[N/m]$', 'Interpreter', 'latex')
xlim([0 9]); xticks((0:1:9)); set(gca, 'fontsize', 15)
end
legend('Random process function','True values','SMC estimates','SEMC estimates','linewidth', 2)

figure;
for i = 1:runs
subplot(3,4,i)
hold on; grid on; box on;
fill([7 7 9 9],[0 2 2 0],'y','facealpha',0.3, 'handlevisibility', 'off')
stairs(time(:,i,2), c(:,i), 'color', '#C0C0C0', 'linewidth',1,'handlevisibility', 'off')
stairs([time(end,i,2),9], [c(end,i),c(end,i)], 'color', '#C0C0C0', 'linewidth',1)
plot((1:1:6)', true_val(:,2,i), 'x', 'color', '#7E2F8E', 'linewidth', 2)
plot((7:1:9)', pre_val(:,2,i), 'x', 'color', '#7E2F8E', 'linewidth', 2,'handlevisibility', 'off')
y_neg_b1 = abs(posterior_mean_smc(:,2,i) - posterior_bounds_c_smc(:,1,i)); % error in the negative y-direction
y_pos_b1 = abs(posterior_mean_smc(:,2,i) - posterior_bounds_c_smc(:,2,i)); % error in the positive y-direction
y_neg_b2 = abs(SMC_prediction_mean(:,2) - prediction_bounds_c_smc(:,1)); % error in the negative y-direction
y_pos_b2 = abs(SMC_prediction_mean(:,2) - prediction_bounds_c_smc(:,2)); % error in the positive y-direction
errorbar((1:1:9)', [posterior_mean_smc(:,2,i);SMC_prediction_mean(:,2)], [y_neg_b1; y_neg_b2], [y_pos_b1; y_pos_b2], '-sb','MarkerSize',5,...
    'MarkerEdgeColor','blue','MarkerFaceColor','blue', 'linewidth',1);
y_neg_b3 = abs(posterior_mean_semc(:,2,i) - posterior_bounds_c_semc(:,1,i)); % error in the negative y-direction
y_pos_b3 = abs(posterior_mean_semc(:,2,i) - posterior_bounds_c_semc(:,2,i)); % error in the positive y-direction
y_neg_b4 = abs(SEMC_prediction_mean(:,2) - prediction_bounds_c_semc(:,1)); % error in the negative y-direction
y_pos_b4 = abs(SEMC_prediction_mean(:,2) - prediction_bounds_c_semc(:,2)); % error in the positive y-direction
errorbar((1:1:9)', [posterior_mean_semc(:,2,i); SEMC_prediction_mean(:,2)], [y_neg_b3; y_neg_b4], [y_pos_b3; y_pos_b4], '-sr','MarkerSize',5,...
    'MarkerEdgeColor','red','MarkerFaceColor','red', 'linewidth',1);
xlabel(['$t$ $[mth]$ $|$ Run $\#$', num2str(i, '%d')], 'Interpreter', 'latex')
ylabel('$c(t)$ $[Ns/m]$', 'Interpreter', 'latex')
xlim([0 9]); xticks((0:1:9));
set(gca, 'fontsize', 15)
end
legend('Random process function','True values','SMC estimates','SEMC estimates','linewidth', 2)

%% SEMC vs SMC Statistics (Plot the acceptance rate values across iterations):

dim = 2; % dimensionality of the problem
target_accept = 0.23 + (0.21./dim);

figure;
for r = 1:10
SEMC_struct = SEMC{r,id}; SMC_struct = SMC{r,id}; 
subplot(3,4,r)
hold on; box on; grid on;
plot([1 size(logL,1)],[target_accept target_accept] , 'c','linewidth', 1.5)
plot([1 size(logL,1)],[0.15 0.15] , 'k','linewidth', 1.5)
plot((1:size(logL,1))', SEMC_struct.acceptance, '--rs', 'MarkerFaceColor','r','linewidth', 1.5)
plot((1:size(logL,1))', SMC_struct.acceptance, '--bs', 'MarkerFaceColor','b','linewidth', 1.5)
plot([1 size(logL,1)],[0.5 0.5] , 'k','linewidth', 1.5)
xlabel(['$j$ $|$ Run $\#$', num2str(r, '%d')], 'Interpreter', 'latex');  ylabel('Acceptance rate');
xlim([1 6]); xticks([1,2,3,4,5,6]); ylim([0 1]); set(gca, 'fontsize', 15)
end
legend('Target acceptance rate', 'Optimum acceptance limits', ['SEMC acceptance rate (T_', num2str(id, '%d'), ')'],...
       ['SMC acceptance rate (T_', num2str(id, '%d'), ')'], 'linewidth', 2)

%% Plot the Log-evidence across runs 
 
% Plot for SEMC:
figure;
hold on; box on; grid on;
for i = 1:10
subplot(3,4,i)
hold on; box on; grid on;
for j = 1:6
log_evidence_semc = SEMC{i,j}.log_evidence;
plot((1:6)', log_evidence_semc(2:end), '--s', 'linewidth', 1.5)
end
xlabel(['$j$ $|$ Run $\#$', num2str(i, '%d')], 'Interpreter', 'latex');  ylabel('Log-evidence');
xlim([1 6]); xticks([1,2,3,4,5,6]); set(gca, 'fontsize', 15)
end
legend('Model 1', 'Model 2', 'Model 3', 'Model 4', 'Model 5', 'Model 6', 'linewidth', 2)

% Plot for SMC:
figure;
hold on; box on; grid on;
for i = 1:10
subplot(3,4,i)
hold on; box on; grid on;
for j = 1:6
log_evidence_smc = SMC{i,j}.log_evidence;
plot((1:6)', log_evidence_smc(2:end), '--s', 'linewidth', 1.5)
end
xlabel(['$j$ $|$ Run $\#$', num2str(i, '%d')], 'Interpreter', 'latex');  ylabel('Log-evidence');
xlim([1 6]); xticks([1,2,3,4,5,6]); set(gca, 'fontsize', 15)
end
legend('Model 1', 'Model 2', 'Model 3', 'Model 4', 'Model 5', 'Model 6', 'linewidth', 2)
