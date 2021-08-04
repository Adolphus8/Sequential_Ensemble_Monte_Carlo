%% The SEMC sampler
%
% The SEMC sampler is based on the original Sequential Monte
% Carlo (SMC) sampling class (see paper by Moral et. al (2006): Sequential 
% Monte Carlo Samplers - https://www.jstor.org/stable/3879283) and employs
% the use of the Affine-invariant Ensemble Sampler (AIES) proposed by 
% Goodman and Weare (2010) to update the samples at each iteration.
%
%% Bayesian Inference of a Time-varying parameter:
%
% In this example, we will evaluate the performance of the SEMC sampler
% in estimating and predicting a time-varying parameter at different time-step. 
% The system of interest will be a simple damped spring-mass system. The 
% dynamical data of the displacement of the mass from its rest position will 
% be obtained between t=0s and t=5s at intervals of dt=0.01s. The system is
% taken to be underdamped.
%
% The simple displacement model is:
%
% x(t) = A*exp(-gamma*t)*cos(omega*t + d); we set d = 0 rad for simplicity;
%
% omega = sqrt((k./m) - (c./(2.*m)).^2); 
% gamma = c./(2.*m); 
% where m = 0.3 kg (mass of the block attached to spring);
%
% The spring stiffness and damper are time-varying and weakens over time
% following a random process.
%
% The spring stiffnes and damper are checked every month and it is assumed 
% that during the process of data-collection at every "inspection" time-step, 
% t_i, the value of k and c are constant.
%
%% Define the parameters and random variables:

m = 0.3;             % Mass of the blocks in [kg]
A = 0.05;            % Displacement of the spring-mass system/Initial amplitude [m]
t_d = (0:0.01:5)';   % Time-input of the dynamical data [s]
t_i = [1,2,3,4,5,6]; % Inspection time [mths]

%% To plot the evolution model of k(t) and c(t):

runs = 10;               % No. of simulated runs;
mu = 4;                  % Mean of the exponential distribution of the time input parameter;
time = zeros(10,runs,2); % Time input parameter of the parameter evolution [mths];
k = zeros(10,runs);      % True values of stiffness [N/m];
c = zeros(10,runs);      % True values of damping [Ns/m];
random_fac_k = @(N) unifrnd(0.85,0.95,N,1); random_fac_c = @(N) unifrnd(0.50,0.95,N,1);

% Generate initial random values of k [N/m] and c [Ns/m]:
k(1,:) = 3.0; c(1,:) = 0.6; 
for i = 1:runs
time(2:10,i,1) = exp_rnd(mu,[0,9],9,1); time(2:10,i,2) = exp_rnd(mu,[0,9],9,1);
for j = 2:10
k(j,i) = k(j-1,i) .* random_fac_k(1); c(j,i) = c(j-1,i) .* random_fac_c(1);
end
end

% Generate noisy measurements of dynamical Displacement:

% Define model for the eigenfrequency, omega:
omega = @(k,c,m) sqrt((k./m) - (c./(2.*m)).^2);

% Define model for simple harmonic oscillator:
delta = 0; % Initial phase term when time starts [rad] 
displacement = @(t_d,k,c) A.*exp(- (c./(2.*m)).*t_d).*cos((omega(k,c,m) .* t_d) + delta); % in [m]
sigma = 0.005; % The noise in the measurement is set at 10% of the amplitude [m]

% To generate noisy data-set of the dynamical system for different time steps, t_i:
data = zeros(length(t_d),length(t_i),runs); 
true_val = zeros(length(t_i),3,runs);
for rdx = 1:runs
for idx = 1:length(t_i)
ins = t_i(idx); % Inspection time
k_idx = find(time(:,rdx,1) <= ins); c_idx = find(time(:,rdx,2) <= ins);
for jdx = 1:length(t_d)
true_val(idx,1,rdx) = k(k_idx(end),rdx); true_val(idx,2,rdx) = c(c_idx(end),rdx); true_val(idx,3,rdx) = 2*sqrt(k(k_idx(end),rdx).*m); 
data(jdx,idx,rdx) = displacement(t_d(jdx),k(k_idx(end),rdx),c(c_idx(end),rdx)) + sigma .* randn(1,1);
end
end
end

% To compute the root-mean-square-error of the data relative to the
% theoretial model of the SHM at each time-step t_i:
rmse = zeros(6,runs);
for rdx = 1:runs
for idx = 1:length(t_i)
ins = t_i(idx); % Inspection time
k_idx = find(time(:,rdx,1) <= ins); c_idx = find(time(:,rdx,2) <= ins);
square_error = (data(:,idx,rdx) - displacement(t_d(1:length(t_d)),k(k_idx(end),rdx),c(c_idx(end),rdx))).^2;
rmse(idx,rdx) = sqrt(mean(square_error));
end
end

%% Plot figures:

t_p = [7,8,9]; % Prediction time
pre_val = zeros(length(t_p),2,runs); % True prediction values of k and c
for rdx = 1:runs
for idx = 1:length(t_p)
pre = t_p(idx); % Prediction time
k_idx = find(time(:,rdx,1) <= pre); c_idx = find(time(:,rdx,2) <= pre);
for jdx = 1:length(t_d)
pre_val(idx,1,rdx) = k(k_idx(end),rdx); pre_val(idx,2,rdx) = c(c_idx(end),rdx); 
end
end
end

figure;
for i = 1:runs
subplot(3,4,i)
hold on; grid on; box on;
stairs(time(:,i,1), k(:,i), 'color', '#C0C0C0', 'linewidth',1)
stairs([time(end,i,1),9], [k(end,i),k(end,i)], 'color', '#C0C0C0', 'linewidth',1, 'handlevisibility', 'off')
plot((1:1:6)', true_val(:,1,i), 'x', 'color', '#7E2F8E', 'linewidth', 2)
plot((7:1:9)', pre_val(:,1,i), 'x', 'color', '#7E2F8E', 'linewidth', 2)
xlabel(['$t$ $[mth]$ $|$ Run $\#$', num2str(i, '%d')], 'Interpreter', 'latex')
ylabel('$k(t)$ $[N/m]$', 'Interpreter', 'latex')
xlim([0 9]); ylim([1,3]); xticks([0,1,2,3,4,5,6,7,8,9]);
set(gca, 'fontsize', 15)
end
legend('Random process function','True values','linewidth', 2)

figure;
for i = 1:runs
subplot(3,4,i)
hold on; grid on; box on;
stairs(time(:,i,2), c(:,i), 'color', '#C0C0C0', 'linewidth',1)
stairs([time(end,i,2),9], [c(end,i),c(end,i)], 'color', '#C0C0C0', 'linewidth',1,'handlevisibility', 'off')
plot((1:1:6)', true_val(:,2,i), 'x', 'color', '#7E2F8E', 'linewidth', 2)
plot((7:1:9)', pre_val(:,2,i), 'x', 'color', '#7E2F8E', 'linewidth', 2)
xlabel(['$t$ $[mth]$ $|$ Run $\#$', num2str(i, '%d')], 'Interpreter', 'latex')
ylabel('$c(t)$ $[Ns/m]$', 'Interpreter', 'latex')
xlim([0 9]); ylim([0,0.55]); xticks([0,1,2,3,4,5,6,7,8,9]);
set(gca, 'fontsize', 15)
end
legend('Random process function','True values','linewidth', 2)

%% Bayesian Model Updating set-up:

% Define the prior distribution:
lowerBound = [0.001, 0.001]; 
upperBound = [100, 10];

prior_pdf_k = @(x) unifpdf(x,lowerBound(1),upperBound(1)); % Prior PDf for k
prior_pdf_c = @(x) unifpdf(x,lowerBound(2),upperBound(2)); % Prior PDF for c
prior_pdf = @(x) prior_pdf_k(x(:,1)) .* prior_pdf_c(x(:,2)); % Overall Prior PDF

prior_rnd = @(N) [unifrnd(lowerBound(1),upperBound(1),N,1),... 
                  unifrnd(lowerBound(2),upperBound(2),N,1)];

% Define the likelihood function:
scale = 7;
measurement_model = @(x,t_d) A.*exp(- (x(:,2)./(2.*m)).*t_d).*...
                            cos((sqrt((x(:,1)./m) - (x(:,2)./(2.*m)).^2).* t_d) + delta); % in [m]
             
loglike = @(x,t_i,r) - 0.5 .* (1./(scale.*rmse(t_i,r))).^2 .*(sum((data(:,t_i,r) - measurement_model(x, t_d)).^2));

logL = cell(length(t_i),runs);
for jdx = 1:runs
for idx = 1:length(t_i)
logL{idx,jdx} = @(x) loglike(x,idx,jdx);
end
end

%% Perform Online Sequential Bayesian Updating via SEMC and SMC:

SEMC = cell(runs,6); SMC = cell(runs,6);
timeSEMC = zeros(runs,6); timeSMC = zeros(runs,6);

% Initialise:
Nsamples = 1000;

markov_model = cell(6,runs);
for rdx = 1:runs
% Define Markov Model 1 for both k and c (Exponential model):
data_time_k = time(1:7,rdx,1); data_k = k(1:7,rdx);
data_time_c = time(1:7,rdx,2); data_c = c(1:7,rdx);
fitob1_k = fit(data_time_k, data_k,'exp1'); fitob1_c = fit(data_time_c, data_c,'exp1');
sigma_k1 = sqrt(mean((((fitob1_k.a).*exp((fitob1_k.b).*(data_time_k))) - data_k).^2)); sigma_k1 = round(sigma_k1,1);
sigma_c1 = sqrt(mean((((fitob1_c.a).*exp((fitob1_c.b).*(data_time_c))) - data_c).^2)); sigma_c1 = round(sigma_c1,2);
a0_k(rdx) = fitob1_k.a; a0_c(rdx) = fitob1_c.a;
a1(1,rdx) = round(fitob1_k.b,1); a2(1,rdx) = round(fitob1_c.b,1); a3(rdx) = sigma_k1; a4(rdx) = sigma_c1;
k_new1 = @(k_old) k_old.*exp(-abs(a1)) + sigma_k1.*randn(Nsamples,1); c_new1 = @(c_old) c_old.*exp(-abs(a2)) + sigma_c1.*randn(Nsamples,1); 
markov_model{1,rdx} = @(x) [k_new1(x(:,1)), c_new1(x(:,2))];
a1(2,rdx) = round(((fitob1_k.b).*0.8),1); a2(2,rdx) = round(((fitob1_c.b).*0.8),1); 
k_new1 = @(k_old) k_old.*exp(-abs(a1)) + sigma_k1.*randn(Nsamples,1); c_new1 = @(c_old) c_old.*exp(-abs(a2)) + sigma_c1.*randn(Nsamples,1); 
markov_model{2,rdx} = @(x) [k_new1(x(:,1)), c_new1(x(:,2))];
a1(3,rdx) = round((fitob1_k.b).*1.5,1); a2(3,rdx) = round((fitob1_c.b).*1.5,1); 
k_new1 = @(k_old) k_old.*exp(-abs(a1)) + sigma_k1.*randn(Nsamples,1); c_new1 = @(c_old) c_old.*exp(-abs(a2)) + sigma_c1.*randn(Nsamples,1); 
markov_model{3,rdx} = @(x) [k_new1(x(:,1)), c_new1(x(:,2))];

% Define Markov Model 2 for both k and c (Linear model):
fitob2_k = polyfit(data_time_k, data_k,1); fitob2_c = polyfit(data_time_c, data_c,1);
sigma_k2 = sqrt(mean(((fitob2_k(1).*(data_time_k) + fitob2_k(2)) - data_k).^2)); sigma_k2 = round(sigma_k2,1);
sigma_c2 = sqrt(mean(((fitob2_c(1).*(data_time_c) + fitob2_c(2)) - data_c).^2)); sigma_c2 = round(sigma_c2,2);
b0_k(rdx) = fitob2_k(2); b0_c(rdx) = fitob2_c(2);
b1(1,rdx) = round(fitob2_k(1),1); b2(1,rdx) = round(fitob2_c(1),2); b3(rdx) = sigma_k2; b4(rdx) = sigma_c2;
k_new2 = @(k_old) k_old - abs(b1) + sigma_k2.*randn(Nsamples,1); c_new2 = @(c_old) c_old - abs(b2) + sigma_c2.*randn(Nsamples,1); 
markov_model{4,rdx} = @(x) [k_new2(x(:,1)), c_new2(x(:,2))];
b1(2,rdx) = round((fitob2_k(1).*0.8),1); b2(2,rdx) = round((fitob2_c(1).*0.8),2); 
k_new2 = @(k_old) k_old - abs(b1) + sigma_k2.*randn(Nsamples,1); c_new2 = @(c_old) c_old - abs(b2) + sigma_c2.*randn(Nsamples,1);  
markov_model{5,rdx} = @(x) [k_new2(x(:,1)), c_new2(x(:,2))];
b1(3,rdx) = round((fitob2_k(1).*1.5),1).*1.5; b2(3,rdx) = round((fitob2_c(1).*1.5),2); 
k_new2 = @(k_old) k_old - abs(b1) + sigma_k2.*randn(Nsamples,1); c_new2 = @(c_old) c_old - abs(b2) + sigma_c2.*randn(Nsamples,1);  
markov_model{6,rdx} = @(x) [k_new2(x(:,1)), c_new2(x(:,2))];
end

%%

%parpool(15);
for r = 1:runs
fprintf('Run no.: %d \n',r)

logl = cell(6,1);
for i = 1:6
logl{i} = logL{i,r};
end

% Start SEMC sampler:
for jdx = 1:6
tic;
SEMC{r,jdx} = SEMCsampler('nsamples',Nsamples,'loglikelihoods',logl,...
                   'dynamic_model',markov_model{jdx,rdx},'priorpdf',prior_pdf,...
                   'priorrnd',prior_rnd,'burnin',0,'stepsize',10);
timeSEMC(r,jdx) = toc;
fprintf('Time elapsed is for the SEMC sampler: %f \n',timeSEMC(r,jdx))
end
end

for r = 1:runs
fprintf('Run no.: %d \n',r)

logl = cell(6,1);
for i = 1:6
logl{i} = logL{i,r};
end

% Start SMC sampler:
for jdx = 1:6
tic;
SMC{r,jdx} = SMCsampler('nsamples',Nsamples,'loglikelihoods',logl,...
                 'dynamic_model',markov_model{jdx,rdx},'priorpdf',prior_pdf,...
                 'priorrnd',prior_rnd,'burnin',0);
timeSMC(r,jdx) = toc;
fprintf('Time elapsed is for the SMC sampler: %f \n',timeSMC(r,jdx))
end
end

%% Save the data:

save('DampedOscillator2D_Part1');
