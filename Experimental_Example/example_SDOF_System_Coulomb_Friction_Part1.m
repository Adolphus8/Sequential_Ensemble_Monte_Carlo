%% Experimental: SDOF System with Coulomb Friction
%
% This set-up is based on the Single Degree-of-Freedom Dynamical System
% subjected to Coulomb Friction Force that is presented in the literature:
%
% L. Marino and A. Cicirello (2020). Experimental investigation of a single-
% degree-of-freedom system with Coulomb friction. Nonlinear Dynamics, 99(3), 
% 1781-1799. doi: 10.1007/s11071-019-05443-2
%
% This set-up is applicable only for the Base motion (with fixed wall) case.
%
% Here, Bayesian Inference will be done using actual experimental data.
%-------------------------------------------------------------------------%
%% Load Numerical Data:
load('experimental_data');

%% Define key parameters:
omega_n = data.omega_n;                % Natural frequency of the structure [rad/s]
driving_force = 2.5;                   % Driving force amplitude by the rotor [N]
r_nom = data.r_nom;                    % Nominal values of the dimensionless input frequency ratios

%% Define the data-set:

coulomb_force_nom = [1.375; 1.000; 0.625; 0.250]; % Nominal values of Coulomb Friction [N]
coulomb_force = data.coulomb_friction;            % Experimental values of Coulomb Frictions [N]
r_exp = data.r_exp;                               % Experimental values of Frequency ratios
phi_exp = data.phi;                               % Experimental values of Phase angles [deg]
beta_exp = data.beta;                             % Experimental values of Force ratio

% Generate Analytical solution through Den-Hartog's solution:
output = DenHartogHarmonic(beta_exp');
r_an = output.frequency_ratios; % The output frequency ratios
phi_an = output.phase_angles;   % The output analytical phase angles
phi_bound = output.phase_bound; % The phase angle bound defined by Den-Hartog's Boundary

%% To plot the Phase angles vs Frequency ratio curves:
colors = [0 0 1; 0 0.5 0; 1 0 0; 1 0 1];

figure; 
hold on; box on; grid on;
for ib = 1:length(beta_exp) % To plot for different Friction Force ratio
% Analytical plot:
plot(r_an, phi_an(ib+1,:), '-', 'color', colors(ib,:), 'linewidth', 1, 'handlevisibility', 'off');
% Numerical scatterplot: 
plot(r_exp(:,ib), phi_exp(:,ib), 'x', 'color', colors(ib,:), 'linewidth', 2); 
end
plot(r_an, phi_bound,'-- k');
legend(['Data for F_{\mu}(t_{1}) = ',num2str(coulomb_force(1), '%.3f'), ' N'], ['Data for F_{\mu}(t_{2}) = ',num2str(coulomb_force(2), '%.3f'), ' N'],...
['Data for F_{\mu}(t_{3}) = ',num2str(coulomb_force(3), '%.3f'), ' N'],['Data for F_{\mu}(t_{4}) = ',num2str(coulomb_force(4), '%.3f'), ' N'],...
'Den-Hartog''s Boundary',...
'linewidth', 2, 'location', 'southeast');
xlabel('$r$', 'Interpreter', 'latex'); 
ylabel('$\phi$ $[deg]$', 'Interpreter', 'latex');
xlim([0, 2]); ylim([0 180]);
set(gca,'FontSize',20);

%% Bayesian Model Updating Set-up:
% The epistemic parameters to be inferred are the following: 
% {Coulomb Force, Natural Frequency, Frequency Ratio Noise, Phase Angle Noise}

% Define the Prior distribution:
lowerbound = [0.01, 0.001, 0.001]; upperbound = [100, 10, 1];
prior_coulomb = @(x) unifpdf(x, lowerbound(1), upperbound(1));   % Prior for Coulomb Friction
prior_omega = @(x) unifpdf(x, lowerbound(1), upperbound(1));     % Prior for Natural Frequency
prior_sigma_phi = @(x) unifpdf(x, lowerbound(2), upperbound(2)); % Prior for Phase Angle Noise
prior_sigma_r = @(x) unifpdf(x, lowerbound(3), upperbound(3));   % Prior for Frequency Ratio Noise

prior_pdf = @(x) prior_coulomb(x(:,1)) .* prior_omega(x(:,2)) .* ...
                 prior_sigma_phi(x(:,3)) .* prior_sigma_r(x(:,4));
prior_rnd = @(N) [unifrnd(lowerbound(1), upperbound(1), N, 1), ...
                  unifrnd(lowerbound(1), upperbound(1), N, 1), ...
                  unifrnd(lowerbound(2), upperbound(2), N, 1), ...
                  unifrnd(lowerbound(3), upperbound(3), N, 1)];          

% Define the loglikelihood function:
model = @(x, r_exp) blackbox_model(x, r_exp, driving_force);

data_frequency = r_exp .* omega_n;       % Experimental driving requency in [rad/s]
loglike = @(x,t) loglikelihood(x, model, phi_exp(:,t), data_frequency(:,t), r_exp(:,t), r_nom);

time_step = size(coulomb_force,1);
logL = cell(time_step,1);
logL{1,1} = @(x) loglike(x,1); logL{2,1} = @(x) loglike(x,2);
logL{3,1} = @(x) loglike(x,3); logL{4,1} = @(x) loglike(x,4);
           
%% Define the exponential nominal evolution function for Dynamic Model 2:

fitob2 = fit((1:4)', coulomb_force, 'exp1');
force_evolution2 = @(t) fitob2.a * exp(fitob2.b * t);
                    
%% Plot the Evolution of F_mu:

t_i = [1;2;3;4]; % Time-step values

% To plot the estimates of F(t_i) with the theoretical model (Model 1 vs Model 2):
figure;
% Subplot for Model 1:
subplot(1,2,1)
hold on; grid on; box on;
plot(t_i, coulomb_force, 'x',  'color', '#7E2F8E', 'linewidth', 2)
plot(t_i, coulomb_force_nom, 'k--', 'linewidth', 1)
xlim([1 4]); ylim([0.2 1.6]);
xlabel('$t_{s}$ $[mth]$','Interpreter','latex')
ylabel('$F_{\mu}(t_{s})$ $[N]$','Interpreter','latex')
legend('True values', 'Evolution model \Gamma_{1}', 'linewidth', 2)
set(gca, 'fontsize', 20)

% Subplot for Model 2:
subplot(1,2,2)
hold on; grid on; box on;
plot(t_i, coulomb_force, 'x',  'color', '#7E2F8E', 'linewidth', 2)
plot((1:0.01:4), force_evolution2((1:0.01:4)), 'k--', 'linewidth', 1)
xlim([1 4]); ylim([0.2 1.6]);
xlabel('$t_{s}$ $[mth]$','Interpreter','latex')
ylabel('$F_{\mu}(t_{s})$ $[N]$','Interpreter','latex')
legend('True values', 'Evolution model \Gamma_{2}', 'linewidth', 2)
set(gca, 'fontsize', 20)
                      
%% Perform Online Sequential Bayesian Updating via SEMC (with Dynamical Model 1):

% Initialise:
Nsamples = 1000;

% Define the Dynamic Model 1:
sigma_1 = sqrt(mean((coulomb_force - coulomb_force_nom).^2));    % Stdev of process noise estimated from rmse
F_new_1 = @(F_old) (F_old - 0.375) + sigma_1.*randn(Nsamples,1); % Evolution model for Coulomb Force
omega_new = @(omega_old) omega_old;             % Evolution model for Natural Frequency
sigma_phi_new = @(sigma_phi_old) sigma_phi_old; % Evolution model for Phase Angle measurement noise
sigma_r_new = @(sigma_r_old) sigma_r_old;       % Evolution model for Frequency Ratio measurement noise

dynamic_model_1 = @(x) [F_new_1(x(:,1)), omega_new(x(:,2)), ...
                          sigma_phi_new(x(:,3)), sigma_r_new(x(:,4))];

% Start SEMC sampler:
tic;
SEMC1 = SEMCsampler('nsamples',Nsamples,'loglikelihoods',logL,...
                   'dynamic_model',dynamic_model_1,'priorpdf',prior_pdf,...
                   'priorrnd',prior_rnd,'burnin',0,'stepsize',30);
semc_allsamples1 = SEMC1.allsamples;
SEMC1.acceptance;
timeSEMC1 = toc;
fprintf('Time elapsed is for the SEMC sampler: %f \n',timeSEMC1)

%% Perform Online Sequential Bayesian Updating via SEMC (with Dynamical Model 2):

% Initialise:
Nsamples = 1000;

% Define the Dynamic Model 2:
sigma_2 = sqrt(mean((coulomb_force - force_evolution2([1,2,3,4]')).^2));  % Stdev of process noise estimated from rmse
F_new_2 = @(F_old) (exp(fitob2.b) .* F_old) + sigma_2.*randn(Nsamples,1); % Evolution model for Coulomb Force
omega_new = @(omega_old) omega_old;             % Evolution model for Natural Frequency
sigma_phi_new = @(sigma_phi_old) sigma_phi_old; % Evolution model for Phase Angle measurement noise
sigma_r_new = @(sigma_r_old) sigma_r_old;       % Evolution model for Frequency Ratio measurement noise
dynamic_model_2 = @(x) [F_new_2(x(:,1)), omega_new(x(:,2)), ...
                          sigma_phi_new(x(:,3)), sigma_r_new(x(:,4))];

% Start SEMC sampler:
tic;
SEMC2 = SEMCsampler('nsamples',Nsamples,'loglikelihoods',logL,...
                   'dynamic_model',dynamic_model_2,'priorpdf',prior_pdf,...
                   'priorrnd',prior_rnd,'burnin',0,'stepsize',30);
semc_allsamples2 = SEMC2.allsamples;
SEMC2.acceptance;
timeSEMC2 = toc;
fprintf('Time elapsed is for the SEMC sampler: %f \n',timeSEMC2)

%% Save Data:

save('experimental_SDOF_System_Coulomb_Friction_Part1_m');
