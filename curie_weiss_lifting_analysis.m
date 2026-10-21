% =========================================================================
% NON-REVERSIBLE MCMC: LIFTED METROPOLIS-HASTINGS ON THE CURIE-WEISS MODEL
% Author: Riccardo Puschetti
% Description: Spectral analysis, error decay, and stochastic trajectories 
%              (Zig-Zag PDMP limit) comparing standard reversible MH and 
%              Lifted MH at the critical temperature (beta = 1.0).
% =========================================================================
clear; clc; close all;

% Critical inverse temperature (phase transition regime)
beta = 1.0;  

%% PART 1: ASYMPTOTIC SCALING ANALYSIS (N^{3/2} vs N^{3/4})
fprintf('1. Computing singular values and relaxation times...\n');

N_vec = [20, 40, 80, 160, 320, 640];  
tau_rev  = zeros(length(N_vec), 1); 
tau_lift = zeros(length(N_vec), 1);

for idx = 1 : length(N_vec)
    N = N_vec(idx);
    k = (0 : N)';  % Number of positive spins (+1)
    
    % Macroscopic energy and log-target distribution (using gammaln to prevent overflow)
    E = - (beta / (2*N)) * (2*k - N).^2;  
    log_pi = (gammaln(N+1) - gammaln(k+1) - gammaln(N-k+1)) - beta * E;  
    
    % --- 1A. Build Reversible Transition Matrix (P_rev) ---
    P_rev = zeros(N+1, N+1);
    for i = 1 : (N+1)
        k_val = i - 1;
        if k_val < N 
            P_rev(i, i+1) = ((N - k_val) / N) * min(1, exp(-beta * (E(i+1) - E(i))));
        end
        if k_val > 0
            P_rev(i, i-1) = (k_val / N) * min(1, exp(-beta * (E(i-1) - E(i))));
        end
        P_rev(i, i) = 1 - sum(P_rev(i, :));
    end
    
    % --- 1B. Build Lifted Non-Reversible Transition Matrix (P_lift) ---
    P_lift = zeros(2*N+2, 2*N+2);  
    for i = 1 : (N+1)
        if i < N+1, T_plus = P_rev(i, i+1); else, T_plus = 0; end
        if i > 1,   T_minus = P_rev(i, i-1); else, T_minus = 0; end
        
        % Turitsyn-Chertkov-Vucelja switching rates
        T_plus_minus = max(0, T_minus - T_plus);
        T_minus_plus = max(0, T_plus - T_minus);
        
        % Forward sector (+1)
        if i < N+1, P_lift(i, i+1) = T_plus; end
        P_lift(i, i + N+1) = T_plus_minus;  
        P_lift(i, i)       = 1 - T_plus - T_plus_minus;
        
        % Backward sector (-1)
        if i > 1, P_lift(i + N+1, i - 1 + N+1) = T_minus; end
        P_lift(i + N+1, i)       = T_minus_plus;
        P_lift(i + N+1, i + N+1) = 1 - T_minus - T_minus_plus;
    end
    
    % --- 1C. SVD Spectral Extraction ---
    A_rev = zeros(N+1, N+1);
    for i = 1 : N+1
        for j = 1 : N+1
            if P_rev(i,j) > 0
                A_rev(i,j) = P_rev(i,j) * exp(0.5 * (log_pi(i) - log_pi(j)));
            end
        end 
    end
    s_rev = svd(eye(N+1) - A_rev); 
    tau_rev(idx) = 1 / s_rev(end-1);
    
    A_lift = zeros(2*N+2, 2*N+2);
    log_pi_lift = [log_pi; log_pi];
    for i = 1 : 2*N + 2
        for j = 1 : 2*N + 2
            if P_lift(i,j) > 0
                A_lift(i,j) = P_lift(i,j) * exp(0.5 * (log_pi_lift(i) - log_pi_lift(j))); 
            end
        end
    end
    s_lift = svd(eye(2*N+2) - A_lift); 
    tau_lift(idx) = 1 / s_lift(end-1);
end

% Empirical Power-Law Fit
p_rev  = polyfit(log(N_vec), log(tau_rev), 1);  exp_rev_emp  = p_rev(1);
p_lift = polyfit(log(N_vec), log(tau_lift), 1); exp_lift_emp = p_lift(1);

fprintf('\n--- NUMERICAL SCALING RESULTS ---\n');
fprintf('Reversible MH : tau ~ N^{%.3f} (Theoretical: N^{1.500})\n', exp_rev_emp);
fprintf('Lifted MH     : tau ~ N^{%.3f} (Theoretical: N^{0.750})\n\n', exp_lift_emp);

% --- FIGURE 1: ASYMPTOTIC RELAXATION TIME SCALING ---
figure('Color', 'w', 'Position', [100, 100, 700, 500]);
loglog(N_vec, tau_rev, 'o-', 'LineWidth', 2, 'Color', [0.8 0.2 0.2]); hold on;
loglog(N_vec, tau_lift, 's-', 'LineWidth', 2, 'Color', [0.2 0.4 0.8]);

% Theoretical reference lines
C_rev  = tau_rev(1)  / N_vec(1)^(3/2);
C_lift = tau_lift(1) / N_vec(1)^(3/4);
loglog(N_vec, C_rev * N_vec.^(3/2), '--', 'LineWidth', 1.4, 'Color', [0.8 0.2 0.2]);
loglog(N_vec, C_lift * N_vec.^(3/4), '--', 'LineWidth', 1.4, 'Color', [0.2 0.4 0.8]);
grid on;
xlabel('System Size $N$', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('Relaxation Time $\tau$', 'Interpreter', 'latex', 'FontSize', 13);
legend(sprintf('Reversible MH (Empirical: $N^{%.2f}$)', exp_rev_emp), ...
       sprintf('Lifted MH (Empirical: $N^{%.2f}$)', exp_lift_emp), ...
       'Theoretical $N^{3/2}$', 'Theoretical $N^{3/4}$', ...
       'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', 11);
title('Critical Curie-Weiss ($\beta = 1$): Relaxation Time Scaling', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'FontSize', 11);

%% PART 2: MEAN SQUARED ERROR DECAY OVER TIME (N = 320)
fprintf('2. Computing temporal evolution for N = 320...\n');
N = 320; 
k = (0:N)'; 
E = - (beta / (2*N)) * (2*k - N).^2;

% Construct transition matrices for N = 320
P_rev_320 = zeros(N+1, N+1);
for i = 1:(N+1)
    k_val = i - 1;
    if k_val < N, P_rev_320(i, i+1) = ((N - k_val) / N) * min(1, exp(-beta * (E(i+1) - E(i)))); end
    if k_val > 0, P_rev_320(i, i-1) = (k_val / N) * min(1, exp(-beta * (E(i-1) - E(i)))); end
    P_rev_320(i, i) = 1 - sum(P_rev_320(i, :));
end

P_lift_320 = zeros(2*N+2, 2*N+2);
for i = 1:(N+1)
    if i < N+1, T_plus = P_rev_320(i, i+1); else, T_plus = 0; end
    if i > 1,   T_minus = P_rev_320(i, i-1); else, T_minus = 0; end
    T_plus_minus = max(0, T_minus - T_plus); 
    T_minus_plus = max(0, T_plus - T_minus);
    
    if i < N+1, P_lift_320(i, i+1) = T_plus; end
    P_lift_320(i, i + N+1) = T_plus_minus; 
    P_lift_320(i, i)       = 1 - T_plus - T_plus_minus;
    
    if i > 1, P_lift_320(i + N+1, i - 1 + N+1) = T_minus; end
    P_lift_320(i + N+1, i)       = T_minus_plus; 
    P_lift_320(i + N+1, i + N+1) = 1 - T_minus - T_minus_plus;
end

M_vec = (2*k - N) / N; 
M_lift_vec = [M_vec; M_vec]; 

% Initialization: fully polarized state (m = +1)
v_rev  = zeros(1, N+1);   v_rev(N+1)  = 1.0; 
v_lift = zeros(1, 2*N+2); v_lift(N+1) = 0.5; v_lift(2*N+2) = 0.5;

T_max = 5000;
MSE_rev  = zeros(T_max, 1); 
MSE_lift = zeros(T_max, 1);

for t = 1:T_max
    v_rev  = v_rev * P_rev_320;
    v_lift = v_lift * P_lift_320;
    MSE_rev(t)  = (v_rev * M_vec)^2;
    MSE_lift(t) = (v_lift * M_lift_vec)^2;
end

% --- FIGURE 2: TEMPORAL ERROR DECAY ---
figure('Color', 'w', 'Position', [750, 100, 600, 450]);
semilogy(1:T_max, MSE_rev, 'r-', 'LineWidth', 2); hold on;
semilogy(1:T_max, MSE_lift, 'b-', 'LineWidth', 2); grid on;
xlabel('MCMC Steps (t)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Squared Expectation $\langle m(t) \rangle^2$', 'Interpreter', 'latex', 'FontSize', 12);
legend('Reversible MH (Diffusive Decay)', 'Lifted MH (Damped Oscillations)', 'FontSize', 11);
title(sprintf('Convergence to Equilibrium (N = %d)', N), 'FontSize', 13);

%% PART 3: SCHEMATIC DIAGRAM OF REVERSIBLE VS. LIFTED DYNAMICS
figure('Color', 'w', 'Position', [100, 550, 1100, 450]);
n_nodes = 9;
x = linspace(-1, 1, n_nodes);

% Left Subplot: Reversible
subplot(1,2,1); hold on; axis equal; axis off;
plot(x, zeros(size(x)), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
for i = 1:n_nodes-1
    quiver(x(i), 0, x(i+1)-x(i), 0, 0, 'k', 'LineWidth', 1.5, 'MaxHeadSize', 0.25);
    quiver(x(i+1), 0, x(i)-x(i+1), 0, 0, 'k', 'LineWidth', 1.5, 'MaxHeadSize', 0.25);
end
text(0, 0.28, 'Reversible Dynamics', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
text(-1, -0.22, '$-1$', 'Interpreter', 'latex', 'HorizontalAlignment', 'center', 'FontSize', 12);
text(1, -0.22, '$+1$', 'Interpreter', 'latex', 'HorizontalAlignment', 'center', 'FontSize', 12);
text(0, -0.42, 'Random walk behavior (No directional persistence)', 'HorizontalAlignment', 'center', 'FontSize', 11);

% Right Subplot: Lifted
subplot(1,2,2); hold on; axis equal; axis off;
y_plus  = 0.35;
y_minus = -0.35;
plot(x, y_plus*ones(size(x)), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
plot(x, y_minus*ones(size(x)), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 8);

for i = 1:n_nodes-1
    quiver(x(i), y_plus, x(i+1)-x(i), 0, 0, 'k', 'LineWidth', 1.5, 'MaxHeadSize', 0.25);
    quiver(x(i+1), y_minus, x(i)-x(i+1), 0, 0, 'k', 'LineWidth', 1.5, 'MaxHeadSize', 0.25);
end

mid = ceil(n_nodes/2);
quiver(x(mid), y_plus, 0, y_minus-y_plus, 0, 'k', 'LineWidth', 1.5, 'MaxHeadSize', 0.35);
quiver(x(mid), y_minus, 0, y_plus-y_minus, 0, 'k', 'LineWidth', 1.5, 'MaxHeadSize', 0.35);

text(0, 0.75, 'Lifted Dynamics', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
text(-1, y_plus+0.18, '$\xi = +1$', 'Interpreter', 'latex', 'HorizontalAlignment', 'center', 'FontSize', 12);
text(-1, y_minus-0.18, '$\xi = -1$', 'Interpreter', 'latex', 'HorizontalAlignment', 'center', 'FontSize', 12);
text(0, -0.78, 'Ballistic persistence + directional flips', 'HorizontalAlignment', 'center', 'FontSize', 11);
sgtitle('State-Space Topologies: Reversible vs. Lifted MCMC', 'FontSize', 15, 'FontWeight', 'bold');

%% PART 4: STOCHASTIC TRAJECTORIES (ZIG-ZAG PDMP LIMIT)
fprintf('3. Simulating stochastic trajectories...\n');
T_sim = 800; 
m_rev_traj  = zeros(T_sim, 1);
m_lift_traj = zeros(T_sim, 1);
xi_lift_traj = zeros(T_sim, 1); 

% Start from equilibrium (m = 0)
state_rev  = N/2 + 1; 
state_lift = N/2 + 1; 

for t = 1:T_sim
    % Reversible step
    p_cum_rev = cumsum(P_rev_320(state_rev, :));
    state_rev = find(rand <= p_cum_rev, 1);
    m_rev_traj(t) = M_vec(state_rev); 
    
    % Lifted step
    p_cum_lift = cumsum(P_lift_320(state_lift, :));
    state_lift = find(rand <= p_cum_lift, 1);
    
    if state_lift > (N + 1)
        k_idx = state_lift - (N + 1);
        xi_lift_traj(t) = -1; 
    else
        k_idx = state_lift;
        xi_lift_traj(t) = +1; 
    end
    m_lift_traj(t) = M_vec(k_idx);
end

switch_idx = find(diff(xi_lift_traj) ~= 0);

% --- FIGURE 4: TRAJECTORY COMPARISON ---
figure('Color', 'w', 'Position', [100, 100, 1100, 450]);
subplot(1,2,1);
plot(1:T_sim, m_rev_traj, 'Color', [0.8 0.2 0.2 0.9], 'LineWidth', 1.2);
title('Reversible MH Trajectory', 'FontSize', 13);
xlabel('Time (steps)', 'FontSize', 12); ylabel('Magnetization m(t)', 'FontSize', 12);
ylim([-0.45 0.45]); xlim([0 400]); grid on;
text(200, -0.38, 'Diffusive Dynamics', 'FontSize', 11, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'Color', [0.6 0.1 0.1]);

subplot(1,2,2);
plot(1:T_sim, m_lift_traj, 'Color', [0.2 0.4 0.8], 'LineWidth', 1.5); hold on;
plot(switch_idx, m_lift_traj(switch_idx), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 5);
title('Lifted MH Trajectory (Zig-Zag Process)', 'FontSize', 13);
xlabel('Time (steps)', 'FontSize', 12);
ylim([-0.45 0.45]); xlim([0 400]); grid on;
legend('Trajectory', 'Direction Flip (\xi \rightarrow -\xi)', 'Location', 'northwest', 'FontSize', 10);
text(200, -0.38, 'Piecewise Deterministic Motion (PDMP)', 'FontSize', 11, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'Color', [0.1 0.2 0.6]);
sgtitle('Phase Space Exploration (\beta = 1.0, N = 320)', 'FontSize', 14, 'FontWeight', 'bold');

%% PART 5: EIGENVALUE SPECTRUM IN THE COMPLEX PLANE
fprintf('4. Computing complex eigenvalue spectrum...\n');
lambda_rev  = eig(P_rev_320);
lambda_lift = eig(P_lift_320);

% --- FIGURE 5: COMPLEX SPECTRUM ---
figure('Color', 'w', 'Position', [750, 100, 600, 500]);
plot(real(lambda_rev), imag(lambda_rev), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6); hold on;
plot(real(lambda_lift), imag(lambda_lift), 'b.', 'MarkerSize', 12);
grid on; axis equal;
theta = linspace(0, 2*pi, 100); 
plot(cos(theta), sin(theta), 'k--', 'LineWidth', 1);
xlabel('Real Part', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Imaginary Part', 'FontSize', 12, 'FontWeight', 'bold');
title('Transition Matrix Eigenvalue Spectrum', 'FontSize', 13);
legend('Reversible MH (Real Spectrum)', 'Lifted MH (Complex Spectrum)', 'Unit Circle', 'Location', 'west');
fprintf('Analysis complete!\n');he
