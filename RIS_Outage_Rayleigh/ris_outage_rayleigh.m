%% ========================================================================
%  RIS-Assisted Outage Probability under Rayleigh Fading (Monte Carlo)
%  ========================================================================
%  System model (single-antenna source S, destination D, N-element RIS):
%
%      y = sqrt(Pt) * ( h_d + sum_i  g_i * exp(j*theta_i) * h_i ) * x + n
%
%      h_i ~ CN(0, beta_sr)   : S  -> RIS element i   (Rayleigh)
%      g_i ~ CN(0, beta_rd)   : RIS element i -> D    (Rayleigh)
%      h_d ~ CN(0, beta_d)    : S  -> D direct link   (Rayleigh, 0 = blocked)
%      n   ~ CN(0, sigma^2)
%
%  Instantaneous SNR:  gamma = rho * |h_d + sum_i h_i g_i e^{j theta_i}|^2,
%  with rho = Pt / sigma^2 (average transmit SNR).
%
%  Outage probability:  P_out = Pr{ log2(1 + gamma) < R_th }
%                             = Pr{ gamma < gamma_th },  gamma_th = 2^R_th - 1
%
%  RIS phase configurations simulated:
%    'ideal'  : continuous optimal phases, theta_i = arg(h_d) - arg(h_i g_i)
%               -> coherent combining, |.| = |h_d| + sum |h_i||g_i|
%    'quant'  : optimal phases quantized to b bits (practical RIS)
%    'random' : uniformly random phases (no RIS optimisation)
%    'none'   : no RIS, direct link only
%
%  Analytical benchmarks (plotted as lines):
%    * Ideal phases: moment-matched Gamma approximation of
%          A = |h_d| + sum |h_i||g_i|
%      E[|h_i||g_i|] = (pi/4) sqrt(beta_sr beta_rd),
%      Var[|h_i||g_i|] = (1 - pi^2/16) beta_sr beta_rd,
%      E[|h_d|] = sqrt(pi beta_d)/2,  Var[|h_d|] = (1 - pi/4) beta_d
%      A ~ Gamma(k, t), k = mu^2/var, t = var/mu
%      P_out = Pr{A < sqrt(gamma_th/rho)} = gammainc(sqrt(gamma_th/rho)/t, k)
%    * Random phases: CLT -> received coefficient ~ CN(0, beta_d + N beta_sr beta_rd)
%      P_out = 1 - exp(-gamma_th / (rho (beta_d + N beta_sr beta_rd)))
%    * Direct link only (exact): P_out = 1 - exp(-gamma_th / (rho beta_d))
%
%  Requires MATLAB R2016b+ (local functions in scripts, implicit expansion).
%  ========================================================================

clear; clc; close all;
rng(2026);                                 % reproducibility

%% ---------------------------- Parameters --------------------------------
P.N_list   = [8 16 32 64];                 % number of RIS elements
P.rho_dB   = -40:1:20;                     % average transmit SNR Pt/sigma^2 [dB]
P.R_th     = 2;                            % target rate [bit/s/Hz]
P.beta_sr  = 1;                            % mean power gain S  -> RIS (per element)
P.beta_rd  = 1;                            % mean power gain RIS -> D  (per element)
P.beta_d   = 1;                            % mean power gain S  -> D  (0 = blocked)
P.nTrials  = 1e6;                          % Monte Carlo channel realisations
P.batch    = 5e4;                          % realisations per batch (memory control)

P.N_cmp    = 32;                           % N used for phase-configuration comparison
P.bits_cmp = [1 2];                        % phase-quantisation resolutions to compare
P.rho_fix_dB = -15;                        % fixed SNR for P_out vs N sweep
P.N_sweep  = 1:2:80;                       % RIS sizes for P_out vs N sweep

rho      = 10.^(P.rho_dB/10);
gamma_th = 2^P.R_th - 1;
thr      = gamma_th ./ rho;                % outage <=> |coef|^2 < gamma_th/rho

fprintf('RIS outage simulation: R_th = %g bps/Hz (gamma_th = %.3f), %g trials\n', ...
        P.R_th, gamma_th, P.nTrials);

%% ------------------ Fig. 1: P_out vs SNR for several N ------------------
nN   = numel(P.N_list);
Psim = zeros(nN, numel(rho));
Pana = zeros(nN, numel(rho));
for n = 1:nN
    N = P.N_list(n);
    fprintf('  ideal phases, N = %2d ...\n', N);
    Psim(n,:) = sim_outage(N, 'ideal', Inf, thr, P);
    Pana(n,:) = ana_outage_ideal(N, thr, P);
end
Pdir_sim = sim_outage(0, 'none', Inf, thr, P);
Pdir_ana = 1 - exp(-thr / P.beta_d);

figure('Name','Outage vs SNR','Color','w'); hold on; grid on; box on;
cols = lines(nN);
leg  = {};
for n = 1:nN
    semilogy(P.rho_dB, Psim(n,:), 'o', 'Color', cols(n,:), 'LineWidth', 1.2);
    leg{end+1} = sprintf('Sim., N = %d', P.N_list(n));        %#ok<SAGROW>
end
for n = 1:nN
    semilogy(P.rho_dB, Pana(n,:), '-', 'Color', cols(n,:), 'LineWidth', 1.2);
    leg{end+1} = sprintf('Gamma approx., N = %d', P.N_list(n)); %#ok<SAGROW>
end
semilogy(P.rho_dB, Pdir_sim, 'ks', 'LineWidth', 1.2);
semilogy(P.rho_dB, Pdir_ana, 'k--', 'LineWidth', 1.2);
leg = [leg, {'Sim., no RIS (direct)', 'Exact, no RIS (direct)'}];
set(gca, 'YScale', 'log');
ylim([1/P.nTrials 1]); xlim([P.rho_dB(1) P.rho_dB(end)]);
xlabel('Average transmit SNR  \rho = P_t/\sigma^2  (dB)');
ylabel('Outage probability  P_{out}');
title(sprintf('RIS-assisted outage, Rayleigh fading, R_{th} = %g bps/Hz', P.R_th));
legend(leg, 'Location', 'southwest', 'NumColumns', 2);

%% --------- Fig. 2: phase configuration comparison at fixed N ------------
N = P.N_cmp;
fprintf('  phase-configuration comparison, N = %d ...\n', N);
P_ideal = sim_outage(N, 'ideal', Inf, thr, P);
P_quant = zeros(numel(P.bits_cmp), numel(rho));
for b = 1:numel(P.bits_cmp)
    P_quant(b,:) = sim_outage(N, 'quant', P.bits_cmp(b), thr, P);
end
P_rand     = sim_outage(N, 'random', Inf, thr, P);
P_rand_ana = 1 - exp(-thr / (P.beta_d + N*P.beta_sr*P.beta_rd));

figure('Name','Phase configurations','Color','w'); hold on; grid on; box on;
semilogy(P.rho_dB, P_ideal, 'b-o', 'LineWidth', 1.2);
mk = {'r-^', 'm-v', 'c-d', 'g-p'};
leg = {'Ideal (continuous) phases'};
for b = 1:numel(P.bits_cmp)
    semilogy(P.rho_dB, P_quant(b,:), mk{mod(b-1,numel(mk))+1}, 'LineWidth', 1.2);
    leg{end+1} = sprintf('%d-bit quantised phases', P.bits_cmp(b)); %#ok<SAGROW>
end
semilogy(P.rho_dB, P_rand, 'ks', 'LineWidth', 1.2);
semilogy(P.rho_dB, P_rand_ana, 'k-', 'LineWidth', 1.0);
semilogy(P.rho_dB, Pdir_sim, 'x--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
leg = [leg, {'Random phases (sim.)', 'Random phases (CLT)', 'No RIS (direct)'}];
set(gca, 'YScale', 'log');
ylim([1/P.nTrials 1]); xlim([P.rho_dB(1) P.rho_dB(end)]);
xlabel('Average transmit SNR  \rho  (dB)');
ylabel('Outage probability  P_{out}');
title(sprintf('Effect of RIS phase configuration, N = %d', N));
legend(leg, 'Location', 'southwest');

%% ------------------- Fig. 3: P_out vs number of elements ----------------
thr_fix = gamma_th / 10^(P.rho_fix_dB/10);
fprintf('  P_out vs N at rho = %g dB ...\n', P.rho_fix_dB);
PN_sim = zeros(size(P.N_sweep));
PN_ana = zeros(size(P.N_sweep));
for n = 1:numel(P.N_sweep)
    PN_sim(n) = sim_outage(P.N_sweep(n), 'ideal', Inf, thr_fix, P);
    PN_ana(n) = ana_outage_ideal(P.N_sweep(n), thr_fix, P);
end

figure('Name','Outage vs N','Color','w'); hold on; grid on; box on;
semilogy(P.N_sweep, PN_sim, 'bo', 'LineWidth', 1.2);
semilogy(P.N_sweep, PN_ana, 'b-', 'LineWidth', 1.2);
set(gca, 'YScale', 'log');
ylim([1/P.nTrials 1]);
xlabel('Number of RIS elements  N');
ylabel('Outage probability  P_{out}');
title(sprintf('Outage vs RIS size, \\rho = %g dB, R_{th} = %g bps/Hz', ...
              P.rho_fix_dB, P.R_th));
legend('Simulation (ideal phases)', 'Gamma approximation', 'Location', 'southwest');

%% ------------------------------ Summary ---------------------------------
fprintf('\nSNR (dB) needed for P_out <= 1e-3 (ideal phases, simulation):\n');
for n = 1:nN
    idx = find(Psim(n,:) <= 1e-3, 1);
    if isempty(idx)
        fprintf('  N = %2d : not reached in simulated range\n', P.N_list(n));
    else
        fprintf('  N = %2d : %g dB\n', P.N_list(n), P.rho_dB(idx));
    end
end

%% =========================== Local functions ============================
function Pout = sim_outage(N, mode, bits, thr, P)
% Monte Carlo outage probability for each threshold in thr (= gamma_th/rho).
    cnt  = zeros(1, numel(thr));
    left = P.nTrials;
    while left > 0
        nB = min(P.batch, left);
        left = left - nB;

        hd = sqrt(P.beta_d/2) * (randn(nB,1) + 1j*randn(nB,1));
        if N == 0 || strcmp(mode, 'none')
            coef = hd;
        else
            h = sqrt(P.beta_sr/2) * (randn(nB,N) + 1j*randn(nB,N));
            g = sqrt(P.beta_rd/2) * (randn(nB,N) + 1j*randn(nB,N));
            c = h .* g;                                   % cascaded channels
            switch mode
                case 'ideal'
                    % Co-phase every path with the direct link
                    coef = abs(hd) + sum(abs(c), 2);
                case 'quant'
                    d     = 2*pi / 2^bits;                % phase step
                    theta = angle(hd) - angle(c);         % optimal phases
                    theta = d * round(theta / d);         % b-bit quantisation
                    coef  = hd + sum(c .* exp(1j*theta), 2);
                case 'random'
                    theta = 2*pi*rand(nB, N);
                    coef  = hd + sum(c .* exp(1j*theta), 2);
                otherwise
                    error('Unknown mode "%s".', mode);
            end
        end
        % Outage when |coef|^2 < gamma_th / rho  (nB x nRho comparison)
        cnt = cnt + sum(abs(coef).^2 < thr, 1);
    end
    Pout = cnt / P.nTrials;
end

function Pout = ana_outage_ideal(N, thr, P)
% Gamma (moment-matching) approximation for ideal RIS phase shifts.
    b   = P.beta_sr * P.beta_rd;
    mu  = sqrt(pi*P.beta_d)/2 + N * (pi/4) * sqrt(b);
    vr  = P.beta_d*(1 - pi/4) + N * (1 - pi^2/16) * b;
    k   = mu^2 / vr;                                  % shape
    t   = vr / mu;                                    % scale
    Pout = gammainc(sqrt(thr) / t, k);                % Pr{A < sqrt(thr)}
end
