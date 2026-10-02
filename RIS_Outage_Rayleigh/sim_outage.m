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
