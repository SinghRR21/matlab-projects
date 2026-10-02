function Pout = ana_outage_ideal(N, thr, P)
% Gamma (moment-matching) approximation for ideal RIS phase shifts.
    b   = P.beta_sr * P.beta_rd;
    mu  = sqrt(pi*P.beta_d)/2 + N * (pi/4) * sqrt(b);
    vr  = P.beta_d*(1 - pi/4) + N * (1 - pi^2/16) * b;
    k   = mu^2 / vr;                                  % shape
    t   = vr / mu;                                    % scale
    Pout = gammainc(sqrt(thr) / t, k);                % Pr{A < sqrt(thr)}
end
