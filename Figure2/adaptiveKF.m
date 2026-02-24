function alpha_d = adaptiveKF(alpha_s, P0, L)
% ADAPTIVEKF  Adaptive Kalman Filter
T        = numel(alpha_s);
alpha_d  = zeros(1,T);   alpha_d(1)=alpha_s(1);
P_vals   = zeros(1,T);   P_vals(1)=P0;

rho_values = logspace(-2,2,100);  % grid of candidate rhos

for t = 2:T
  % 1) define window
  w0  = max(1, t-L+1);
  idx = w0:t;
  ys  = alpha_s(idx);

  % 2) fit linear trend 
  p    = polyfit(idx, ys, 1);
  y_fit= polyval(p, idx);

  % 3) estimate R from residuals
  res  = ys - y_fit;
  R_est= var(res);

  % 4) search best rho by restarting filter at w0 for each trial
  bestSSE = Inf;
  bestRho = rho_values(1);

  for rho = rho_values
    Q = rho * R_est;
    % --- restart at window start ---
    if w0==1
      a_tmp = alpha_s(1);
      P_tmp = P0;
    else
      a_tmp = alpha_d(w0-1);
      P_tmp = P_vals(w0-1);
    end

    % run KF over the window, measure SSE vs trend
    errs = zeros(size(idx));
    for k = 1:numel(idx)
      % prediction
      P_tmp = P_tmp + Q;
      K     = P_tmp/(P_tmp + R_est);
      % update with static obs
      a_tmp = a_tmp + K*(ys(k) - a_tmp);
      P_tmp = (1-K)*P_tmp;
      % error to linear trend
      errs(k) = a_tmp - y_fit(k);
    end

    SSE = sum(errs.^2);
    if SSE < bestSSE
      bestSSE = SSE;
      bestRho = rho;
    end
  end

  % Store best Q
  Q_est = bestRho * R_est;

  % 5) final single KF step at t, **chained** from t-1
  P_pred       = P_vals(t-1) + Q_est;
  K_final      = P_pred/(P_pred + R_est);
  alpha_d(t)   = alpha_d(t-1) + K_final*(alpha_s(t) - alpha_d(t-1));
  P_vals(t)    = (1-K_final)*P_pred;
end
end
