function theta_St_MSh = St_MSh_estimator(x)
    [Tn, N_T] = size(x);
    x_bar = mean(x);
    a1T_hat = sum(x_bar.^2) / N_T;
    trace_Sn = trace(cov(x));
    a2T_hat = a1T_hat + trace_Sn / (Tn * N_T);
    c_star_hat = a1T_hat / a2T_hat;
    theta_St_MSh = c_star_hat * x_bar;
    theta_St_MSh = theta_St_MSh(:);
end


