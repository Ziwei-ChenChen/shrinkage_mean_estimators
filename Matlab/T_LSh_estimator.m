function theta_T_LSh = T_LSh_estimator(x)
    [Tn, N_T] = size(x);
    x_bar = mean(x);
    a1T_hat = sum(x_bar.^2) / N_T;
    trace_Sn = trace(cov(x));
    a2T_hat = a1T_hat + trace_Sn / (Tn * N_T);
    s_i = sum(x, 2);
    term1 = sum(s_i.^2);
    term2 = (sum(s_i)^2 - term1) / (Tn - 1);
    d2_hat = (term1 - term2) / (Tn^2 * N_T^2);
    grand_mean = mean(x_bar);
    d3_hat = sum((x_bar - grand_mean).^2) / N_T;
    d4_hat = grand_mean^2;
    tilde_a1_hat = a2T_hat - a1T_hat - d2_hat;
    tilde_a2_hat = tilde_a1_hat + d3_hat;
    delta_star_hat = tilde_a1_hat / tilde_a2_hat;
    xi_star_hat = 1 - (d2_hat / (d2_hat + d4_hat)) / delta_star_hat;
    theta_T_LSh = (1 - delta_star_hat) * x_bar + delta_star_hat * xi_star_hat * grand_mean;
    theta_T_LSh = theta_T_LSh(:);
end


