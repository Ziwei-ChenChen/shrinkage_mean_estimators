function theta_D_MSh = D_MSh_estimator(x)
    [Tn, ~] = size(x);
    x_bar = mean(x);
    Sn = cov(x);
    c_star_hat = (x_bar.^2) ./ (x_bar.^2 + diag(Sn)' / Tn);
    theta_D_MSh = c_star_hat .* x_bar;
    theta_D_MSh = theta_D_MSh(:);
end


