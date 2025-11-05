
import numpy as np

# =============================================================================

def St_MSh_estimator(x):
    Tn, N_T = x.shape
    x_bar = np.mean(x, axis=0)
    a1T_hat = np.sum(x_bar**2) / N_T
    trace_Sn = np.trace(np.cov(x, rowvar=False))
    a2T_hat = a1T_hat + trace_Sn / (Tn * N_T)
    c_star_hat = a1T_hat / a2T_hat
    theta_St_MSh = c_star_hat * x_bar
    return theta_St_MSh

def D_MSh_estimator(x):
    Tn, N_T = x.shape
    x_bar = np.mean(x, axis=0)
    Sn = np.cov(x, rowvar=False)
    c_star_hat = (x_bar**2) / (x_bar**2 + (1 / Tn) * np.diag(Sn))
    theta_D_MSh = c_star_hat * x_bar
    return theta_D_MSh

def O_LSh_estimator(x):
    Tn, N_T = x.shape
    x_bar = np.mean(x, axis=0)
    a1T_hat = np.sum(x_bar**2) / N_T
    trace_Sn = np.trace(np.cov(x, rowvar=False))
    a2T_hat = a1T_hat + trace_Sn / (Tn * N_T)
    s_i = np.sum(x, axis=1)
    term1 = np.sum(s_i**2)
    term2 = (np.sum(s_i)**2 - term1) / (Tn - 1)
    d2_hat = (term1 - term2) / (Tn**2 * N_T**2)
    grand_mean = np.mean(x_bar)
    d3_hat = np.sum((x_bar - grand_mean)**2) / N_T
    tilde_a1_hat = a2T_hat - a1T_hat - d2_hat
    tilde_a2_hat = tilde_a1_hat + d3_hat
    delta_star_hat = tilde_a1_hat / tilde_a2_hat
    theta_O_LSh = (1 - delta_star_hat) * x_bar + delta_star_hat * grand_mean
    return theta_O_LSh

def T_LSh_estimator(x):
    Tn, N_T = x.shape
    x_bar = np.mean(x, axis=0)
    a1T_hat = np.sum(x_bar**2) / N_T
    trace_Sn = np.trace(np.cov(x, rowvar=False))
    a2T_hat = a1T_hat + trace_Sn / (Tn * N_T)
    s_i = np.sum(x, axis=1)
    term1 = np.sum(s_i**2)
    term2 = (np.sum(s_i)**2 - term1) / (Tn - 1)
    d2_hat = (term1 - term2) / (Tn**2 * N_T**2)
    grand_mean = np.mean(x_bar)
    d3_hat = np.sum((x_bar - grand_mean)**2) / N_T
    d4_hat = grand_mean**2
    tilde_a1_hat = a2T_hat - a1T_hat - d2_hat
    tilde_a2_hat = tilde_a1_hat + d3_hat
    delta_star_hat = tilde_a1_hat / tilde_a2_hat
    xi_star_hat = 1 - (d2_hat / (d2_hat + d4_hat)) / delta_star_hat
    theta_T_LSh = (1 - delta_star_hat) * x_bar + delta_star_hat * xi_star_hat * grand_mean
    return theta_T_LSh


