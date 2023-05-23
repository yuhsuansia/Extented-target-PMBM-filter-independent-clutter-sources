function lik_intensity = clutter_lik(z,stationary_clutter)
%CLUTTER_LIK computes the sum of log-likelihood of Poisson intensity of stationary source

lik_intensity = size(z,2)*log(stationary_clutter.lambda_clutter) + sum(log_mvnpdf(z',stationary_clutter.x_clutter',stationary_clutter.X_clutter));

end