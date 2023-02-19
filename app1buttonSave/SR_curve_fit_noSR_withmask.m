function J = SR_curve_fit_noSR_withmask(x,X,SR,SE)

% Extract this iteration's returned values
m_hat = x(1);
s_hat = x(2);
B_hat = x(3);

% Generate the SR curve using these data
F = B_hat + m_hat.*(X - s_hat).*(X >= s_hat);

% The cost function
J = sum((1./(SE.^2)).*(SR - F).^2);
end