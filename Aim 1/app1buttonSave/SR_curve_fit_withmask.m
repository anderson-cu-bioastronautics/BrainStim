function J = SR_curve_fit_withmask(x, X, SR, SE)

% Extract the previous iteration's returned values
A0_hat = x(1);
l_hat = x(2);
w0_hat = x(3);
m_hat = x(4);
s_hat = x(5);
f_hat = x(6);
B_hat = x(7);

% Generate the SR curve using these data
r = (l_hat./(sqrt(2)*pi)).*exp(-l_hat^2./(2*((X-f_hat).^2)));
F = B_hat - (A0_hat.*l_hat./(X - f_hat).^2).*r./sqrt(4.*r.^2 + w0_hat^2).*(X-f_hat).*(X >=f_hat);
mask = m_hat*(X-s_hat).*(X >= s_hat);
F = F + mask;

% The cost function
J = sum((1./(SE.^2)).*(SR- F).^2);

end