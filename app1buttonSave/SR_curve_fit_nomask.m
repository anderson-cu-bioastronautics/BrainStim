function J = SR_curve_fit_nomask(x, levels, SR, SE)

% Extract the previous iteration's returned values
A0_hat = x(1);
l_hat = x(2);
w0_hat = x(3);
f_hat = x(4);
B_hat = x(5);

% levels = levels';

% Generate the SR curve using these data
r = (l_hat./(sqrt(2)*pi)).*exp(-l_hat^2./(2*((levels-f_hat).^2)));
F = B_hat - (A0_hat.*l_hat./(levels - f_hat).^2).*r./sqrt(4.*r.^2 + w0_hat^2).*(levels-f_hat).*(levels >=f_hat);

% The cost function
J = sum((1./(SE.^2)).*(SR- F).^2);

end