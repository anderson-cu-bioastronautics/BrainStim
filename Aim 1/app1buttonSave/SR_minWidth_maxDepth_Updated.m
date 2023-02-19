function [c, ceq] = SR_minWidth_maxDepth_Updated(x)
global TestType SRType levels thisSubjectMeasuredThresholds method widthThreshold1 widthThreshold2

if (TestType == 2 && SRType == 1)
    % Extract the fitted values and generate the SR curve
    A0_hat = x(1);
    l_hat = x(2);
    w0_hat = x(3);
    m_hat = x(4);
    s_hat = x(5);
    f_hat = x(6);
    B_hat = x(7);
    
else
    A0_hat = x(1);
    l_hat = x(2);
    w0_hat = x(3);
    f_hat = x(4);
    B_hat = x(5);
end


% Generate the SR curve using these values
X = linspace(min(levels),max(levels)*3);
r = (l_hat./(sqrt(2)*pi)).*exp(-l_hat^2./(2*((X-f_hat).^2)));
F = B_hat - (A0_hat.*l_hat./(X - f_hat).^2).*r./sqrt(4.*r.^2 + w0_hat^2).*(X-f_hat).*(X >=f_hat);



%% The nonlinear constraints. Each constraint function in "c" must be less
% than or equal to zero.

% switch method
%     
%     case 'Derivative'
        
        % Take the derivative of the SR curve
        dy_dx = diff(F)./diff(X);
        
        % Find the x value of this derivative's minimum and maximum - roughly
        % correlates to the downsloping and upsloping portions of the dip
        x1 = min(X(find(dy_dx == min(dy_dx))));
        x2 = max(X(find(dy_dx == max(dy_dx))));
        
%     case 'Percentage'
%         
%         % The depth of the dip
%         depth = max(F) - min(F);
%         
%         % Prevent division by zero
%         if depth == 0
%             depth = 1e-6;
%         end
%         
%         % The curve value that is the lower threshold of the dip depth
%         Yval1 = min(F) + widthThreshold1*depth;
%         
%         % The curve value at 90% of the dip depth
%         Yval2 = min(F) + widthThreshold2*depth;
%         
%         % The x value at the minimum fit value
%         min_x = find(F == min(F));
%         
%         % Break the curve into the portions before and after the minimum
%         F1 = F(1:min_x);
%         F2 = F(min_x:end);
%         X1 = X(1:min_x);
%         X2 = X(min_x:end);
%         
%         % The value of X where the curve drops below 95% of the dip depth
%         x1 = X1(max(find(min(abs(F1 - Yval1)) == abs(F1-Yval1))));
%         
%         % The X value where the curve comes back up to 90% of the dip depth
%         x2 = X2(min(find(min(abs(F2 - Yval2)) == abs(F2-Yval2))));
%         
%     case 'Integral'
%         
%         % define x1 as the first place the derivative is less than 0
%         x1 = X(min(find(diff(F) < 0)));
%         
%         % If it doesn't exist, make x1 equal to zero
%         if isempty(x1)
%             x1 = 0;
%         end
%             
%         
%         % The depth of the curve
%         depth = max(F) - min(F);
%         
%         % Prevent division by zero
%         if depth == 0
%             depth = 1e-6;
%         end
%         
%         % Integrate the area under the curve (subtract the baseline
%         % threshold)
%         x2 = x1 -trapz(X,F-B_hat)./depth;
%         
%     case 'None'
%         
%         % No constraints; thererefore these make no sense
%         x1 = 0;
%         x2 = 0;
        
%end
if SRType == 1
    levelSpacing = 5;
else
    levelSpacing = 0.1;
end

% The minimum of the dip cannot go below 90% of the minimum observed
% threshold.

% switch method
%     case 'None'
%         
%         % Do not impose any constraints
%         c = [];
%         
%     otherwise
        
        c(1) = x2 - max(levels)*1.5;
        
        c(2) = 0.9*min(thisSubjectMeasuredThresholds) - min(F);
        
        c(3) = 0.5*levelSpacing + x1 - x2;
        
        c(4) = x1 - x2 - 5*levelSpacing;
        
% end

%% The ceq return has constraints that must be equal to zero. None exist
% here.
ceq = [];
end