function [slope,intercept,STAT]=myregr(x,y,varargin)
%MYREGR Perform a least-squares linear regression with rich output.
%This function computes a least-squares linear regression supplying several
%informative statistics and plots.
%
% Syntax:
%       myregr(x,y)
%       myregr(x,y,verbose)
%      
% Inputs:
%   X - Array of the independent variable (row or column vector).
%   Y - Dependent variable. If Y is a matrix, the i-th Y row is a
%       repeated measure of i-th X point. The mean value will be used.
%   verbose - Flag to display all information and plots (default = 1).
%
% Outputs:
%   - Slope with standard error and 95% C.I.
%   - Intercept with standard error and 95% C.I.
%   - Pearson's correlation coefficient with 95% C.I. and its adjusted form
%   - Spearman's correlation coefficient
%   - Regression Standard Error
%   - Total variability, regression variability, residual variability
%   - Student's t-test on slope (H0: slope = 0)
%   - Student's t-test on intercept (H0: intercept = 0)
%   - Modified Levene's test for residual homoscedasticity
%   - Power of the regression
%   - "Reverse" regression (Y as independent variable)
%   - Regression on principal standardized component (when meaningful)
%   - Deming regression
%   - Plots:
%       o Data points + LS regression line
%       o 95% CI for regression and for a new observation
%       o Residual plot
%       o Comparison of different regressions
%         (X independent, Y independent, principal component, Deming)
%         with vertical distances (OLS) and perpendicular distances (Deming)
%
% Example:
%   x = [1.0 2.3 3.1 4.8 5.6 6.3];
%   y = [2.6 2.8 3.1 4.7 4.1 5.3];
%   myregr(x,y)
%
% SEE ALSO: myregrinv, myregrcomp
%
% Created by Giuseppe Cardillo
% giuseppe.cardillo.75@gmail.com
% GitHub repository: https://github.com/dnafinder/myregression
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2007) MyRegression: a simple function on LS linear
% regression with many informative outputs.
% GitHub repository: https://github.com/dnafinder/myregression
%
% This code is released under the GNU GPL-3.0 license.

%% Input error handling
p = inputParser;
% Allow X to be any numeric vector (row or column)
addRequired(p,'x',@(x) validateattributes(x,{'numeric'}, ...
    {'vector','real','finite','nonnan','nonempty'}));
% Y can be a vector or a 2D matrix
addRequired(p,'y',@(x) validateattributes(x,{'numeric'}, ...
    {'2d','real','finite','nonnan','nonempty'}));
% Verbose flag: 0 or 1
addOptional(p,'verbose',1, @(x) isnumeric(x) && isreal(x) && ...
    isfinite(x) && isscalar(x) && (x==0 || x==1));

parse(p,x,y,varargin{:});
verbose = p.Results.verbose;
alpha   = 0.05;
clear p

%% Prepare X and Y as column vectors
x = x(:);                         % force column vector
if isvector(y)
    yt = y(:);                    % simple column vector
else
    yt = mean(y)';                % mean over rows if repeated measures
end
assert(length(x)==length(yt), ...
    'X and Y must have the same number of elements.');

%% Collapse duplicated X values (average the corresponding Y values)
ux = unique(x);
if length(ux)~=length(x)
    uy = zeros(size(ux));
    for I = 1:length(ux)
        c = sum(x==ux(I));
        if c==1
            uy(I) = yt(x==ux(I));
        else
            uy(I) = mean(yt(x==ux(I)));
        end
    end
    x  = ux(:);
    yt = uy(:);
    clear uy I
end
clear ux

%% Build design matrix for regress (X and intercept)
xtmp = [x ones(length(x),1)];
ytmp = yt;

%% Least-squares regression coefficients using MATLAB's regress
[p,pINT,R,Rint] = regress(ytmp,xtmp);

%% Outlier detection based on confidence intervals of residuals
% Outliers are points whose residual intervals do not include zero
outl  = find(ismember(sign(Rint),[-1 1],'rows')==0);
reply = ''; % ensure reply is always defined

if ~isempty(outl)
    if verbose==1
        % Interactive mode: show outliers and ask user
        disp('These points are outliers at 95% fiducial level')
        disp(array2table([xtmp(outl) ytmp(outl)], ...
            'VariableNames',{'X' 'Y'}))
        reply = input('Do you want to delete outliers? Y/N [Y]: ', 's');
        disp(' ')
    else
        % Silent mode: by default remove outliers (no user interaction)
        reply = 'Y';
    end
    
    if isempty(reply) || upper(reply)=='Y'
        % Remove outliers from data and recompute regression
        ytmp(outl)   = [];
        xtmp(outl,:) = [];
        [p,pINT,R]   = regress(ytmp,xtmp);
    end
end

%% Remove intercept column from xtmp (keep only X for later computations)
xtmp(:,2) = [];

%% Save slope (m) and intercept (q) from LS regression
m(1) = p(1);   % slope
q(1) = p(2);   % intercept

n   = length(xtmp);       % number of points used in regression
xm  = mean(xtmp);         % mean of X
xsd = std(xtmp);          % standard deviation of X

%% Standard error of regression coefficients
% Student's t critical value for 95% CI
if isvector(y)
    cv = tinv(0.975,n-2);
else
    % If Y is a matrix, degrees of freedom are adjusted
    cv = tinv(0.975,sum(size(y))-3);
end

% Slope standard error and confidence interval
m(2) = (pINT(3)-p(1))/cv;    % std error of slope
m    = [m pINT(1,:)];        % add lower and upper 95% CI

% Intercept standard error and confidence interval
q(2) = (pINT(4)-p(2))/cv;    % std error of intercept
q    = [q pINT(2,:)];        % add lower and upper 95% CI

% Build output structures
slope.value = m(1);
slope.se    = m(2);
slope.lv    = m(3);
slope.uv    = m(4);

intercept.value = q(1);
intercept.se    = q(2);
intercept.lv    = q(3);
intercept.uv    = q(4);

%% Pearson correlation coefficient and adjusted version
[rp,pr,rlo,rup] = corrcoef(xtmp,ytmp);
r(1) = rp(2);                             % Pearson r
r(2) = realsqrt((1-r(1)^2)/(n-2));        % std error of r
r(3) = rlo(2);                            % lower CI
r(4) = rup(2);                            % upper CI
% Adjusted Pearson correlation (removing small-sample bias)
r(5) = sign(r(1))*(abs(r(1))-((1-abs(r(1)))/(n-2)));

%% Spearman correlation coefficient
rx = tiedrank(xtmp);
ry = tiedrank(ytmp);
d  = rx-ry;
sp = 1-(6*sum(d.^2)/(n^3-n));
% Spearman CI via Fisher z-transform
rs = [sp NaN tanh(atanh(sp)+[-1 1].*(1.96/realsqrt(n-3))) NaN];

%% Variability components
% Total variability (around mean predicted value at xm)
ym   = polyval(p,xm);
vtot = sum((ytmp-ym).^2);

% Regression variability (explained by regression line)
ystar = ytmp-R;                  % fitted values yhat = y - residuals
vreg  = sum((ystar-ym).^2);

% Residual variability (unexplained)
vres = sum(R.^2);

%% Regression standard error (RSE)
if isvector(y)
    % Classical formula for simple linear regression
    RSE = realsqrt(vres/(n-2));
else
    % When Y is a matrix with repeated measures, take into account
    % within-row variability as in original implementation
    if ~isempty(outl) && (isempty(reply) || upper(reply)=='Y')
        y2       = y;
        y2(outl) = [];
        RSE = realsqrt((vres+sum(sum((y2-repmat(ytmp',size(y,1),1)).^2))) / ...
                       (sum(size(y2))-3));
    else
        RSE = realsqrt((vres+sum(sum((y-repmat(yt',size(y,1),1)).^2))) / ...
                       (sum(size(y))-3));
    end
end

%% Confidence intervals of regression and new observation
% CI of mean regression line
sy   = RSE*realsqrt(1/n+(((xtmp-xm).^2)/((n-1)*xsd^2)));
cir  = [ystar+cv*sy ystar-cv*sy];

% CI for a new observation (includes residual variance)
sy2  = realsqrt(sy.^2+RSE^2);
cir2 = [ystar+cv*sy2 ystar-cv*sy2];

%% STAT structure
STAT.rse = RSE;
STAT.cv  = cv;
STAT.n   = n;
STAT.xm  = mean(x);
STAT.ym  = ym;
STAT.sse = sum((xtmp-xm).^2);
STAT.r   = r;

%% Display results and compute additional regressions if verbose == 1
if verbose==1
    tr = repmat('-',1,80);

    % Basic statistics for "other regressions" and Deming
    n   = length(xtmp);
    mx  = mean(xtmp);
    my  = mean(ytmp);
    z   = cov(xtmp,ytmp);
    vx  = z(1);     % variance of X
    vy  = z(4);     % variance of Y
    lambda  = vx/vy;
    devx    = vx*(n-1);          % sum of squares of X (about the mean)
    devy    = vy*(n-1);          % sum of squares of Y (about the mean)
    covarxy = z(2);              % covariance
    codevxy = covarxy*(n-1);     % sum of cross-products

    %% Main regression summary
    disp(' ')
    disp('REGRESSION SETTING X AS INDEPENDENT VARIABLE')
    disp(tr)
    disp(array2table([m;q],'RowNames',{'Slope','Intercept'}, ...
        'VariableNames',{'Value','Standard_Error','Lower_bound','Upper_bound'}))
    
    fprintf('\t\t\tCorrelation Coefficients\n')
    disp(tr)
    disp(array2table([r;rs],'RowNames',{'Pearson','Spearman'}, ...
        'VariableNames',{'Value','Standard_Error','Lower_bound','Upper_bound','Adjusted'}))
    
    fprintf('\t\t\tVariability\n')
    disp(tr)
    disp(array2table([RSE vtot vreg vreg/vtot*100 vres vres/vtot*100], ...
        'RowNames',{'Regr_SE'}, ...
        'VariableNames',{'Value','Total','By_Regression','Percent1','Residual','Percent2'}))
    disp(' ')
 
    %% Student's t-tests on slope and intercept
    sturegr    = cell(3,5);
    % Slope
    sturegr{1,1} = abs(m(1)/m(2));     % t-statistic
    sturegr{1,2} = cv;                 % critical t
    sturegr{1,3} = pr(2);              % p-value from corrcoef
    if sturegr{1,1}>sturegr{1,2}
        sturegr{1,5} = 'slope ~= 0';
        sturegr{1,4} = 1 - tcdf(tinv(1-alpha,n-2) - sturegr{1,1},n-2);
    else
        sturegr{1,5} = 'slope = 0';
        sturegr{1,4} = tcdf(sturegr{1,1} - tinv(1-alpha,n-2),n-2);
        % m(1) could be set to 0 here for display, but we keep slope.value intact
    end
    % Intercept
    sturegr{2,1} = abs(q(1)/q(2));
    sturegr{2,2} = cv;
    sturegr{2,3} = (1-tcdf(sturegr{2,1},n-2))*2;
    if sturegr{2,1}>sturegr{2,2}
        sturegr{2,5} = 'intercept ~= 0';
        sturegr{2,4} = 1 - tcdf(tinv(1-alpha,n-2) - sturegr{2,1},n-2);
    else
        sturegr{2,5} = 'intercept = 0';
        sturegr{2,4} = tcdf(sturegr{2,1} - tinv(1-alpha,n-2),n-2);
        % q(1) could be set to 0 here for display, but intercept.value stays as LS estimate
    end
    
    %% Modified Levene's test for residual homoscedasticity
    xme = median(xtmp);
    e1  = R(xtmp<=xme);  me1 = median(e1);
    d1  = abs(e1-me1);   dm1 = mean(d1);  l1 = length(e1);
    e2  = R(xtmp>xme);   me2 = median(e2);
    d2  = abs(e2-me2);   dm2 = mean(d2);  l2 = length(e2);
    gl  = (l1+l2-2);
    S2p = (sum((d1-dm1).^2)+sum((d2-dm2).^2))/gl;
    sturegr{3,1} = abs(dm1-dm2)/realsqrt(S2p*(1/l1+1/l2));
    sturegr{3,2} = tinv(1-alpha/2,gl);
    sturegr{3,3} = (1-tcdf(sturegr{3,1},gl));
    sturegr{3,4} = NaN;
    if sturegr{3,1}>sturegr{3,2}
        sturegr{3,5} = 'Heteroschedastic';
    else
        sturegr{3,5} = 'Homoschedastic';
    end
    
    disp('Statistical tests')
    disp(tr)
    disp(cell2table(sturegr,'VariableNames', ...
        {'t','critical_t','p_value','Power','Comment'}, ...
        'RowNames',{'Slope','Intercept','Residuals'}))

    %% Power of regression (based on Fisher's Z-transform of Pearson r)
    Zrho = 0.5*reallog((1+abs(r(1)))/(1-abs(r(1))));
    sZ   = realsqrt(1/(n-3));
    pwr  = 1 - tcdf(1.96-Zrho/sZ,n-2)*2;
    disp('Power of regression')
    disp(tr)
    disp(array2table([alpha n Zrho sZ pwr], ...
        'VariableNames',{'alpha','points','Z','Sd','Power_of_regression'}))
     
    %% OTHER REGRESSIONS
    disp(' ')
    disp('OTHER REGRESSIONS')
    disp('REGRESSION SETTING Y AS INDEPENDENT VARIABLE')
    disp(tr)

    % Regression setting Y as independent variable (reverse regression)
    % slope = devy / codevxy; intercept accordingly
    regrym = devy/codevxy;
    regryq = -regrym*(mx - codevxy/devy*my);
    disp(array2table([regrym regryq],'VariableNames', ...
        {'slope','intercept'}))
    disp(' ')

    %% Regression on principal standardized component
    % Geometric mean of the two slopes, if both > 0
    cpsm = NaN;
    cpsq = NaN;
    if slope.value>0 && regrym>0
        disp('REGRESSION ON PRINCIPAL STANDARDIZED COMPONENT')
        cpsm = geomean([slope.value regrym]);
        cpsq = my - cpsm*mx;
        disp(array2table([cpsm cpsq],'VariableNames', ...
            {'slope','intercept'}))
    end

    %% Deming regression
    % Deming minimizes squared perpendicular distances to the line,
    % accounting for measurement error in both X and Y.
    z_ = (vy - (lambda*vx));
    demingm = (z_ + realsqrt(z_^2 + 4*lambda*covarxy^2)) / (2*covarxy);
    demingq = my - demingm*mx;

    disp(' ')
    disp('DEMING''S REGRESSION')
    disp(tr)
    disp(array2table([lambda demingm demingq],'VariableNames', ...
        {'Lambda','slope','intercept'}))
    
    %% Plot: main regression and residuals
    figure('Color',[1 1 1],'outerposition',get(groot,'ScreenSize'));
    
    % Left: data + LS regression + CI
    subplot(1,2,1);
    if isvector(y)        
        plot(x,yt,'bo',xtmp,ystar,'k',xtmp,cir,'r:',xtmp,cir2,'g:');
    else       
        hold on
        plot(x',yt,'LineStyle','none','Marker','o','MarkerEdgeColor','b')
        plot(xtmp,ystar,'k',xtmp,cir,'r:',xtmp,cir2,'g:');
        hold off
    end
    axis square
    txt = sprintf(['Red dotted lines: 95%% Confidence interval of regression\n', ...
                   'Green dotted lines: 95%% Confidence interval of new y evaluation using this regression']);
    title(txt)
    
    % Right: residuals with RSE bands and approximate normal cut-offs
    subplot(1,2,2);
    xl = [min(x) max(x)];
    plot(xtmp,R,'bo', ...
         xl,[0 0],'k-', ...
         xl,[RSE RSE],'g--', ...
         xl,[-RSE -RSE],'g--', ...
         xl,1.96.*[RSE RSE],'m--', ...
         xl,-1.96.*[RSE RSE],'m--', ...
         xl,2.58.*[RSE RSE],'r--', ...
         xl,-2.58.*[RSE RSE],'r--');
    axis square
    title('Residuals and reference bands')

       %% Plot: comparison of different regression lines
    % Here we illustrate:
    %  - OLS (X independent): vertical distances are minimized
    %  - "Y independent" regression (shown with vertical distances)
    %  - Principal standardized component (with vertical distances)
    %  - Deming regression: perpendicular distances are minimized
    figure('Color',[1 1 1],'outerposition',get(groot,'ScreenSize'));
    
    % Use the cleaned data used for regression (xtmp, ytmp)
    xg = linspace(min(xtmp),max(xtmp),100);
    
    % ---------- 1) OLS: X as independent (vertical residuals) ----------
    subplot(2,2,1);
    hold on
    % OLS line
    yhat_ols = slope.value.*xg + intercept.value;
    plot(xtmp,ytmp,'bo',xg,yhat_ols,'b-','LineWidth',2);
    
    % Vertical distances from each point to OLS line:
    % (xi, yi) --> (xi, m*xi+q)
    yhat_points_ols = slope.value.*xtmp + intercept.value;
    for i = 1:numel(xtmp)
        plot([xtmp(i) xtmp(i)], [ytmp(i) yhat_points_ols(i)], 'k--'); % vertical segment
    end
    hold off
    axis equal
    title({'OLS regression (X independent)', ...
           'Vertical distances minimized'})
    xlabel('X'); ylabel('Y');

    % ---------- 2) Regression with Y as independent (vertical segments) ----------
    subplot(2,2,2);
    hold on
    % Reverse regression line
    yhat_rev = regrym.*xg + regryq;
    plot(xtmp,ytmp,'ro',xg,yhat_rev,'r-','LineWidth',2);
    
    % Vertical distances to the reverse-regression line:
    % (xi, yi) --> (xi, regrym*xi + regryq)
    yhat_points_rev = regrym.*xtmp + regryq;
    for i = 1:numel(xtmp)
        plot([xtmp(i) xtmp(i)], [ytmp(i) yhat_points_rev(i)], 'k--');
    end
    hold off
    axis equal
    title({'Regression with Y as independent', ...
           'Vertical distances (for comparison)'})
    xlabel('X'); ylabel('Y');

    % ---------- 3) Principal standardized component (vertical segments) ----------
    subplot(2,2,3);
    hold on
    plot(xtmp,ytmp,'go');
    if ~isnan(cpsm)
        % Principal component line
        yhat_cps = cpsm.*xg + cpsq;
        plot(xg,yhat_cps,'g-','LineWidth',2);
        
        % Vertical distances to the PSC line:
        % (xi, yi) --> (xi, cpsm*xi + cpsq)
        yhat_points_cps = cpsm.*xtmp + cpsq;
        for i = 1:numel(xtmp)
            plot([xtmp(i) xtmp(i)], [ytmp(i) yhat_points_cps(i)], 'k--');
        end
        
        title({'Principal standardized component', ...
               'Vertical distances (for comparison)'})
    else
        title('Principal standardized component (not defined)')
    end
    hold off
    axis equal
    xlabel('X'); ylabel('Y');

    % ---------- 4) Deming regression (perpendicular distances) ----------
    subplot(2,2,4);
    hold on
    yhat_dem = demingm.*xg + demingq;
    plot(xtmp,ytmp,'mo',xg,yhat_dem,'m-','LineWidth',2);

    % Perpendicular distances from each point to Deming line:
    % For each point (x0,y0) and line y = a*x + b:
    %   x_perp = (x0 + a*(y0 - b)) / (1 + a^2)
    %   y_perp = a*x_perp + b
    a = demingm;
    b = demingq;
    for i = 1:numel(xtmp)
        x0 = xtmp(i);
        y0 = ytmp(i);
        if isfinite(a)
            x_perp = (x0 + a*(y0 - b)) / (1 + a^2);
            y_perp = a*x_perp + b;
        else
            % Degenerate case: nearly vertical line, project horizontally
            x_perp = mx;
            y_perp = y0;
        end
        plot([x0 x_perp],[y0 y_perp],'k--'); % perpendicular segment
    end
    hold off
    axis equal
    title({'Deming regression', ...
           'Perpendicular distances minimized'})
    xlabel('X'); ylabel('Y');
end

