function myregrcomp(x1,y1,x2,y2)
% MYREGRCOMP Compare two linear regressions.
%
%   This function compares two least-squares linear regressions fitted to
%   two independent datasets (x1,y1) and (x2,y2). It performs:
%     - a global test for equality of the two regressions
%     - a test on slopes (parallelism)
%     - a test on intercepts (coincidence)
%   according to the approach reported by Stanton A. Glantz.
%
%   The function relies on MYREGR to perform the underlying linear
%   regressions.
%
%   Syntax:
%       myregrcomp(x1,y1,x2,y2)
%
%   Inputs:
%       x1  - Independent variable of the first regression (numeric vector).
%       y1  - Dependent variable of the first regression (numeric vector).
%       x2  - Independent variable of the second regression (numeric vector).
%       y2  - Dependent variable of the second regression (numeric vector).
%
%     Note:
%       For repeated-measure designs, average replicates for each x-level
%       before calling MYREGRCOMP, so that y1 and y2 are vectors.
%
%   Outputs:
%       This function does not return variables, but:
%         - prints a summary table of the three regressions
%           (regression 1, regression 2, combined regression)
%         - prints:
%             * the global F-test
%             * the t-test on slopes
%             * the t-test on intercepts
%             * the intersection point (if slopes differ)
%         - produces a figure with:
%             * raw data of both regressions
%             * fitted lines for the two regressions
%             * optional marker at the point of intersection
%
%   Example:
%       x1 = [1 2 3 4 5 6 7 8 9 10];
%       y1 = [-0.5052 0.2045 0.9586 1.3357 1.1463 2.8586 4.0651 4.2444 5.2673 5.8634];
%       x2 = [3 5 7 9 11 13 15];
%       y2 = [-3.6517 -6.0197 -9.2270 -11.0558 -14.7402 -17.3678 -19.9047];
%       myregrcomp(x1,y1,x2,y2)
%
%   SEE ALSO: myregr, myregrinv
%
%   Created by Giuseppe Cardillo
%   giuseppe.cardillo.75@gmail.com
%   GitHub repository: https://github.com/dnafinder/myregression
%
%   To cite this file, this would be an appropriate format:
%   Cardillo G. (2007) MyRegressionCOMP: a simple routine to compare two
%   LS regressions.
%   GitHub repository: https://github.com/dnafinder/myregression
%
%   This code is released under the GNU GPL-3.0 license.

%% Input error handling
p = inputParser;

addRequired(p,'x1',@(x) validateattributes(x,{'numeric'}, ...
    {'vector','real','finite','nonnan','nonempty'}));
addRequired(p,'y1',@(x) validateattributes(x,{'numeric'}, ...
    {'vector','real','finite','nonnan','nonempty'}));

addRequired(p,'x2',@(x) validateattributes(x,{'numeric'}, ...
    {'vector','real','finite','nonnan','nonempty'}));
addRequired(p,'y2',@(x) validateattributes(x,{'numeric'}, ...
    {'vector','real','finite','nonnan','nonempty'}));

parse(p,x1,y1,x2,y2);
clear p

% Ensure column vectors for internal consistency
x1 = x1(:).';
y1 = y1(:).';
x2 = x2(:).';
y2 = y2(:).';

% Check lengths consistency
assert(numel(x1)==numel(y1), ...
    'myregrcomp:SizeMismatch', ...
    'x1 and y1 must have the same number of elements.');
assert(numel(x2)==numel(y2), ...
    'myregrcomp:SizeMismatch', ...
    'x2 and y2 must have the same number of elements.');

%% Check that MYREGR is available
assert(exist('myregr.m','file')~=0, ...
    'myregrcomp:myregrNotFound', ...
    ['myregr.m must be on the MATLAB path.\n' ...
     'Download it from: https://github.com/dnafinder/myregr']);

%% First regression
disp('First regression parameters estimation')
[m1,q1,stat1] = myregr(x1,y1,0); % silent MYREGR, myregrcomp will print

%% Second regression
disp('Second regression parameters estimation')
[m2,q2,stat2] = myregr(x2,y2,0);

% Collect basic quantities
n   = [stat1.n stat2.n];
m   = [m1.value m2.value];
mse = [m1.se    m2.se];
q   = [q1.value q2.value];
qse = [q1.se    q2.se];
rse = [stat1.rse stat2.rse];

%% Combined regression (on pooled data)
disp('Combined regression parameters estimation')
x_all = sortrows([x1' y1'; x2' y2'],1);
[m3,q3,stat3] = myregr(x_all(:,1).',x_all(:,2).',0);

%% Assemble summary table
matrix = [n          stat3.n;   ... % points
          m          m3.value; ... % slopes
          mse        m3.se;    ... % slope standard errors
          q          q3.value; ... % intercepts
          qse        q3.se;    ... % intercept standard errors
          rse        stat3.rse];   % regression standard errors

disp(array2table(matrix, ...
    'VariableNames',{'Regr1','Regr2','Regr_Tot'}, ...
    'RowNames',{'Points','Slope','Slope_se','Intercept','intercept_se','Regre_se'}))

%% Global test
disp('GLOBAL TEST')

vd = (sum(n) - 4); % denominator degrees of freedom

% Combined variance of separated regressions
cs = (sum((n-2).*rse.^2)) / vd;

% Variation in variance between combined and separated regressions
vs = ((vd+2)*stat3.rse^2 - vd*cs) / 2;

% F-test
F = abs(vs/cs);
p = 1 - fcdf(F,2,vd);

disp(array2table([F 2 vd p], ...
    'VariableNames',{'F','DF_numerator','DF_denominator','p_value'}));

if p >= 0.05
    disp('These regressions are equal')
else
    disp('These regressions are different')
    disp(' ')
    
    %% Test on slopes
    disp('TEST ON SLOPES')
    cs  = sum((n-2).*mse.^2) / vd;
    cse = realsqrt(cs .* sum(1./((n-1).*mse.^2)));
    t   = abs(diff(m)) / cse;
    p1  = 1 - tcdf(t,vd);
    
    disp(array2table([t vd p1], ...
        'VariableNames',{'t','DF','p_value'}))
    
    if p1 >= 0.05
        disp('The slopes are equal and these regressions should be parallel')
    else
        disp('The slopes are different and these regressions should not be parallel')
        disp('Evaluate if the product of slopes = -1 to check perpendicularity')
    end
    disp(' ')
    
    %% Test on intercepts
    disp('TEST ON INTERCEPTS')
    cs  = sum((n-2).*qse.^2) / vd;
    cse = realsqrt(cs .* sum(1./((n-1).*qse.^2)));
    t   = abs(diff(q)) / cse;
    p2  = 1 - tcdf(t,vd);
    
    disp(array2table([t vd p2], ...
        'VariableNames',{'t','DF','p_value'}))
    
    if p2 >= 0.05
        disp('The intercepts are equal')
    else
        disp('The intercepts are different')
    end
    
    %% Intersection point (only if slopes differ)
    if p1 < 0.05
        xc = (q2.value - q1.value) / (m1.value - m2.value);
        yc = m1.value * xc + q1.value;
        disp(' ')
        disp('These regressions cross at the point:')
        disp(table(xc,yc))
    end
end

%% Plot regressions and (optionally) intersection
figure('Color',[1 1 1],'outerposition',get(groot,'ScreenSize'));
hold on
plot(x1,y1,'ro','DisplayName','Data 1')
plot(x2,y2,'bo','DisplayName','Data 2')

if exist('xc','var')
    plot(xc,yc,'k+','MarkerSize',12,'DisplayName','Intersection')
    xg = [min([x1 x2 xc]) max([x1 x2 xc])];
else
    xg = [min([x1 x2]) max([x1 x2])];
end

plot(xg,polyval([m1.value q1.value],xg),'r-','LineWidth',2,'DisplayName','Regression 1')
plot(xg,polyval([m2.value q2.value],xg),'b-','LineWidth',2,'DisplayName','Regression 2')
xlabel('X')
ylabel('Y')
legend('Location','best')
hold off

