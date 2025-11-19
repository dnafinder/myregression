function x = myregrinv(xc,yc,yo,varargin)
% MYREGRINV Inverse regression and calibration analysis.
%
%   This function resolves a calibration (inverse regression) problem:
%   given a linear calibration line y = m*x + q obtained from calibration
%   points (xc, yc), it estimates the corresponding x values for a set of
%   observed responses yo, together with their 95% confidence intervals.
%
%   The function also computes:
%     - calibration quality index (Sokal & Rohlf)
%     - limit of detection (LOD)
%     - limit of quantification (LOQ)
%
%   It relies on MYREGR to perform the underlying linear regression.
%
%   Syntax:
%       x = myregrinv(xc,yc,yo)
%       x = myregrinv(xc,yc,yo,verbose)
%
%   Inputs:
%       xc      - Calibration x-values (independent variable).
%                 Numeric vector (row or column), strictly increasing.
%
%       yc      - Calibration y-values (dependent variable).
%                 If yc is a vector, single measurement per xc.
%                 If yc is a matrix, each row is a repeated measure
%                 of the corresponding xc point; row means are used.
%
%       yo      - Observed responses (unknown samples).
%                 If yo is a vector, single measurement per sample.
%                 If yo is a matrix, each row is a repeated measure
%                 of a given unknown sample; row means are used.
%
%       verbose - Flag to display calibration and inverse regression
%                 information (default = 1).
%
%   Outputs:
%       x       - N x 3 matrix, with:
%                   1) inverse predicted x
%                   2) lower 95% CI
%                   3) upper 95% CI
%
%   Example:
%       xc = 0:2:12;
%       yc = [1.2 5 9 12.6 17.3 21 24.7];
%       yo = [1 2 6 22 30];
%       x  = myregrinv(xc,yc,yo);
%
%   References:
%       Sokal R.R. & Rohlf F.J. (2003)
%       Biometry. The Principles and Practice of Statistics in Biological
%       Research (3rd ed., 8th printing). Freeman and Company, New York,
%       pp. 491–493.
%
%   SEE ALSO: myregr, myregrcomp
%
%   Created by Giuseppe Cardillo
%   giuseppe.cardillo.75@gmail.com
%
%   To cite this file, this would be an appropriate format:
%   Cardillo G. (2007) MyRegressionINV: resolve a calibration problem that
%   is: to estimate mean value and confidence interval of x since y.
%   https://github.com/dnafinder/myregrinv
%
%   This code is released under the GNU GPL-3.0 license.


%% Input error handling
p = inputParser;

% xc: numeric vector, increasing, row or column (like myregr uses 'vector')
addRequired(p,'xc',@(x) validateattributes(x,{'numeric'}, ...
    {'vector','real','finite','nonnan','nonempty','increasing'}));

% yc: numeric 2D array (vector or matrix, as in myregr)
addRequired(p,'yc',@(x) validateattributes(x,{'numeric'}, ...
    {'2d','real','finite','nonnan','nonempty'}));

% yo: numeric 2D array (vector or matrix of replicates)
addRequired(p,'yo',@(x) validateattributes(x,{'numeric'}, ...
    {'2d','real','finite','nonnan','nonempty'}));

% verbose flag (same style as myregr)
addOptional(p,'verbose',1, @(x) isnumeric(x) && isreal(x) && isfinite(x) && ...
    isscalar(x) && (x==0 || x==1));

parse(p,xc,yc,yo,varargin{:});
verbose = p.Results.verbose;
clear p

% Ensure xc is a column vector (like x in myregr)
xc = xc(:);

%% Check that MYREGR is available
assert(exist('myregr.m','file')~=0, ...
    'myregrinv:myregrNotFound', ...
    ['myregr.m must be on the MATLAB path.\n' ...
     'Download it from: https://github.com/dnafinder/myregr']);

%% Perform calibration regression using MYREGR
% Same calling style: [slope,intercept,STAT] = myregr(...)
[m,q,stat] = myregr(xc,yc,verbose);

% Calibration quality index (Sokal & Rohlf)
% quality = ( (t_crit * RSE / slope)^2 ) / sum of squares of X about the mean
quality = ((stat.cv * stat.rse / m.value)^2) / stat.sse;

if verbose==1
    disp(' ')
    if quality >= 0.1
        fprintf('quality = %0.4f >= 0.1\t This is not a good calibrator\n',quality)
    else
        fprintf('quality = %0.4f < 0.1\t This is a good calibrator\n',quality)
    end
end

%% Limits of detection and quantification
% LOD and LOQ defined on the Y scale and then back-calculated to X
lod = q.value + 3*q.se;          % limit of detection (on Y)
loq = lod + 7*q.se;              % limit of quantification (on Y)

%% Prepare yo (unknown responses)
if isvector(yo)
    yo = yo(:);                  % column vector
else
    % If yo is a matrix, each row is a repeated measure
    yo = mean(yo,2);             % mean over rows
end

%% Inverse prediction (point estimates)
% Direct inverse of calibration line: x_hat = (y - q) / m
xo = (yo - q.value) ./ m.value;

%% Confidence interval of inverse prediction
% Following Zar's formulation (Biostatistical Analysis), the CI for x given y
% uses:
%   - distance of y from the calibration mean (a)
%   - regression sum of squares (stat.sse)
%   - residual standard error (stat.rse)
%   - t critical value (stat.cv)
%
% K is a function of slope and its uncertainty:
K = m.value^2 - stat.cv^2 * m.se^2;

% a: scaled squared distance of each yo from mean response
a = ((yo - stat.ym).^2) ./ stat.sse;

% b: term including sample size and K
b = K * (1 + 1/stat.n);

% c: RSE * sqrt(a + b)
c = stat.rse * realsqrt(a + b);

% d: half-width of the CI in x-space, scaled by t and K
d = (stat.cv ./ K) .* c;

% e: center of the CI in x-space
e = stat.xm + (m.value .* (yo - stat.ym) ./ K);

% Final [lower, upper] confidence bounds for each yo
f = repmat(e,1,2) + repmat([-1 1],length(yo),1) .* repmat(d,1,2);

% Assemble output matrix: [x_hat, x_lower, x_upper]
x = [xo f];

%% Display LOD / LOQ and classification table
if verbose==1
    disp(' ')
    fprintf('Limit of detection (LOD): %0.4f\t x_LOD = %0.4f\n', ...
        lod, (lod - q.value) ./ m.value)
    fprintf('Limit of quantification (LOQ): %0.4f\t x_LOQ = %0.4f\n', ...
        loq, (loq - q.value) ./ m.value)
    disp(' ')
    
    tr = repmat('-',1,60);
    disp('                  Inverse Prediction Values'); 
    disp(tr)
end

% Calibration interval on the Y scale (for "within calibration" checks)
if isvector(yc)
    calint = [min(yc) max(yc)];
else
    calint = [min(mean(yc,2)) max(mean(yc,2))];
end

%% Classify each yo and adjust output x when needed
for I = 1:length(yo)
    if yo(I) > loq
        % Above LOQ
        if yo(I) >= calint(1) && yo(I) <= calint(2)
            % Within calibration interval: full inverse prediction
            if verbose==1
                fprintf('%10.4f     %10.4f     %10.4f   %10.4f\n', ...
                    yo(I), x(I,:));
            end
        else
            % Out of calibration interval
            if verbose==1
                fprintf('%10.4f     Out of calibration interval\n', yo(I));
            end
            x(I,:) = NaN(1,3);
        end
    elseif yo(I) > lod && yo(I) <= loq
        % Between LOD and LOQ (quantification not reliable)
        if verbose==1
            fprintf('%10.4f     LOD<=x<=LOQ\n', yo(I));
        end
        x(I,:) = NaN(1,3);
    elseif yo(I) <= lod
        % Below LOD (signal indistinguishable from blank)
        if verbose==1
            fprintf('%10.4f          x<LOD\n', yo(I));
        end
        x(I,:) = NaN(1,3);
    end
end

if verbose==1
    disp(tr)
end
