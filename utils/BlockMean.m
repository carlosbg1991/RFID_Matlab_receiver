function Y = BlockMean(X, V, W)
% 2D block mean over 1st and 2nd dim [MEX]
% The mean of V*W elements along the 1st and 2nd dimensions is calculated.
% Y = BlockMean(X, V, W)
% INPUT:
%   X: UINT8 or DOUBLE array of any size.
%   V, W: Scalar numerical with integer value. Each element of the output is
%      the mean over V*W neighbouring elements of the input.
%      V and W are limited to 256 to limit memory usage.
%      A square V*V block is used, if W is omitted.
% OUTPUT:
%   Y: UINT8 or DOUBLE array, the 1st and 2nd dimensions are V and W times
%      shorter: [FLOOR(X / V) x FLOOR(Y / W) x (further dims...)]
%      If the size of the 1st or 2nd dimension is not a multiple of V and W,
%      the remaining elements at the end are skipped.
%      The empty matrix is replied for empty inputs or if the 1st or 2nd
%      dimension is shorter than V or W.
%
% NOTE: This is implemented for DOUBLE and UINT8, because I used it for
%   anti-aliasing of images stored as 3D RGB arrays.
%
% COMPILE: See BlockMean.c
%
% TESTING: Run uTest_BlockMean to check validity and speed.
%
% Tested: Matlab 6.5, 7.7, 7.8, WinXP
% Author: Jan Simon, Heidelberg, (C) 2009-2010 matlab.THISYEAR(a)nMINUSsimon.de

% $JRev: R0r V:013 Sum:DzKusULxGsEb Date:22-Sep-2010 02:30:25 $
% $License: BSD (see Docs\BSD_License.txt) $
% $UnitTest: uTest_BlockMean $
% $File: Tools\GLMath\BlockMean.m $
% History:
% 001: 20-Jul-2009 23:22, Generalized ChunkMeanRGB: Free trailing dimensions.
% 010: 12-Mar-2010 09:34, Rectangular blocks. TestBlockMean fixed.
%      The TestBlockMean published on the FEX at 21-Jul-2009 missed the function
%      isEqualTol - sorry.

% ==============================================================================
% This M-function is just a proof of concept and used for testing the results
% of the MEX version. The MEX version has a better check of the inputs.
% If M- and MEX-version are found in the path, the MEX is preferred.

% Get size of X and calculate the reduced sizes for the 1st and 2nd dimension:
if nargin < 3
   W = V;
end
S = size(X);
M = S(1) - mod(S(1), V);
N = S(2) - mod(S(2), W);
if M * N == 0
   Y = X([]);  % Copy type of X
   return;
end
MV = M / V;
NW = N / W;

% Cut and reshape input such that the 1st and 3rd dimension have the lengths V
% and W:
XM = reshape(X(1:M, 1:N, :), V, MV, W, NW, []);

% Different methods depending on the type of the input:
if isa(X, 'double')
   Y = sum(sum(XM, 1), 3) .* (1.0 / (V * W));
elseif uint8(0.8) == 1  % UINT8 of Matlab7 rounds:
   % NOT .*1/(V*W) !!! Otherwise the rounding differs from C-Mex
   Y = uint8(sum(sum(XM, 1), 3) ./ (V * W));
else                    % UINT8 of Matlab6 truncates, so round manually:
   Y = uint8(round(sum(sum(XM, 1), 3) ./ (V * W)));
end

% Remove singleton dimensions:
S(1) = MV;
S(2) = NW;
Y    = reshape(Y, S);

return;
