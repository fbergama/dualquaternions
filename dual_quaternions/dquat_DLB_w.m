function qm = dquat_DLB_w( dquats, weights )
%DQUAT_DLB_W 
% Computes the Dual-quaternions Linear Blending
%   
%   Reference:
%   Ladislav Kavan, Steven Collins, Carol O'Sullivan, Jiri Zara,
%   "Dual Quaternions for Rigid Transformation Blending", 
%   Tech. rep. TCD-CS-2006-46, Trinity College Dublin
% 
%
%   Syntax:
%   qm = dquat_DLB_w( dquats, weights )
%
%   In:
%   dquats  -  input dual quaternions as an 8xN matrix
%   weights -  (optional) input weights as 1xN vector (if omitted, uniform
%                         weights distribution is assumed).
%
%   Out:
%   qm  - weighted average of input dual quaternions
%
%
%
% Copyright 2013 Filippo Bergamasco
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
% http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%

if nargin<2
    weights = ones( 1, size(dquats,2)  )/size(dquats,2);   % assume uniform weights
end

if size(weights,2)<size(weights,1)
    weights=weights';
end

if size(weights,2) ~= size(dquats,2)
    error('Size of weights vector differs from number of given dual quaternions.');
end

if abs(sum(weights)-1)>1E-9
    error('sum(weights) ~= 1');
end

qm = dquat_normalize( sum(dquats .* repmat( weights, 8, 1),2) );


end

