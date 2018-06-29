function l = dquat_screw_length( q, p )
%DQUAT_SCREW_LENGTH Returns the length of the screw path of a point p under
%                   the transformation q
%
%
%   The path is computed as:
%   norm( log(q) p - p log(q)* )  (See Torsello's technical report)
%
%   Syntax:
%   l = dquat_screw_length( q, p )
%
%   In:
%   q  - Input dual quaternion
%   p  - Input point in R3 (3x1 Matrix)
%
%   Out:
%   l  - Screw path length
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

pq = [1;0;0;0;0;p];
logq = dquat_LogI(q);
d = dquat_mult( logq, pq ) - dquat_mult( pq, dquat_conj_dual(logq) );

l = sqrt(d' * d );


end

