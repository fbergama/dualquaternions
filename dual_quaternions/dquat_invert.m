function qs = dquat_invert( q )
%DQUAT_GEODESIC
% Returns the unit dual quaternion that represents the inverse of the rigid
% motion defined by the input unit dual quaternion q
%
% Note: this is the same operation than calling dquat_conj_quat(q) but it's
%       also checked if q has unitary norm
%
%   Syntax:
%   qs = dquat_invert( q )
%
%   In:
%   q  - input unitary dual quaternion as 8x1 matrix
%
%   Out:
%   qs  - the inverse of q
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

n = dquat_norm(q);
if abs(n.a0-1)>1E-8 || abs(n.ae)>1E-8
    error('Unable to invert a non unit dual quaternion');
end


qs = dquat_conj_quat(q);

end

