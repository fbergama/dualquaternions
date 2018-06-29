function n = quat_dot( q1, q2 )
%QUAT_DOT 
%   Returns the dot product between q1 and q2
%
%   Syntax:
%   n = quat_norm( q )
%
%   In:
%   q1  - First quaternion expressed as a 4x1 matrix
%   q2  - Second quaternion expressed as a 4x1 matrix
%
%   Out:
%   n  - <q1,q2>
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

n = q1'*q2;

end

