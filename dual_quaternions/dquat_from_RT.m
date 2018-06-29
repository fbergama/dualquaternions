function q = dquat_from_RT( RT )
%DQUAT_FROM_RT Returns the dual quaternion associated to the rigid
%              motion defined by the 4x4 matrix in homogeneous coordinates
%
%   Syntax:
%   q = dquat_from_RT( RT )
%
%   In:
%   RT    - Input rigid motion in homogeneous coordinates as a 4x4 matrix
%
%   Out:
%   q     -  resulting dual quaternion
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

q=zeros(8,1);

r = quat_from_R(RT(1:3,1:3));
q(5:8) = quat_mult( [0;RT(1:3,4)], r) .* 0.5;

q(1:4)=r;

end

