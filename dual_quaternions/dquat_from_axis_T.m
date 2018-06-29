function q = dquat_from_axis_T( axis, angle, T )
%DQUAT_FROM_AXIS_T Returns the dual quaternion associated to the rigid
%                  motion defined by a rotation about an arbitrary axis and
%                  a translation T.
%
%   Syntax:
%   q = dquat_from_axis_T( axis, angle, T )
%
%   In:
%   axis  -  Rotation axis as a 3x1 matrix
%   angle -  Rotation angle in radians
%   T     -  Translation as a 3x1 matrix
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

r = quat_from_axis_angle( axis, angle );

tq = [0;T];

q = zeros(8,1);
q(1:4) = r;
q(5:8) = 0.5*quat_mult( tq, r );


end

