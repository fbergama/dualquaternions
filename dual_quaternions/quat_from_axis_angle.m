function q = quat_from_axis_angle( axis, angle )
%QUAT_FROM_AXIS_ANGLE
%   Returns the unit quaternion representing a rotation around an arbitrary
%   axis
%
%   Syntax:
%   q = quat_from_axis_angle( axis, angle )
%
%   In:
%   axis  -  3x1 vector
%   angle -  rotation angle in radians
%
%   Out:
%   q  - Resulting unit quaternion
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


if abs(norm(axis) - 1) > 1E-8
    error('Axis must have unitary norm');
end

sinht = sin(angle/2);
q = [ cos(angle/2 ) ; sinht.*axis ];


end

