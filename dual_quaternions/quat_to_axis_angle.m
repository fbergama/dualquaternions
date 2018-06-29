function [ axis , angle ] = quat_to_axis_angle( q )
%QUAT_TO_AXIS_ANGLE
%   Recovers the rotation axis and angle from an unitary quaternion
%
%   Syntax:
%   [ axis , angle ] = quat_to_axis_angle( q )
%
%   In:
%   q  - Input quaternion expressed as a 4x1 matrix
%
%   Out:
%   axis  - rotation axis as a 3x1 matrix
%   angle - rotation angle in radians 
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

angle = 2*atan2( norm(q(2:4)),q(1,1) );

if angle < 1E-8
    angle = 0;
    axis = [0;0;0];
else
    axis = q(2:4)/sin(angle/2);
end

end

