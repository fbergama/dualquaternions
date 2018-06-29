function q = dquat_from_screw_motion( th0, the, s0, r )
%DQUAT_FROM_SCREW_MOTION Returns the dual quaternion associated to the
%                        screw motion defined by th0, th3, s0, se
%
%   Syntax:
%   q = dquat_from_screw_motion( th0, the, s0, se )
%
%   In:
%   th0  -  Rotation angle in radians (scalar)
%   the  -  Amount of translation along the axis (scalar)
%   s0   -  Axis direction (3x1 matrix)
%   r    -  Axis point     (3x1 matrix)
%
%   Out:
%   q     -  resulting dual quaternion
%
%
%
% Copyright 2014 Filippo Bergamasco
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

s0 = s0/norm(s0);
se = cross(r,s0);
%dot( s0, se ) %<- debug

costh0m = cos(th0*0.5);
sinth0m = sin(th0*0.5);

S0 = [0;s0];
Se = [0;se];


q = zeros(8,1);
q(1:4) = [costh0m;0;0;0] + S0.*sinth0m;
q(5:8) = [-the*0.5*sinth0m;0;0;0] + S0.*(the*0.5)*costh0m + Se*sinth0m;

q = dquat_normalize(q);

end

