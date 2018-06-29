function [ th_h, s ] = dquat_to_screw_form( dq )
%DQUAT_TO_SCREW_FORM  Returns the screw form of the input unit dual
%                    quaternion dq
%
%   Any unit dual quaternion can be represented as:
%       cos( th_h ) + i s sin( th_h )
%
%   where:
%       i = [i;j;k]             is the row vector of quaternions imaginary 
%                               bases
%       th_h = th0 + eps the    is a dual number describing the quantity of 
%                               motion,  specifically:
%                                  th0 is the angle in radians
%                                  the is the amount of translation
%
%       i s = i s0 + eps i se   is a dual quaternion (with zero scalar
%                               part) that describes the axis:
%                                  i s0 is the axis direction
%                                  i se is the position of the axis such
%                                       that se=cross(d,s0) for any point d
%                                       on the axis
%                                 
%
%
%   Syntax:
%   [ th_h, s ] = dual_to_screw_form( dq )
%
%   In:
%   dq    - input unit dual quaternion (8x1 matrix)
%
%   Out:
%   th_h  -  dual number describing the quantity of motion
%   s     -  dual number describing the axis
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

[q0, t]=dquat_to_q0t(dq);
[axis, angle]=quat_to_axis_angle( q0 );

theta0 = angle;
s0 = [0;axis];

if sum(abs(q0-[1;0;0;0])) < 1E-8
    theta0=0;
    thetae=quat_norm(t);
    s0=t/quat_norm(t);
    se = [0;0;0;0];
elseif sum(abs(q0-[-1;0;0;0])) < 1E-8
    theta0=2*pi;
    thetae=quat_norm(t);
    s0=t/quat_norm(t);
    se = [0;0;0;0];
else
    thetae = quat_dot(t,s0);
    
    se = 0.5*cross( cross(s0(2:4),t(2:4))*cot(theta0/2)+t(2:4), s0(2:4));
    se = [0;se];    
end


th_h = dual_number(theta0,thetae);
s = [s0;se];

end

