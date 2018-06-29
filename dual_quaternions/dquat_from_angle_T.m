function q = dquat_from_angle_T( alpha, beta, gamma, t )
%DQUAT_FROM_ANGLE_T
%   Creates the unit dual quaternion representing the rigid motion with
%   translation t and rotation with euler angles: alpha, beta and gamma
%
%   Syntax:
%   q = dquat_from_angle_T( alpha, beta, gamma, t )
%
%   In:
%   alpha  -  rotation angle in radians with respect to x-axis
%   beta   -  rotation angle in radians with respect to y-axis
%   gamma  -  rotation angle in radians with respect to z-axis
%   t      -  translation vector as a 3x1 matrix
%
%   Out:
%   q  - resulting unit dual quaternion
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

Rx = [  1       0           0          0;...
        0   cos(alpha)  -sin(alpha)    0;...
        0   sin(alpha)   cos(alpha)    0;...
        0       0           0          1];
        
       
Ry = [  cos(beta)       0           sin(beta)     0;...
           0            1               0         0;...
        -sin(beta)      0           cos(beta)     0;...
           0            0               0         1];
       
        
Rz = [  cos(gamma)    -sin(gamma)         0         0;...
        sin(gamma)     cos(gamma)         0         0;...
           0            0                 1         0;...
           0            0                 0         1];
       
Tm=eye(4,4);
Tm(1:3,4)=t;
q=dquat_from_RT( Tm*Rz*Ry*Rx );

end

