function R = quat_to_R( q )
%QUAT_TO_R 
%   Returns the rotation matrix represented by the unit quaternion q
%
%   Syntax:
%   R = quat_to_R( q )
%
%   In:
%   q  - Input quaternion as a 4x1 matrix
%
%   Out:
%   R  - 3x3 rotation matrix represented by q
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

w = q(1,1);
i = q(2,1);
j = q(3,1);
k = q(4,1);

R(1,1) = 1 - 2*j*j - 2*k*k;
R(1,2) = 2*i*j - 2*k*w;
R(1,3) = 2*i*k + 2*j*w;

R(2,1) = 2*i*j + 2*k*w;
R(2,2) = 1 - 2*i*i - 2*k*k;
R(2,3) = 2*j*k - 2*i*w;

R(3,1) = 2*i*k - 2*j*w;
R(3,2) = 2*j*k + 2*i*w;
R(3,3) = 1 - 2*i*i - 2*j*j;


end

