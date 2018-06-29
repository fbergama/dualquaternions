function q = quat_from_R( R )
%QUAT_FROM_R
%   Returns the unit quaternion representing a rotation matrix R
%
%   Syntax:
%   q = quat_from_R( R )
%
%   In:
%   R  - 3x3 rotation matrix
%
%   Out:
%   q  - Unit quaternion that represents the rotation R
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

q=zeros(4,1);

m00 = R(1,1);
m11 = R(2,2);
m22 = R(3,3);
m01 = R(1,2);
m02 = R(1,3);
m10 = R(2,1);
m12 = R(2,3);
m20 = R(3,1);
m21 = R(3,2);
tr = m00 + m11 + m22;

if tr > 0
    
    s = sqrt(tr+1) * 2;
    q(1,1) = 0.25 * s;
    q(2,1) = (m21 - m12) / s;
    q(3,1) = (m02 - m20) / s;
    q(4,1) = (m10 - m01) / s;

elseif ((m00 > m11)&&(m00 > m22))

    s = sqrt(1 + m00 - m11 - m22) * 2; 
    q(1,1) = (m21 - m12) / s;
    q(2,1) = 0.25 * s;
    q(3,1) = (m01 + m10) / s;
    q(4,1) = (m02 + m20) / s;
 
elseif (m11 > m22)

    s = sqrt(1 + m11 - m00 - m22) * 2;
    q(1,1) = (m02 - m20) / s;
    q(2,1) = (m01 + m10) / s;
    q(3,1) = 0.25 * s;
    q(4,1) = (m12 + m21) / s;

else

    s = sqrt(1 + m22 - m00 - m11) * 2;
    q(1,1) = (m10 - m01) / s;
    q(2,1) = (m02 + m20) / s;
    q(3,1) = (m12 + m21) / s;
    q(4,1) = 0.25 * s;
 
end



end

