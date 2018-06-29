function  qm  = dquat_mult( q1, q2 )
%DQUAT_MULT
% Multiplies q1 and q2
%
%   In case of unit dual quaternions, the resulting dual quaternion
%   represents the rigid transformation that is obtained by applying
%   q2 and then q1.
% 
%
%   Syntax:
%   qm  = dquat_mult( q1, q2 )
%
%   In:
%   q1      -  First operand
%   q2      -  Second operand
%
%   Out:
%   qm  - q1*q2
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

qm=zeros(8,1);
qm(1:4)=quat_mult(q1(1:4),q2(1:4));
qm(5:8)=quat_mult(q1(1:4),q2(5:8)) + quat_mult(q1(5:8),q2(1:4));

end

