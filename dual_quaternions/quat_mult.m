function qm = quat_mult( q1, q2 )
%QUAT_MULT
%   Returns the product of q1 and q2
%
%   Syntax:
%   qm = quat_mult( q1, q2 )
%
%   In:
%   q1  - frist quaternion
%   q2  - second quaternion
%
%   Out:
%   q  - the product between q1 and q2 (q1*q2)
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

qm = zeros(4,1);

aw=q1(1,1);
ai=q1(2,1);
aj=q1(3,1);
ak=q1(4,1);

bw=q2(1,1);
bi=q2(2,1);
bj=q2(3,1);
bk=q2(4,1);


qm(1,1) = aw*bw - ai*bi - aj*bj - ak*bk;
qm(2,1) = ai*bw + aw*bi - ak*bj + aj*bk;
qm(3,1) = aj*bw + ak*bi + aw*bj - ai*bk;
qm(4,1) = ak*bw - aj*bi + ai*bj + aw*bk;

      
end

