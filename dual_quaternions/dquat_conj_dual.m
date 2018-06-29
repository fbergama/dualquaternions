function qc = dquat_conj_dual( q )
%DQUAT_CONJ_DUAL 
% Returns the dual number conjugate of the dual quaternion q
% 
%   if  q=q0 + ?qe,
%   dquat_conj_dual( q ) = q0 - ?qe
%
%   Syntax:
%   qc = dquat_conj_dual( q )
%
%   In:
%   q   - dual quaternion to conjugate
%
%   Out:
%   qc  - output conjugated dual quaternion
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

qc = q;
qc(5:8) = qc(5:8)*-1;

end

