function qc = quat_conj( q )
%QUAT_CONJ 
%   Returns the conjugate of a quaternion q
%
%   Syntax:
%   qc = quat_conj( q )
%
%   In:
%   q  - Quaternion expressed as a 4x1 matrix
%
%   Out:
%   qc - q conjugate
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

qc = q.*[1;-1;-1;-1];


end

