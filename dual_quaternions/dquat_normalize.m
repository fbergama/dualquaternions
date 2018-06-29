function qn = dquat_normalize( q )
%DQUAT_NORMALIZE 
%   Normalizes the input dual quaternion to unitary norm
%
%   Syntax:
%   qn = dquat_normalize( q )
%
%   In:
%   q  - Input dual quaternion expressed as a 8x1 matrix
%
%   Out:
%   qn - Output normalized dual quaternion
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



n0 = quat_norm(  q(1:4) );
qn = zeros(8,1);

qn(1:4) = q(1:4)/n0;
qn(5:8) = -q(1:4)*quat_dot( q(1:4), q(5:8) )/(n0^3) + q(5:8)/n0;



%{
% this is also ok (found online)
nondual = q(1:4);
dual = q(5:8);

lenNondual = sqrt( quat_dot( nondual, nondual ) );
if abs(lenNondual)>1E-10
    lenRecip = 1.0 / lenNondual;
    nondual = nondual * lenRecip;
    dual = dual * lenRecip;

    r = nondual;
    r = r * quat_dot( nondual, dual );
    r = r * -1;
    dual = dual + r;

    qn = [nondual;dual];
else
    error('Unable to normalize, non dual part is almost zero');
end
%}

end

