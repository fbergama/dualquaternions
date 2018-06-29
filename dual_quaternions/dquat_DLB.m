function m = dquat_DLB( p, q, t )
%DQUAT_DLB 
% Computes the Dual-quaternions Linear Blending between p and q
%   
%   Reference:
%   Ladislav Kavan, Steven Collins, Carol O'Sullivan, Jiri Zara,
%   "Dual Quaternions for Rigid Transformation Blending", 
%   Tech. rep. TCD-CS-2006-46, Trinity College Dublin
% 
%
%   Syntax:
%   m = dquat_DLB( p, q, t )
%
%   In:
%   p  -  first dual quaternion (8x1 matrix)
%   q  -  second dual quaternion (8x1 matrix)
%   t  -  blend factor (0..1)
%
%   Out:
%   m  - interpolated dual quaternion
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

m = dquat_DLB_w(  [p,q], [1-t,t] );

end

