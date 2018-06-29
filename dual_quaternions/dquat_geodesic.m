function x = dquat_geodesic( p, q, t )
%DQUAT_GEODESIC
% Returns the point k lying on the geodesic path between p and q at ratio t
% For example:
%    dquat_geodesic( p,q,0 )=p
%    dquat_geodesic( p,q,1 )=q
%    dquat_geodesic( p,q,0.5 )=k (exactly in between p and q)
%
%   Syntax:
%   x = dquat_geodesic( p, q, t )
%
%   In:
%   p  - path starting point
%   q  - path end point
%   t  - path length
%
%   Out:
%   x  - output manifold point
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

x = dquat_Exp( p, dquat_Log(p,q).*t );

end

