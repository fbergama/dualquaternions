function v = dquat_parallel_transport_v( p,q,v0,N )
%DQUAT_PARALLEL_TRANSPORT_V
% Transports the vector v from TpM to TqM along
% the geodesic path between p and q using Schild's Ladder.
%
%   WARNING: this needs more testing!!!
%
%   Syntax:
%   v = dquat_parallel_transport_v( p,q,v0,N )
%
%   In:
%   p  - Path starting point
%   q  - Path end point
%   v0 - Vector to transport
%   N  - Number of steps (optional, default 20)
%
%   Out:
%   v  - transported vector
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

% Arguments check
if nargin<3
    error('Too few arguments');
end
if nargin<4
    N=20;
end

% do the parallel transport
%v = dquat_mult(q,v0);
%v = dquat_mult( v0, dquat_mult( dquat_invert(p), q ) );
%v = dquat_mult(dquat_mult( q, dquat_invert(p) ),v0 );


%v = dquat_mult(dquat_mult( q, dquat_invert(p) ),v0 );
%v = dquat_mult(v0,dquat_mult( q,dquat_invert(p) ) );

%qps= dquat_mult( q, dquat_invert(p)) ;
% pqs=dquat_mult( p, dquat_invert(q)) ;
% v = dquat_mult( dquat_mult( qps, v0), pqs );



qps= dquat_mult( q, dquat_invert(p)) ;
d0 = qps(1:4);
de = qps(5:8);
vr=v0(1:4);
ve=v0(5:8);

Tvr = quat_mult( quat_mult(d0,vr), quat_conj(d0));
Tve = quat_mult( quat_mult(d0,ve), quat_conj(d0));

v = [Tvr;Tve];
%}


% X=p;
% 
% vlen = norm(v0);
% scale=1;
% v=v0/vlen*scale; % scale v to reduce numerical errors
% 
% a = dquat_Exp(p,v);
% 
% for ii=1:N;
%     
%     Xp = X;
%     X = dquat_geodesic(p,q,ii/N);      
%     b = dquat_geodesic(a,X,0.5);
%     
%     %vd = sph2_Log(Xp,b);
%     %a = sph2_geodesic_v(Xp,vd,2.0);    
%     a = dquat_geodesic(Xp,b,2);
%     
%     v = dquat_Log(X,a);
%     
% end
% 
% v = v*vlen/scale;


end

