function Po = dquat_parallel_transport( p,q,P,N )
%DQUAT_PARALLEL_TRANSPORT
% Transports the bilinear mapping P on TpM along
% the geodesic curve between p and q.
%
%   Syntax:
%   Po = dquat_parallel_transport( p,q,P,N )
%
%   In:
%   p  - Path starting point
%   q  - Path end point
%   P  - Bilinear mapping (must be a symmetric matrix)
%   N  - Number of steps (optional, default 20)
%
%   Out:
%   Po - transported mapping
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

% Since P (MxM) is symmetric, it can be written as
% P = sum_{m=1}^{M} L_m V_m V_m'
% where V_1 ... V_m are the eigenvectors and L_1 .. L_m 
% are the eigenvalues of P
%
% Hence, we can obtain an orthonormal base of TpM by
% computing the eigenvectors of P, parallel transport each
% eigenvector and recompose the transported P

[V,L]=eig(P,'nobalance');   
%[V,L]=eig(P);   

Po=zeros(size(P));
for ii=1:size(P,2)
    vt = dquat_parallel_transport_v(p,q,V(:,ii),N);        

    %vt = vt/norm(vt);
    %vt = V(:,ii);
    Po = Po+L(ii,ii)*(vt*vt');
end


end

