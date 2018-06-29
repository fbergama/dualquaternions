function qm = dquat_DIB( dquats, weights, max_iterations )
%DQUAT_DIB 
% Computes the Dual-quaternions Iterative Blending
%   
%   Reference:
%   Ladislav Kavan, Steven Collins, Carol O'Sullivan, Jiri Zara,
%   "Dual Quaternions for Rigid Transformation Blending", 
%   Tech. rep. TCD-CS-2006-46, Trinity College Dublin
% 
%
%   Syntax:
%   qm = dquat_DIB( dquats, weights, max_iterations )
%
%   In:
%   dquats  -  input dual quaternions as an 8xN matrix
%   weights -  (optional) input weights as 1xN vector (if omitted, uniform
%                         weights distribution is assumed).
%   max_iterations - (optional) maximum number of iterations
%
%   Out:
%   qm  - weighted average of input dual quaternions
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

if nargin<2
    weights = ones( 1, size(dquats,2) ) / size(dquats,2);   % assume uniform weights
end

if nargin<3
    max_iterations = 20;
end

if size(weights,2)<size(weights,1)
    weights=weights';
end

if size(weights,2) ~= size(dquats,2)
    error('Size of weights vector differs from number of given dual quaternions.');
end

nquats = size(dquats,2);
b = dquat_DLB_w(dquats,weights);

%b = dquats(:,1);
%fprintf('fix DIB\n')
%%
%test
%{
figure;
axis equal;
grid on;
for ii=1:size(dquats,2)
    render_frame(dquats(:,ii),0.1);
end
render_frame(b,0.8);
%}
%%

for ii=1:max_iterations
    x=zeros(8,1);
    
    for jj=1:size(dquats,2)
        x = x + weights(jj).*dquat_LogI( dquat_mult( dquat_conj_quat(b), dquats(:,jj) ) );
    end
    
    biter = dquat_mult( b, dquat_ExpI(x) );    
    
    % Flip check and recover:
    % we must ensure that the dot between biter and all other dual 
    % quaternions given is always positive
    dots = dquats'*biter;
    posneg = sign(dots);
    if sum(posneg) == -nquats
        % all dots are negative, so we just have to flip biter
        biter = biter *-1;
        %fprintf('!!!!!! Flipping!\n');
    end
    if sum(posneg) ~= nquats
        % some positives, some negatives
        error('Bad mean.');
    end
    
    % all dots are positive, we can carry on
    b = biter;
    
    xl = dquat_norm(x);
    if (xl.a0^2 + xl.ae^2) < 1E-20
        %fprintf('Exiting after %d/%d iterations\n',ii,max_iterations);
        break
    end
end

qm = b;


end

