function qs = dquat_shortest_screw_path( q )
%DQUAT_SHORTEST_SCREW_PATH Returns q or -q to ensure to minimize the screw
%                          path of any point p through q.
%
%   As of quaternions, dual quaternion ring is a dual cover of SE(3). As
%   such, both q and -q represent the same transformation but with
%   different screw paths lengths.
%   This function chooses between q and -q to keep the one that produces 
%   the shortest screw path when applyed to any point p
%
%   Syntax:
%   qs = dquat_shortest_screw_path( q )
%
%   In:
%   q  - Input dual quaternion (8x1 matrix)
%
%   Out:
%   qs  - Resulting shortest path quaternion
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

d1=0;
d2=0;

while abs(d1-d2)<1E-5
    
    % we choose a random point and test
    p = rand(3,1);
    d1 = dquat_screw_length(q, p);
    d2 = dquat_screw_length(-q, p);    
    
end

if d1<d2
    qs = q;
else
    qs = -q;
end

end

