function  X=rays_intersec(Rv, Re)
%M num of rays
%N 3
m = size(Rv,1);

A = zeros(m*3,3);
A(1:3:end,:) =  [zeros(m,1) -Rv(:,3)    Rv(:,2)];
A(2:3:end,:) =  [Rv(:,3)    zeros(m,1)  -Rv(:,1)];
A(3:3:end,:) =  [-Rv(:,2)   Rv(:,1)     zeros(m,1)];
% A(1:3:end,:) =  [zeros(m,1) Rv(:,3)    -Rv(:,2)   ];
% A(2:3:end,:) =  [-Rv(:,3)   zeros(m,1)  Rv(:,1)   ];
% A(3:3:end,:) =  [Rv(:,2)    -Rv(:,1)    zeros(m,1)];

Re= Re';
B = Re(:);

X = A\B;

% X = -(A'*A)\(A'*B)
% X = lsqnonneg(A,B);
end
