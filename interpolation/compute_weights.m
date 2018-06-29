function w = compute_weights( C, cs, lambda )
%COMPUTE_WEIGHTS 

% Add conditions to ensure that sum(w)=1
C = [C; ones(1,size(C,2))];
cs = [cs;1];

I = eye(size(C,2));
D = eye(size(C,2));

KK = C-repmat(cs,1,size(C,2));
D = (KK'*KK).*eye(size(C,2));

Z = -0.5*inv(D+lambda*I)*C';
A = C*Z;
lambdas =linsolve(A,cs);

w = Z*lambdas;

end

