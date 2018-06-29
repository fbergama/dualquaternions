function [qmean, P] = dquat_empirical_covarianceI( qts )
%DQUAT_EMPIRICAL_COVARIANCE

K = size(qts,2);

qmean = dquat_DIB(qts);

% move qts back to 1 to compute covariances
qmean_c = dquat_invert( qmean );
qts1 = qts;
for ii=1:size(qts,2)
    qts1(:,ii) = dquat_mult( qmean_c, qts(:,ii) ); 
end

P = zeros( size(qts1,1), size(qts1,1) );
for ii=1:K
    qq_tm = dquat_LogI(qts1(:,ii));
    P = P + qq_tm*qq_tm';    
end


end

