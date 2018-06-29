function P = dquat_empirical_covariance( qts, qmean )
%DQUAT_EMPIRICAL_COVARIANCE

K = size(qts,2);

P = zeros( size(qts,1), size(qts,1) );
for ii=1:K
    qq_tm = dquat_Log(qmean,qts(:,ii));
    P = P + qq_tm*qq_tm';    
end


end

