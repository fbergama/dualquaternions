function P = dquat_make_covariance( q, rx, ry, rz, tx, ty, tz )
%DQUAT_MAKE_COVARIANCE

if q ~= dquat_eye()
    error('do not use this, parallel transport is not ready yet')
end


% Generate covariance at identity's tangent space
Pi = diag( [0 rx ry rz 0 tx ty tz ] );

% Parallel transport P to q
P = dquat_parallel_transport(dquat_eye(),q,Pi);

%P=Pi;
end

