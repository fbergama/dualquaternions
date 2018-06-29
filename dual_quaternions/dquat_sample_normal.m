function q = dquat_sample_normal( p, P )
%DQUAT_SAMPLE_NORMAL 
%   Samples a dual quaternion from a normal distribution with mean p
%   and covariance P (expressed in dual quaternion identity's tangent
%   space)
%

qt = mvnrnd(zeros(8,1),P)';         % sample a vector in 1 tangent space
q = dquat_mult(p, dquat_ExpI(qt) ); % move to manifold and transform 

end

