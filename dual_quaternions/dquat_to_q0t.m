function [ q0 , t ] = dquat_to_q0t( q )
%DQUAT_TO_Q0T 

q0 = q(1:4);
t = quat_mult( q(5:8,1) , quat_conj(q(1:4,1)) );
t = [0;2*t(2:4)];

end

