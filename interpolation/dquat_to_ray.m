function [ Rp Rv ] = dquat_to_ray( dq )
%DQUAT_TO_RAY 

RT = dquat_to_RT(dq);
Rv = RT(1:3,1);
Rp = RT(1:3,4);

end

