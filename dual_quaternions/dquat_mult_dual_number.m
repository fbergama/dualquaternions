function dqout = dquat_mult_dual_number( dq, dn )
%DQUAT_MULT_DUAL_NUMBER 

dq0 = dq(1:4);
dqe = dq(5:8);

dqout=zeros(8,1);
dqout(1:4) = dq0*dn.a0;
dqout(5:8) = dq0*dn.ae + dqe*dn.a0;

end

