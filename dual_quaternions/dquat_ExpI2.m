function  q  = dquat_ExpI2( qt )
%DQUAT_EXPI2

error('do not use this function')
k = inv( dquat_norm(qt) );

s = dquat_mult_dual_number(qt,k);
th = dquat_norm(qt);
th.a0=th.a0*2;
th.ae=th.ae*2;

cc = cos(th);
ss = sin(th);

sss = dquat_mult_dual_number(s,ss);


q = zeros(8,1);
q(1:4) = [cc.a0;0;0;0]+sss(1:4);
q(5:8) = [cc.ae;0;0;0]+sss(5:8);

end

