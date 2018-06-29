function q = dquat_ExpI( qt )
%DQUAT_EXPI Summary of this function goes here
%   Detailed explanation goes here

q = qt;
h0=2.0*quat_norm(q(1:4));


if (h0*h0 < 1.0E-12)
    
    %small angle approximation: sin(h0)=h0, cos(h0)=1
    q(1:4) = q(1:4)*2;
    q(1,1) = 1.0;
    q(5:8) = q(5:8)*2;
    
else
    
    he=4.0*quat_dot(qt(5:8),qt(1:4))/h0;
    sh0=sin(h0);
    ch0=cos(h0);
    Rp = q(1:4) * -quat_dot(q(1:4),q(5:8))/quat_dot(q(1:4),q(1:4));
    q(5:8) = q(5:8) + Rp;
    q(5:8) = q(5:8) * 2.0/h0;
    
    q(5:8) = q(5:8) * sh0;    
    Rp = q(1:4) * he*ch0*2.0/h0;
    q(5:8) = q(5:8) + Rp;
    q(5,1) = -he*sh0;
    
    q(1:4) = q(1:4) * sh0*2.0/h0;
    q(1,1) = ch0;
    
end
  
q = dquat_normalize(q);
    
end

