function qt = dquat_LogI( q )
%DQUAT_LOGI

%{
if q(1,1)<0
    q(1:8) = q(1:8)*-1;    
end
%}

qt = q;
h0 = acos( qt(1,1) );

if h0*h0<1E-12
    % Small angle approximation: sin(h0)=h0, cos(h0)=1
    qt(1,1)=0;
    qt(5,1)=0;
    qt=qt.*0.5;
else
    qt(1,1)=0;
    ish0=1.0/quat_norm(qt(1:4));    
    qt(1:4) = quat_normalize( qt(1:4) );
    he=-qt(5,1)*ish0;
    qt(5,1)=0;
    
    Rp = qt(1:4) * -quat_dot(qt(1:4),qt(5:8))/quat_dot(qt(1:4),qt(1:4));
    qt(5:8) = qt(5:8) + Rp;
    qt(5:8) = qt(5:8)*ish0;
    qt(5:8) = qt(5:8)*h0;
    
    Rp = qt(1:4);
    Rp = Rp*he;
    qt(5:8) = qt(5:8) + Rp;
    qt(5:8) = qt(5:8) * 0.5;
    
    qt(1:4) = qt(1:4)*h0*0.5;
end

       
    
end

