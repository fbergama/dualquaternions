function qt = dquat_Log( p, q )
%DQUAT_LOG

% Log_p (q) = Log_1 ( q p* ) p
%qt = dquat_mult( dquat_LogI( dquat_mult(q,dquat_conj_quat(p)) ) , p);
%qt = dquat_mult( p, dquat_LogI( dquat_mult(dquat_conj_quat(p),q) ));
qt = dquat_LogI( dquat_mult(q,dquat_conj_quat(p)) ) ;

end

