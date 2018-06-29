function q = dquat_Exp( p, qt )
%DQUAT_EXP

% Exp_p (qt) = Exp_1 ( qt p* ) p

%q = dquat_mult( dquat_ExpI( dquat_mult(qt, dquat_conj_quat(p)) ) , p);
q = dquat_mult( dquat_ExpI( qt ) , p);

end

