function qt = dquat_LogI2( q )
%DQUAT_LOGI2

[th_h, s]=dual_to_screw_form(q);

qt = dquat_mult_dual_number( s./4, th_h );

end

