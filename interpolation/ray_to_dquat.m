function dq = ray_to_dquat( Rp, Rv )
%RAY_TO_DQUAT

% choose a frame
v2 = (rand(3,1)-[0.5 ; 0.5 ; 0.5]); v2 = v2/norm(v2);
while abs(dot( v2, Rv ))<1E-3
    v2 = rand(3,1); v2 = v2/norm(v2);
end
v3 = cross(Rv,v2);
v3=v3/norm(v3);
v2=cross(v3,Rv);

% compute RT
RT = [ Rv , v2 , v3 ,  Rp;
        0   0    0     1 ];

% create dual quaternion
dq = dquat_from_RT(RT);

end

