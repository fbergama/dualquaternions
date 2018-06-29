function dq = ray_frame_to_dquat( Rp, Rv1, Rv2 )

% choose a frame
v2 = Rv2/norm(Rv2);
while abs(dot( v2, Rv1 ))>1E-3
    v2 = rand(3,1); v2 = v2/norm(v2);
end
v3 = cross(Rv1,v2);
v3=v3/norm(v3);
v2=cross(v3,Rv1);

% compute RT
RT = [ Rv1 , v2 , v3 ,  Rp;
         0   0    0     1 ];

% create dual quaternion
dq = dquat_from_RT(RT);

end

