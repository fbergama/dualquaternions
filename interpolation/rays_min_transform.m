function [dq, midpoint, pt_ra, pt_rb] = rays_min_transform( Rps, Rvs )
%RAYS_MIN_TRANSFORM 

vaxis_r = cross( Rvs(1,:), Rvs(2,:) );
norm_vaxis_r = norm(vaxis_r);

if norm_vaxis_r < 1E-5
    % Rays are almost parallel, so there is only
    % a translation component
    fprintf('Rays are parallel\n');
    
    pt_ra = Rps(1,:);
    pt_rb = Rps(2,:) + Rvs(2,:)*dot( Rps(1,:)-Rps(2,:), Rvs(2,:));

    tvec = pt_rb-pt_ra;   
    midpoint = 0.5*( pt_rb+pt_ra );
    dq = [ quat_eye(); 0; 0.5.*tvec' ]; % output dual quaternion
    return;
end


% Since rays are not parallel, there exists a single well defined point
% 'midpoint' minimizing the distance to both the rays
%
Res = cross(Rvs,Rps);
midpoint = rays_intersec(Rvs,Res);

pt_ra = Rps(1,:) + Rvs(1,:)*dot( midpoint-Rps(1,:)', Rvs(1,:));
pt_rb = Rps(2,:) + Rvs(2,:)*dot( midpoint-Rps(2,:)', Rvs(2,:));

vaxis_t = pt_rb-pt_ra;
norm_vaxis_t = norm(vaxis_t);


axis_point = midpoint;
axis_direction = vaxis_r';
anglesign=1;
if norm_vaxis_t > norm_vaxis_r
    axis_direction = vaxis_t';
    if dot( vaxis_t, vaxis_r ) < 0
        anglesign = -1;
    end
    %fprintf('Choosing Pa-Pb as axis\n');
end

axis_direction = axis_direction/norm(axis_direction);
translation_amount = norm_vaxis_t;
angle = anglesign * asin( norm_vaxis_r );
dq = dquat_from_screw_motion(angle,translation_amount,axis_direction,axis_point);

