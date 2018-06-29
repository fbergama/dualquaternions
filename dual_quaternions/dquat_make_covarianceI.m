function P = dquat_make_covarianceI( center, anglespan, tspan )
%DQUAT_MAKE_COVARIANCEI

N=100;
qts = zeros(8,N);
for ii=1:N
    T=center+mvnrnd([0;0;0],diag(tspan))';
    angles=mvnrnd([0;0;0],diag(anglespan))';
    qts(:,ii) = dquat_from_angle_T(angles(1),...
                                   angles(2),...
                                   angles(3),...
                                   T);
    %render_frame(qts(:,ii),1);
end

[~, P]=dquat_empirical_covarianceI(qts);

end

