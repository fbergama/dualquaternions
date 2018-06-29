%%
close all;clc;clear;
addpath('../../extlibs');

size = 10;
axis([-size size -size size -size size]);
axis equal;
hold on;
xlabel('x');
ylabel('y');
zlabel('z');
grid on;



%% 
N=20;
qts = zeros(8,N);
anglespan=[pi;pi;pi];
tc=[3;4;1];
tspan=[pi;pi;pi];
for ii=1:N
    T=tc+mvnrnd([0;0;0],diag(tspan))';
    angles=mvnrnd([0;0;0],diag(tspan))';
    qts(:,ii) = dquat_from_angle_T(angles(1),...
                                   angles(2),...
                                   angles(3),...
                                   T);
    render_frame(qts(:,ii),1);
end

[mu, P]=dquat_empirical_covarianceI(qts);
render_frame(mu,10);
P


