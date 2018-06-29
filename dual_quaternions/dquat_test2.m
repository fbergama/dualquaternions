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

qz = dquat_from_angle_T(-0.5,0,0,[-5;-3;3]);
render_frame(qz,3);
qdv = dquat_from_angle_T(0,0,0,[1;2;0]);

for ii=1:30
    qz = dquat_mult( qdv, qz)
    render_frame(qz,0.8);
end

