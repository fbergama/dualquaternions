%%
close all;
figure;
fsize=15;
axis([-fsize fsize -fsize  fsize -fsize fsize ]);
grid on;
axis equal;


%%
q1 = dquat_from_angle_T(0,0,0,[0;1;1]);
q2 = dquat_from_angle_T(0,0.5,-0.5,[-3;-5;-2]);

%%

render_frame(q1,5.0,2.0);
render_frame(q2,5.0,1.0);

for t=linspace(0,1,30)
    render_frame( dquat_geodesic(q1,q2,t), 0.2, 1.0 );
end


%%
q12 = dquat_mult( q2, dquat_invert(q1) );
q12 = q12*-1;
dquat_norm(q12)
for t=linspace(0,1,30)
    qp = dquat_geodesic(dquat_eye(),q12,t) ;
    render_frame(dquat_mult(qp,q1),1.0,2.0);
end


