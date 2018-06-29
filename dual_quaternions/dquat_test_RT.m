
%%


q = dquat_from_angle_T( 3,0,0,[0;0;0]);


%%
figure;
axis([-5 5 -5 5 -5 5 ]);
grid on;

for ii=linspace(0,6.28,20)
    q = dquat_normalize( dquat_from_angle_T( ii,0,0,[3;3;0]) );
    dquat_norm(q)
    render_frame(q,1);
end
