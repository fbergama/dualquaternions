%%
close all;clc;clear;
addpath('../../extlibs');

size = 20;
axis([-size size -size size -size size]);
axis equal;
hold on;
xlabel('x');
ylabel('y');
zlabel('z');
grid on;



%%
% geodesic
% q = dquat_from_angle_T(0,0,0,[-5;0;-6]);
% p = dquat_from_angle_T(0,pi,pi,[5;0;6]);
% 
% render_frame(q,3);
% render_frame(p,3);
% 
% 
% N=10;
% tvals = linspace(0,1,N);
% for t=tvals
%     pt = dquat_geodesic(p,q,t);
%     render_frame(pt,1);
%     
%     m = dquat_DLB(p,q,t);
%     render_frame(m,2);
% end

%{
qz = dquat_from_angle_T(0.1,0.0,0.0,[30;30;0]);
%qz = dquat_eye();
P = dquat_make_covariance(qz,0.1,0.1,0.1,...
                              0.5,0.5,0.5)
eig(P)
N=50;
for ii=1:N
    v = mvnrnd(zeros(8,1),P)';
    qc = dquat_Exp(qz, v );
    
    RT = dquat_to_RT(qc);
    T = RT(1:3,4);
    if norm(T-[30;30;0])>50 
        render_frame( qc , 30 );
    else
        render_frame( qc , 3 );
    end
end
render_frame( qz , 6 );
%}

%%
% sample distribution

qz = dquat_from_angle_T(0,1.9,0.8,[0;15;0]);
P = dquat_make_covariance(dquat_eye(),...
                                0,0,0,...
                              0.5,0.2,0.2)
                          
render_frame( dquat_eye() , 6 );
render_frame( qz , 6 );
N=50;
for ii=1:N
    q = dquat_sample_normal(qz,P);
    render_frame( q , 3 );    
end


%%
% empirical covariance
%{
N=30;
dqs = zeros(8,N);
w=ones(N,1)/N;
for ii=1:N
    dqs(:,ii)=dquat_normalize( dquat_from_angle_T(rand(),rand(),rand(),[rand();rand();rand()]) );
    render_frame( dqs(:,ii) , 2 );
end

% dm = dquat_DIB(dqs,w,10);
% render_frame( dm , 8 );
% P = dquat_empirical_covariance(dqs, dm);
% [V,D] = eig(P,'nobalance')
% for ii=1:N
%     qt = mvnrnd(zeros(8,1),P)';    
%     qz = dquat_Exp(dm,qt);
%     %render_frame( qz , 3 );
% end

dm = dquat_DIB(dqs,w,30);
render_frame( dm , 8 );

P = dquat_empirical_covariance(dqs, dm)
%[V,D] = eig(P,'nobalance')

z = dquat_from_angle_T(pi,pi/2,0.2,[-5;0;-7]);
render_frame( z , 8 );
Pt = dquat_parallel_transport(dm,z,P,30);

for ii=1:N
    qt = mvnrnd(zeros(8,1),Pt)';
    qz = dquat_Exp(z,qt);
    render_frame( qz , 3 );
end

%}
%%
% rotq = dquat_from_angle_T(0,0,0.8,[0;0;0]);
% traxq = dquat_from_angle_T(0,0,0,[5;0;0]);
% trayq = dquat_from_angle_T(0,0,0,[0;5;0]);
% 
% %render_frame(rotq,5);
% 
% render_frame(dquat_mult(traxq,rotq),5);
% 
% RT_rotq = dquat_to_RT(rotq);
% RT_traxq = dquat_to_RT(traxq);
% 
% render_frame(RT_traxq*RT_rotq,10);
