%%
addpath('D:\\projects\\WavesCNR\\matlab\\dynamic\\dual_quaternions');
close all;
clear;

%%
Rp = [1    1  -1;
      1   1  -1];
  
% Rv = [0   0  1;
%      -1   0  1];
Rv = [0   0  1;
      0   1  1];
  
nrays = size( Rp, 1);
for ii=1:nrays
    Rv(ii,:)=Rv(ii,:)/norm(Rv(ii,:));
end 

%%
% close all;
% q1 = dquat_eye();
% q2 = dquat_from_screw_motion(pi/2,3,[1;0;0],[0;0;0])
% 
% for t=linspace(0,1,10);
%     render_frame( dquat_geodesic( q1, dquat_mult( q2, q1), t ), 0.2, 1.0 );
% end
% grid on
% hold on
% axis equal
% view(-21,72); 
% xlabel('X');
% ylabel('Y');
% zlabel('Z');


%%
% render all rays
filename = 'test2.gif';
for t=linspace(0,1,50)
    disp(t)
    close all;
    fig = figure('position',[100 100 800 500]);
    set(fig,'Renderer','zbuffer')
    
    Rp(2,1)=t;
    Rv(2,1)=sin(t*pi);
    for ii=1:nrays
    Rv(ii,:)=Rv(ii,:)/norm(Rv(ii,:));
    end 

    for ii=1:nrays
        render_ray(Rp(ii,:),Rv(ii,:),1.5,'b');
    end
    grid on
    hold on
    axis equal
    view(-33,16);    
    axis( [-2 2 -2 2 -2 2] );

    [mindq,midpoint,pta,ptb] = rays_min_transform( Rp, Rv );
    dq_init = ray_frame_to_dquat( pta', Rv(1,:)', (ptb-pta)' );
    for t2=linspace(0,1,20);
        render_frame( dquat_geodesic( dq_init, dquat_mult( mindq, dq_init), t2 ), 0.2, 1.0 );
    end

    scatter3(midpoint(1),midpoint(2),midpoint(3));
    drawnow;
    
    frame = getframe(fig);
    im = frame2im(frame);    
    if t==0
        [imind,cm] = rgb2ind(im,256,'nodither');
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime', 0);        
    else
        imind = rgb2ind(im,cm);
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
end

%%
