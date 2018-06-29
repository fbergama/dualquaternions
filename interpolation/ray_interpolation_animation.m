%%
addpath('D:\\projects\\WavesCNR\\matlab\\dynamic\\dual_quaternions\\');
close all;
clear;
filename = 'test.gif'
%%

C = [100 200 230 180;
     100 130 200 300]


%%
Rp = [1    1  0;
      2   1.3  0
      2.8   2  0
      1.8   3  0];
  
nrays = size( Rp, 1);
Rk =  [ -1 2 7];
Rv = repmat(Rk,nrays,1) - Rp;
%Rv(2,:) = [0 0 1];
for ii=1:nrays
    Rv(ii,:)=(Rv(ii,:))/norm(Rv(ii,:));
end 
Rp(:,1) = Rp(:,1)*1.3;

    
Rpinit = mean( Rp );
Rvinit = mean( Rv ); 
Rvinit = Rvinit/norm(Rvinit);
dq_init = ray_to_dquat( Rpinit', Rvinit' );
firstloop = true;

for t=linspace(0,1,100)
    cs = [180;210] + [70*t*sin(t*10*pi);90*t*cos(t*10*pi)];
    w = compute_weights(C,cs,1);
    
    % Find all transforms
    dquats = zeros(8,nrays);
    for ii=1:nrays
        dquats(:,ii)=rays_min_transform( [Rpinit;Rp(ii,:)], [Rvinit;Rv(ii,:)] ) ;
    end

    % Interpolate
    dq_interp = dquat_mult( dquat_DIB(dquats,w,10), dq_init);
        
    [Rpi, Rvi] = dquat_to_ray(dq_interp);
    %render_frame(dq_interp,0.2,2.0);    
    
    % Render code
    if ~firstloop
        delete(hray);
        delete(hcode);
    end
    
    close all;
    fig = figure('position',[100 100 800 500]);
    set(fig,'Renderer','zbuffer')
    subplot(1,2,1);
    Cplot = [C C(:,1)];
    plot( Cplot(1,:),Cplot(2,:) );
    hold on
    grid on
    hcode = scatter( cs(1),cs(2),'oR');
    axis( [50 280 50 350] );
    title('Code');
    
    % Render rays
    subplot(1,2,2);
    for ii=1:nrays
        render_ray(Rp(ii,:),Rv(ii,:),1.5,'b');
    end
    hray = render_ray(Rpi',Rvi',2.0,'r');
    view(-82,52);
    grid on
    hold on
    axis( [0 3.3 0 3.3 0 3] );
    title('Ray');
    
    drawnow    
    % Save gif
    firstloop = false;
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
%{
Re = cross(Rv,Rp);
midpoint = rays_intersec(Rv,Re);
scatter3( midpoint(1), midpoint(2), midpoint(3) );


pt_ra = Rp(1,:) + Rv(1,:)*dot( midpoint-Rp(1,:)', Rv(1,:));
pt_rb = Rp(2,:) + Rv(2,:)*dot( midpoint-Rp(2,:)', Rv(2,:));
scatter3( pt_ra(1), pt_ra(2), pt_ra(3) );
scatter3( pt_rb(1), pt_rb(2), pt_rb(3) );


vaxis = pt_rb-pt_ra; vaxis=-vaxis'/norm(vaxis);
RT_ra = [ Rv(1,:)' , vaxis , cross(Rv(1,:)',vaxis) ,  pt_ra';
             0         0                0               1 ];

dquat_ra = dquat_from_RT(RT_ra);

RT_rb = [ Rv(2,:)' , vaxis , cross(Rv(2,:)',vaxis) ,  pt_rb';
             0         0                0               1 ];
          
dquat_rb = dquat_from_RT(RT_rb);         


render_frame(dquat_ra,0.2,2.0);
render_frame(dquat_rb,0.2,2.0);


quat_a2b = dquat_mult( dquat_rb, dquat_invert(dquat_ra) );

for t=linspace(0,1,20)
    render_frame(dquat_geodesic(dquat_ra,dquat_mult(quat_a2b,dquat_ra),t),0.2,1.0);
end




RT_ra = [ Rv(1,:)' , vaxis , cross(Rv(1,:)',vaxis) ,  (pt_ra + Rv(1,:))';
             0         0                0               1 ];

dquat_ra = dquat_from_RT(RT_ra);

for t=linspace(0,1,20)
    render_frame(dquat_geodesic(dquat_ra,dquat_mult(quat_a2b,dquat_ra),t),0.2,1.0);
end
%}