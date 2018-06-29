%%
addpath('D:\\projects\\WavesCNR\\matlab\\dynamic\\dual_quaternions');
close all;
clear;
%%

C = [100 200 230 180;
     100 130 200 300]

hc = figure
Cplot = [C C(:,1)];
plot( Cplot(1,:),Cplot(2,:) )


%%
Rp = [1    1  0;
      2   1.3  0
      2.3   2  0
      1.8   3  0];
  
nrays = size( Rp, 1);
Rk =  [ 1 2 7];
%Rv = repmat(Rk,nrays,1) - Rp;
%for ii=1:nrays
%    Rv(ii,:)=Rv(ii,:)/norm(Rv(ii,:));
%end 
Rv = [0.0976    0.1952    0.9759;
      -0.0976    0.1952    0.9759;
      0 0 1;
      0 0 1];
%Rp(:,1) = Rp(:,1)*1.3;
 x=200;
 y=200;
 w = compute_weights(C,[x;y],1)
   
    

%%
% render all rays
hlerp = figure;
for ii=1:nrays
    render_ray(Rp(ii,:),Rv(ii,:),1.5,'b');
end
grid on
hold on
axis equal
view(-21,72);    
Rpinit = mean( Rp );
Rvinit = mean( Rv ); Rvinit = Rvinit/norm(Rvinit);
render_ray(Rpinit,Rvinit,2.5,'k');

dq_init = ray_to_dquat( Rpinit', Rvinit' );

%while true
    figure(hc)
    x=200;
    y=200;
    w = compute_weights(C,[x;y],1)
    %w=[1/4 1/4 1/4 1/4];
    figure( hlerp )
    
    %w=[1 1 1];
    %w=w/norm(w);  

    % Find all transforms
    dquats = zeros(8,nrays);
    for ii=1:nrays
        dquats(:,ii)=rays_min_transform( [Rpinit;Rp(ii,:)], [Rvinit;Rv(ii,:)] );
        
        dqa = ray_to_dquat( Rp(ii,:)', Rv(ii,:)' );
        for t=linspace(0,1,50);            
            render_frame( dquat_geodesic( dq_init, dquat_mult( dquats(:,ii), dq_init), t ), 0.2, 1.0 );
        end
    end

    % Interpolate
    dq_interp = dquat_mult( dquat_DIB(dquats,w,30), dq_init);
        
    [Rpi, Rvi] = dquat_to_ray(dq_interp);
        
    if exist('hlerpray','var')==1
        %delete(hlerpray)
        clear hlerpray
    end
        
    hlerpray = render_ray(Rpi',Rvi',2.0,'r');
    %hlerpray =render_frame(dq_interp,0.2,2.0);    
    drawnow
%end

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