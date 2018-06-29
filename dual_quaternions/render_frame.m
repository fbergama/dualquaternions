function [ h ] = render_frame( q, scale, linewidth )
%RENDER_FRAME

if nargin<3
    linewidth=2;
end

if size(q,2)==1
    RT = dquat_to_RT( q );
else
    RT = q;
end

src = RT(1:3,4);
dst1 = src+RT(1:3,1)*scale;
dst2 = src+RT(1:3,2)*scale;
dst3 = src+RT(1:3,3)*scale;

line( [src(1,1) dst1(1,1)],[src(2,1) dst1(2,1)],[src(3,1) dst1(3,1)], 'Color','r','LineWidth',linewidth,'LineSmoothing','on' );
line( [src(1,1) dst2(1,1)],[src(2,1) dst2(2,1)],[src(3,1) dst2(3,1)], 'Color','g','LineWidth',linewidth,'LineSmoothing','on' );
line( [src(1,1) dst3(1,1)],[src(2,1) dst3(2,1)],[src(3,1) dst3(3,1)], 'Color','b','LineWidth',linewidth,'LineSmoothing','on' );


h=0;

end

