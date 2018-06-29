function h = render_ray( Rp, Rv, linewidth, color, raylen )

if nargin<3
    linewidth = 1.0;
end

if nargin<4
    color = 'b';
end

if nargin<5
    raylen = 5;
end

src = Rp - Rv*raylen*0.5;
dst = Rp + Rv*raylen*0.5;
h = line( [src(1,1) dst(1,1)],[src(1,2) dst(1,2)],[src(1,3) dst(1,3)], 'Color',color,'LineWidth',linewidth,'LineSmoothing','on' );


end

