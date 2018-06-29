function render_cov( q, P, sz )
%RENDER_COV

for ii=1:50
    xt = mvnrnd(zeros(8,1),P)';
    xx = dquat_ExpI(xt);    
    
    % move xx back to q reference frame
    render_frame(dquat_mult(q,xx),sz);
end

end

