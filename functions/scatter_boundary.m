function ss = scatter_boundary(oi_x,oi_y,Size)
if nargin <3
    ss = scatter(oi_x,oi_y,'.');
else
    ss = scatter(oi_x,oi_y,Size,'.');
end
ss.MarkerEdgeColor = [0.4 0.4 0.4];
% ss.LineWidth=0.01;