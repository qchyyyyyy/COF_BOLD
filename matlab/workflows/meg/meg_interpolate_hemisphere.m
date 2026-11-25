function meshData = meg_interpolate_hemisphere(sourceData, x, y, mask)
%MEG_INTERPOLATE_HEMISPHERE Map one hemisphere to the Project 01 201x251 grid.

x = double(x(:)); y = double(y(:)); sourceData = double(sourceData);
if size(sourceData,1) ~= numel(x) || numel(x) ~= numel(y)
    error('MEG:WorkflowMeshShape', 'Source data and coordinates do not align.');
end
xTarget = -250:2:250;
yTarget = -200:2:200;
if ~isequal(size(mask), [numel(yTarget) numel(xTarget)])
    error('MEG:WorkflowMaskShape', 'Expected a 201x251 hemisphere mask.');
end
[xi,yi] = meshgrid(xTarget,yTarget);
[uniquePoints,~,group] = unique([x y],'rows');
counts = accumarray(group,1,[size(uniquePoints,1) 1]);
averager = sparse(group,(1:numel(x))',1./counts(group),size(uniquePoints,1),numel(x));
values = averager * sourceData;
active = find(mask(:) ~= 0);
tri = delaunayTriangulation(uniquePoints);
[simplex,bary] = pointLocation(tri,[xi(active) yi(active)]);
inside = ~isnan(simplex);
vertices = tri.ConnectivityList(simplex(inside),:);
rows = repmat(active(inside),size(vertices,2),1);
weights = bary(inside,:) .* mask(active(inside));
interpolator = sparse(rows,vertices(:),weights(:),numel(mask),size(uniquePoints,1));
flat = interpolator * values;
outside = active(~inside);
if ~isempty(outside)
    extrap = scatteredInterpolant(uniquePoints(:,1),uniquePoints(:,2),values(:,1));
    for t = 1:size(values,2)
        extrap.Values = values(:,t);
        flat(outside,t) = extrap(xi(outside),yi(outside)) .* mask(outside);
    end
end
meshData = reshape(flat,size(mask,1),size(mask,2),size(sourceData,2));
end
