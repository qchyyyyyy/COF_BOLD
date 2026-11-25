function mesh = meg_grid_zscore(mesh, mask)
%MEG_GRID_ZSCORE Match Project 01 temporal grid-point z-scoring.

[nY,nX,nT] = size(mesh);
flat = reshape(double(mesh),nY*nX,nT);
active = mask(:) ~= 0;
mu = mean(flat(active,:),2);
sd = std(flat(active,:),0,2);
if any(~isfinite(sd) | sd <= eps)
    error('MEG:WorkflowGridVariance', 'An active grid point has zero variance.');
end
flat(active,:) = bsxfun(@rdivide,bsxfun(@minus,flat(active,:),mu),sd);
flat(~active,:) = NaN;
mesh = single(reshape(flat,nY,nX,nT));
end
