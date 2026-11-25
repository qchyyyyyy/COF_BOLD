function map = grid_to_vector(map,mask)
[M,N] = size(mask);
map = reshape(map,M*N,[]);
map(mask(:)~=1,:)=[];