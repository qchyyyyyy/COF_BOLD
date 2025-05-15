function map = map2vec(map,mask)
[M,N] = size(mask);
map = reshape(map,M*N,[]);
map(mask(:)==0,:)=[];