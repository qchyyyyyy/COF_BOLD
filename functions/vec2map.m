function Map = vec2map(vec,mask)
if length(size(mask))==2
    [M,N] = size(mask);
    Map = zeros(M*N,size(vec,2));
    Map(mask(:)==1,:)=vec;
    Map = reshape(Map,M,N,size(vec,2));

elseif length(size(mask))==3

    [M,N,S] = size(mask);
    Map = zeros(M*N*S,size(vec,2));
    Map(mask(:)==1,:)=vec;
    Map = reshape(Map,M,N,S,size(vec,2));
end