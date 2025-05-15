function q = FDR_correction(p)
[~,index] = sort(p,'descend');
p = p(index);  
q = zeros(size(p));
q(1:end) = length(p):-1:1;
q = p*length(p)./q;
for i = 1:length(p)-1
    if q(i)<q(i+1)
        q(i+1)=q(i);
    end
end
q(index)=q;