function [DivMap, CurlMap] = velocity_derivatives(Ux,Uy,V_mask,filtSigma,do_zscore,do_binarize,n_std)

if nargin <5
    do_zscore = 0;
end
if nargin <6
    do_binarize = 0;
end
assert(isequal(size(Ux),size(Uy)) && isequal(size(Ux,[1 2]),size(V_mask)), 'cof:Shape','Velocity/mask dimensions disagree.');
if nargin >= 4
    validateattributes(filtSigma,{'numeric'},{'scalar','positive','finite'});
    filtWidth = round(3*filtSigma);
    imageFilter=fspecial('gaussian',filtWidth,filtSigma);
end

T = size(Ux,3);
for t = T:-1:1
    U = sqrt(Ux(:,:,t).^2+Uy(:,:,t).^2);

    div_oneframe= divergence(Ux(:,:,t),Uy(:,:,t));
    div_oneframe(V_mask==0)=nan;

    [curl_oneframe,~]= curl(Ux(:,:,t)./U,Uy(:,:,t)./U);
    curl_oneframe(V_mask==0)=nan;

    if nargin<4
        DivMap(:,:,t) = div_oneframe;
        CurlMap(:,:,t) = curl_oneframe;
    else
        DivMap(:,:,t) = nanconv(div_oneframe,imageFilter);
        CurlMap(:,:,t) = nanconv(curl_oneframe,imageFilter);
    end
end

CurlMap(repmat(V_mask,1,1,T)==0)=nan;
DivMap(repmat(V_mask,1,1,T)==0)=nan;


if nargin>4 && do_zscore
    DivMean = mean(DivMap,'all','omitnan');
    DivStd = std(DivMap,[],'all','omitnan');
    DivMap = (DivMap - DivMean)/DivStd;
    CurlMean = mean(CurlMap,'all','omitnan');
    CurlStd = std(CurlMap,[],'all','omitnan');
    CurlMap = (CurlMap - CurlMean)/CurlStd;
end

if nargin>5 && ~exist('n_std','var')
    n_std = 2;
end

if do_zscore && do_binarize
    DivMap(DivMap<n_std & DivMap>-n_std)=0;
    DivMap(DivMap>0)=1;
    DivMap(DivMap<0)=-1;
    CurlMap(CurlMap<n_std & CurlMap>-n_std)=0;
    CurlMap(CurlMap>0)=1;
    CurlMap(CurlMap<0)=-1;
elseif ~do_zscore && do_binarize
    DivMean = mean(DivMap,'all','omitnan');
    DivStd = std(DivMap,[],'all','omitnan');
    DivMap(DivMap<DivMean+n_std*DivStd & DivMap>DivMean-n_std*DivStd)=0;
    DivMap(DivMap>0)=1;
    DivMap(DivMap<0)=-1;

    CurlMean = mean(CurlMap,'all','omitnan');
    CurlStd = std(CurlMap,[],'all','omitnan');
    CurlMap(CurlMap<CurlMean+n_std*CurlStd & CurlMap>CurlMean-n_std*CurlStd)=0;
    CurlMap(CurlMap>0)=1;
    CurlMap(CurlMap<0)=-1;
end