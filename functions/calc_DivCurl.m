function [DivMap, CurlMap] = calc_DivCurl(Ux,Uy,V_mask,filtSigma,do_zscore,do_binarize)

if nargin <5
    do_zscore = 0;
end
if nargin <6
    do_binarize = 0;
end
filtWidth = round(3*filtSigma);
imageFilter=fspecial('gaussian',filtWidth,filtSigma);

T = size(Ux,3);
for t = T:-1:1
    U = sqrt(Ux(:,:,t).^2+Uy(:,:,t).^2);

    div_oneframe= divergence(Ux(:,:,t),Uy(:,:,t));
    div_oneframe(V_mask==0)=nan;

    [curl_oneframe,~]= curl(Ux(:,:,t)./U,Uy(:,:,t)./U);
    curl_oneframe(V_mask==0)=nan;

    if nargin <4
        DivMap(:,:,t) = div_oneframe;
        CurlMap(:,:,t) = curl_oneframe;
    else
        DivMap(:,:,t) = nanconv(div_oneframe,imageFilter);
        CurlMap(:,:,t) = nanconv(curl_oneframe,imageFilter);
    end
end

CurlMap(repmat(V_mask,1,1,T)==0)=nan;
DivMap(repmat(V_mask,1,1,T)==0)=nan;

if do_zscore
    DivMean = mean(DivMap,'all','omitnan');
    DivStd = std(DivMap,[],'all','omitnan');
    DivMap = (DivMap - DivMean)/DivStd;
    CurlMean = mean(CurlMap,'all','omitnan');
    CurlStd = std(CurlMap,[],'all','omitnan');
    CurlMap = (CurlMap - CurlMean)/CurlStd;
end

if do_binarize
    DivMap(DivMap<k & DivMap>-k)=0;
    DivMap(DivMap>0)=1;
    DivMap(DivMap<0)=-1;
    CurlMap(CurlMap<k & CurlMap>-k)=0;
    CurlMap(CurlMap>0)=1;
    CurlMap(CurlMap<0)=-1;
end