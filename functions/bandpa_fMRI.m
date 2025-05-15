function [bandpasSig, hilbertSig, ampSig, phaseSig] = bandpa_fMRI(sigValid,fs,fLow,fHigh,isdemean) 
% bandpass the fMRI data
[b, a]=butter(4,[fLow fHigh]/(fs/2),'bandpass');

badChannels = find(sum(sigValid,2)==0) ;
bandpasSig = zeros(size(sigValid)) ;
hilbertSig = nan(size(sigValid)) ;
phaseSig = nan(size(sigValid)) ;
ampSig = nan(size(sigValid)) ;

for iNode = setdiff(1:size(sigValid,1),badChannels)
    sigOriMean = mean(sigValid(iNode,:));
    sigOriDemean = sigValid(iNode,:)- sigOriMean;
    bandpasSig(iNode,:) = filter(b,a,double(sigOriDemean));
    if ~isdemean
        bandpasSig(iNode,:) = bandpasSig(iNode,:) + sigOriMean;
    end
    hilbertSig(iNode,:) = hilbert(bandpasSig(iNode,:)) ;
    phaseSig(iNode,:) = angle(hilbertSig(iNode,:)) ;
    ampSig(iNode,:) = abs(hilbertSig(iNode,:)) ;
end
