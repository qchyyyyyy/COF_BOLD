function [orderedScores, lengths, subjects] = reorder_participant_blocks(scores, blockLengths, subjectList, sourceOrder)
%REORDER_SUBJECT_BLOCKS Align concatenated scores to an explicit subject list.
% blockLengths is indexed in subjectList order. sourceOrder specifies whether
% scores was concatenated in 'forward' or 'reverse' subjectList order. Within-
% subject time order is preserved. No data-driven permutation is estimated.
validateattributes(scores,{'numeric'},{'vector','finite'});
validateattributes(blockLengths,{'numeric'},{'vector','integer','positive','finite'});
lengths=double(blockLengths(:));subjects=cellstr(string(subjectList(:)));
assert(numel(lengths)==numel(subjects),'cof:BlockCount','Lengths and subjects differ in count.');
assert(numel(unique(subjects))==numel(subjects),'cof:DuplicateSubject','Subject list must be unique.');
assert(sum(lengths)==numel(scores),'cof:BlockLength','Lengths do not partition scores.');
switch char(sourceOrder)
    case 'forward'
        orderedScores=scores(:);
    case 'reverse'
        blocks=mat2cell(scores(:),flipud(lengths));
        orderedScores=vertcat(blocks{end:-1:1});
    otherwise
        error('cof:SourceOrder','Specify forward or reverse explicitly.');
end
end
