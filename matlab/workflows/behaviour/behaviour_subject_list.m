function subjects = behaviour_subject_list(file)
% Read historical/current participant-list spelling without changing order.
S = load(file);
upper = isfield(S,'SubjectList');
lower = isfield(S,'subjectList');
assert(upper || lower,'Behaviour:SubjectField','Expected SubjectList or subjectList in %s.',file);
if upper && lower
    assert(isequal(string(S.SubjectList(:)),string(S.subjectList(:))), ...
        'Behaviour:SubjectConflict','SubjectList and subjectList disagree in %s.',file);
end
if upper,subjects=S.SubjectList;else,subjects=S.subjectList;end
end
