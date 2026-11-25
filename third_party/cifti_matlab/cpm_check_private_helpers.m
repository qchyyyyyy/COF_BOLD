function cpm_check_private_helpers()
%CPM_CHECK_PRIVATE_HELPERS Test repaired FieldTrip warning dependencies.
% Release integration test, covered by the repository's author-code MIT licence.
% Located beside ft_read_cifti so MATLAB can resolve that reader's private API.
s = struct('nested',struct('value',42));
assert(issubfield(s,'nested.value') && ~issubfield(s,'nested.absent'));
assert(getsubfield(s,'nested.value') == 42);
s = setsubfield(s,'nested.value',43);
assert(s.nested.value == 43);
global ft_default
oldDefault = ft_default;
oldWarning = warning('off','CPM:DependencyProbe');
cleanup = onCleanup(@() restore_warning(oldDefault,oldWarning)); %#ok<NASGU>
ft_warning('CPM:DependencyProbe','Synthetic dependency probe',0);
end

function restore_warning(oldDefault,oldWarning)
global ft_default
ft_default = oldDefault;
warning(oldWarning);
end
