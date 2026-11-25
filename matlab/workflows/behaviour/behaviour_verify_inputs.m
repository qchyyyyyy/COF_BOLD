function behaviour_verify_inputs(groupDir,manifestFile)
% Require the exact frozen group inputs when replaying historical estimates.
manifest = jsondecode(fileread(manifestFile));
for k = 1:numel(manifest.files)
    item = manifest.files(k);
    file = fullfile(groupDir,item.name);
    assert(isfile(file),'Behaviour:MissingInput','Missing %s.',file);
    info = dir(file);
    assert(info.bytes==item.bytes,'Behaviour:InputVersion','Size differs for %s.',item.name);
    md = java.security.MessageDigest.getInstance('SHA-256');
    fid = fopen(file,'rb');
    assert(fid>=0,'Behaviour:InputRead','Cannot read %s.',file);
    cleanup = onCleanup(@()fclose(fid));
    while ~feof(fid)
        block = fread(fid,1048576,'*uint8');
        md.update(block);
    end
    digest = lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
    assert(strcmp(digest,item.sha256),'Behaviour:InputVersion','Hash differs for %s.',item.name);
    clear cleanup
end
end
