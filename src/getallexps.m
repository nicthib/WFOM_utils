function expnames = getallexps(root,vars)
% expnames = getallexps(root, vars) returns the experiment names in the given
% directory root that meets the variable requirements in the cell string vars.
tmp = dir(fullfile(root,'**/*.mat'));
files = {tmp.name};
expnames = {};
for i = 1:numel(files)
    try
        if sum(cellfun(@(x) ismember(x,who('-file',files{i})), vars)) == numel(vars)
            expnames{end+1} = files{i};
        end
    catch
    end
end
