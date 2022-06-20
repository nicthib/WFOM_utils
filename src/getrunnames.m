function runnames = getrunnames(mouse,runs,allrunnames)
% Run filter function. Input is a cell string vector.

runnames = allrunnames(cellfun(@(s) ~isempty(regexp(s,runs)),allrunnames)); % Parse out runs that end with any of the letters in runs
runnames = runnames(cellfun(@(s) ~isempty(regexp(s,mouse)),runnames)); % Parse out runs that contain mouse ID.