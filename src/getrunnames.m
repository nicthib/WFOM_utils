function runnames = getrunnames(mouse,runs,allrunnames)
runnames = allrunnames(cellfun(@(s) ~isempty(regexp(s,runs)),allrunnames)); % Parse out runs that end with B,C,D
runnames = runnames(cellfun(@(s) ~isempty(regexp(s,mouse)),runnames)); % Parse out runs that are 'mouse'