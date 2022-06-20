% add_group_lines(networks) draws boundaries on a plot at the specified intervals on a square FC. 
% The input variable, 'groups', is a 1 x n integer matrix  which
% specifies how many components are in each cluster. For example, 
% add_group_lines([3 3 3]) would draw lines at 3 and 6 component divisions.
function add_group_lines(groups)
regions_all = cumsum(groups);
r = sum(groups);
clines = regions_all+.5;
for i = 1:numel(clines)
    line([0 r+1],[clines(i) clines(i)],'Color','k')
    line([clines(i) clines(i)],[0 r+1],'Color','k')
end