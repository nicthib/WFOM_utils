function [stj,sth] = state_dist(results)
if numel(results) == 5
    stj = cell(5,1); sth = cell(5,1);
    for i = 1:numel(results)
        for j = 1:numel(results{i})
            if isfield(results{i}{j},'stj') && numel(results{i}{j}.stj)>1000
                stj{i} = [stj{i};results{i}{j}.stj];
                sth{i} = [sth{i};results{i}{j}.sth];
            end
        end
    end
else
    stj = []; sth = [];
    for i = 1:numel(results)        
        if isfield(results{i},'stj') 
            stj = [stj;results{i}.stj];
            sth = [sth;results{i}.sth];
        end
    end
end

if iscell(stj)
    cj = cellfun(@(x) histcounts(x,1:9)'/sum(x),stj,'UniformOutput',false);
    ch = cellfun(@(x) histcounts(x,1:9)'/sum(x),sth,'UniformOutput',false);
else
    cj = histcounts(stj);
    ch = histcounts(sth);
end
close all; figure
subplot(211); bar([cj{:}]')
subplot(212); bar([ch{:}]')