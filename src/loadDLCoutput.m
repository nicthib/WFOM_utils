function [DLCoutputmod] = loadDLCoutput(inputfile,threshold)

% load raw DLCoutput.mat file created by 'h5_to_mat_BR.ipynb'
% specify p-value threshold for likelihood of feature
% specify type of mod: replace with Nan or interpolate using linspace

raw = load(inputfile);

DLCoutputmod = raw; % won't modify input file

raw_fn = fieldnames(raw);
for k = 1:numel(raw_fn) % iterate through each h5 file in output structure
    run = raw.(raw_fn{k});
    runmod = run;
    features = fieldnames(run);
    for m = 2:numel(features) % iterate through each body part
        pvals = run.(features{m})(:,3); % list of all pvals for each part
        pvals_mod = zeros(length(pvals),1);
        for n = 1:length(pvals) % iterate through and replace with nan
            if pvals(n) < threshold
                pvals_mod(n) = nan;
            else
                pvals_mod(n) = pvals(n);
            end
        end
        runmod.(features{m})(:,3) = pvals_mod; 
    end
    DLCoutputmod.(raw_fn{k}) = runmod;
end