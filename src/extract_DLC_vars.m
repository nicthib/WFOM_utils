function DLC = extract_DLC_vars(file_path,var_name)
% When you point to a folder, this loads each variable matching the
% filename variable, assigns the id using the filename's first word (up to
% an underscore), then appends it to the DLC variable.
DLC = [];
files = dir(fullfile(file_path,'*.mat'));
for i = 1:numel(files)
    eval(['clear ' var_name])
    load(fullfile(file_path,files(i).name),var_name)
    id = files(i).name(1:(table(strfind(files(i).name,'_')).Var1(1)-1));
    DLC.(id) = [];
    if exist(var_name,'var')
        fn = fieldnames(eval(var_name));
        for j = 1:numel(fn)
            DLC.(id).(strrep(fn{j},'-','_')) = eval(var_name).(fn{j});
        end        
    end
end