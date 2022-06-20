% General putpose script that checks files for a bad aux field in metadata. Sometimes, aux fields contain cells instead of matrices. This fixes them into matrixes.
mice = {'cm124','cm125','cm126','cm127','cm128'};
day = '1234567';
for a = 4:numel(mice)
    for b = 1:numel(day)
        auxpath = dir(['/local_mount/space/dingus/1/cmdata/' mice{a} '_' day(b) '/auxillary/**/*_b.mat']);
        if ~isempty(auxpath)
            for c = 1:numel(auxpath)
                fn = fullfile(auxpath(c).folder,auxpath(c).name);
                m = matfile(fn);
                if iscell(m.aux)
                    load(fn)
                    tmpaux = [];
                    disp('aux is bad!!')
                    for d = 1:numel(m.aux)
                        if numel(aux{d}) > 1
                            tmpaux(:,d) = aux{d}(1:2);
                        else
                            tmpaux(:,d) = 0;
                        end
                    end
                    aux = tmpaux;
                    save(fn,'aux')
                end
            end
        end
    end
end

