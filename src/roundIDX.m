function IDX = roundIDX(IDX)
sqsz = 2; IDX(isnan(IDX)) = 0;
for i = sqsz+1:size(IDX,1)-sqsz
    for j = sqsz+1:size(IDX,2)-sqsz
        if round(IDX(i,j)) == IDX(i,j)
        elseif isnan(IDX(i,j))
            IDX(i,j) = NaN;
        else
            tmp = IDX(i-sqsz:i+sqsz,j-sqsz:j+sqsz);
            IDX(i,j) = round(mode(tmp(:)));
        end
    end
end
IDX = round(IDX);