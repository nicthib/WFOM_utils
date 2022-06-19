
% t_mat = gettransitionmatrix(st,ord,show) calculates the transition matrix
% for the 1 x t dataset st.
% row --> column

function t_mat = gettransitionmatrix(st1,st2,k,show)
t_mat = zeros(k,k);
for i = 1:numel(st1)
    t_mat(st1(i),st2(i)) = t_mat(st1(i),st2(i)) + 1;
end
%t_mat = t_mat/sum(t_mat(:));
%t_mat = round(t_mat,3);
if show

    vis_mat = ones(size(t_mat));
    for i = 1:size(t_mat,1)
       vis_mat(i,i) = NaN;
    end

    im = imagesc(t_mat); axis image
    im.AlphaData = ~isnan(vis_mat);
    caxis([-6 -2]); colorbar; 
    for i = 1:k
        for j = 1:k
            if i ~= j
                text(i,j,mat2str(t_mat(i,j)),'HorizontalAlignment','center')
            end
        end
    end
    ylabel('Current state')
    xlabel('Next state')
    title(sprintf('Transition matrix'))
end