function sp_inds = sp_extent(idx,sp_sz,fig_sz)
v = reshape(1:prod(fig_sz),fliplr(fig_sz))';
sp_inds = v(idx(1):idx(1)+sp_sz(1)-1,idx(2):idx(2)+sp_sz(2)-1);
sp_inds = [min(sp_inds) max(sp_inds)];