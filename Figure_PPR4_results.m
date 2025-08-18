function Figure_PPR4_results(model_best)
%%

figure('position',[100   100   1500   300]);
t = tiledlayout(2, 9);
t.TileSpacing = 'compact';
t.Padding = 'compact';

filters=model_best.filters;
res2=size(filters,1);
res=sqrt(res2);
nfilters=size(filters,2);
filters=reshape(filters,res,res,nfilters);
clim=max( abs(max(filters(:))),abs(min(filters(:))) );

for j=1:nfilters
    nexttile

    % set(gcf,'position',[100   100   1500   500]);
    imagesc(filters(:,:,j))
    colormap gray
    axis equal
    axis off

end

nexttile
axis off
axis([0 1 0 1])
y=1; dy=0.9/6;
y=y-dy; text(0,y,['cell id=' strrep(model_best.cell_id,'_',' ')])

cc=nanmean(model_best.ccs);
if cc<0.1
    y=y-dy; text(0,y,['Vali CC=' num2str(cc)],'Color','r')
else
    y=y-dy; text(0,y,['Vali CC=' num2str(cc)])
end

test_cc = model_best.test_cc;
if test_cc<0.1
    y=y-dy; text(0,y,['Test CC=' num2str(test_cc)],'Color','r')
else
    y=y-dy; text(0,y,['Test CC=' num2str(test_cc)])
end

y=y-dy; text(0,y,['Lamda=' num2str(model_best.lamda)])

nonlin=model_best.nonlin;
for j=1:nfilters
    nexttile(9+j)
    errorbar(nonlin(j).x, nonlin(j).y, nonlin(j).e)
end
