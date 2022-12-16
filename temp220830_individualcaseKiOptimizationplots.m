%open ki comparison values for given algorithm
netmeanerr_umb=squeeze(mean(netmeanerr_um,3));
k=1;
figure
for i=1:length(cases)
    for j=1:length(Dall) 
        subplot(length(Dall),length(cases),k)
        imagesc(squeeze(netmeanerr_umb(:,:,i,j))')
        stdaxisformat(ki_xy,ki_z)
%         title(['D = ' num2str(Dall(j)./1e12) ' s, b = ' num2str(cases(i,1)./1e3) ', ' num2str(cases(i,2)./1e3) ' kcps'])
        k=k+1;
    end
end

%% local functions must remain at bottom of file
function stdaxisformat(xtxt,ytext)
%generate axes labels
for j=1:length(xtxt)
    exylabs(j)={num2str(xtxt(j))};
end

for j=1:length(ytext)
    yislabs(j)={num2str(ytext(j))};
end
%put axis directions back to normal
set(gca,'YDir','normal')

%the formatting of the plots
xlabel('ki xy')
xticks(1:length(xtxt))
xticklabels(exylabs)
xtickangle(90)
ylabel('ki z')
yticks(1:length(ytext))
yticklabels(yislabs)

% c=colorbar;
% c.Label.String = 'stage position error (um)';
caxis([0,5])
% axis square

fh = findall(gcf,'Type','Figure');
txt_obj = findall(fh,'Type','text');
set(txt_obj,'FontName','Arial','FontSize',10);

ah = findall(gcf,'Type','Axes');
set(ah,'FontSize',6);
set(ah,'FontName','Arial');
end