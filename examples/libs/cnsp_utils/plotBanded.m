function h=plotBanded(stats,names)
lambdas = unique(stats.lambdas(:,1));

if size(stats.lambdas,2)>2
    error('Plotting can only be done for two-dimensional banding')
else
for ii = 1:length(lambdas)
    for jj = 1:length(lambdas)
        I = stats.lambdas(:,1)==lambdas(ii) & stats.lambdas(:,2)==lambdas(jj);
        r_mat(ii,jj) = mean(mean(stats.r(:,I,:),1),3);
    end
end
end

[~,bestLambda] = max(squeeze(mean(mean(stats.r,1),3)));

Ibanded = [find(lambdas==stats.lambdas(bestLambda,1)') find(lambdas==stats.lambdas(bestLambda,2)')];
[~,I] = max(r_mat(logical(eye(size(r_mat)))));


h=figure;
subplot(3,1,1:2)
imagesc(r_mat')
axis xy
xticklabels(lambdas)
yticklabels(lambdas)
xlabel(['Lambda Feature ' names(1)])
ylabel(['Lambda Feature ' names(2)])
hold on
ax = axis(gca);
plot(ax(1:2),ax(1:2),'k')
plot(Ibanded(1),Ibanded(2),'k*')
plot(I,I,'ko')
colorbar

subplot(3,1,3)
hold on
plot(r_mat(:,Ibanded(2)))
plot(r_mat(Ibanded(1),:))
plot(r_mat(find(eye(size(r_mat)))))
xticklabels(lambdas)
xlabel('Lambda')
ylabel('pred.acc.')
legend({[names(1)' ' tuning curve'],[names(2) ' tuning curve'],'ridge tuning curve'},'Location','southwest')

