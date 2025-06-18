close all
n = 1000;
sampleMax = 1000;
iterbMax = 7;
iteraMax = 4;
m = 840;
factm = m./(8:-1:2);
err = zeros(iteraMax,iterbMax,sampleMax);
for sample = 1:sampleMax
    rng(sample)
    sample
    X0 = randn(n,m/2);
    X0 = orth(X0);
for itera = 1:iteraMax
    for iterb = 1:iterbMax

        alpha = 2^(1-itera);
        aa = [(0:m-1)*alpha,m+1:n]';
        b = factm(iterb);
        A = @(X) X.*aa;
        X = X0(:,1:b);
        Q = Blanczos(A,X,m);
        err(itera,iterb,sample) = norm(Q(m+1:n,:)/Q(1:m,:));

    end
end
end

stat = genStat(err);
xx = 1:7;
linw = 2;
figure
hold on
for itera = 1:iteraMax
fill([9-xx, fliplr(9-xx)], [stat(itera,xx,4), fliplr(stat(itera,xx,3))],[0.5 0.5 0.5], 'EdgeColor', 'none','HandleVisibility','off');
plot(9-xx,stat(itera,xx,1),'-','LineWidth',linw,'DisplayName',['$\alpha=',num2str(2^(1-itera)),'$'])
end
xlabel('$d=m/b$','FontSize',18,'Interpreter','latex');
hold off
set(gca,'yscale','log')
legend('FontSize',18,'Interpreter','latex','Location','northwest','Box','off')
set(gcf, 'Color', 'w');

export_fig('fig/EPa.pdf')
export_fig('fig/EPa.eps')



function stat = genStat(data)

[m1,m2,n] = size(data);
stat = zeros(m1,m2,4);
alpha = 0.05;
for ii = 1:m1
    for jj = 1:m2
    outputs = data(ii,jj,:);
    mean_value = mean(outputs);
    std_dev = std(outputs);
    t_value = tinv(1 - alpha/2, n - 1);
    sem = std_dev / sqrt(n);
    stat(ii,jj,1) = mean_value;
    stat(ii,jj,2) = std_dev;
    stat(ii,jj,3) = mean_value - t_value * sem;
    stat(ii,jj,4) = mean_value + t_value * sem;

    end
end
end