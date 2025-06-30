close all
n = 1000;
sampleMax = 1000;
iterdMax = 4;
iterbMax = 12;
m = 60;
err = zeros(iterdMax,iterbMax,sampleMax);
parfor sample = 1:sampleMax
    rng(sample)
    sample
    X0 = randn(n,m);
    X0 = orth(X0);
    for iterd = 1:iterdMax
        d = iterd+1;
        b = m/d;
        for iterbeta = 1:iterbMax

            beta = 2^(-iterbeta);
            aa = [];
            for ii = 1:d
                alpha = ii*(1+2*beta)/m;
                aa = [aa,linspace(alpha-beta,alpha+beta,b)];
            end
            aa = [aa,(-n+m:-1)/(n-m)-1]';
            A = @(X) X.*aa;
            X = X0(:,1:b);
            Q = Blanczos(A,X,m);
            err(iterd,iterbeta,sample) = norm(Q(m+1:n,:)/Q(1:m,:));

        end
    end
end

stat = genStat(err);

xx = 2.^(-1:-1:-iterbMax);
linw = 2;
figure
hold on
color = {'r','b','k','g'};
for iterd = 1:iterdMax
    fill([xx, fliplr(xx)], [stat(iterd,:,3), fliplr(stat(iterd,:,2))],[0.5 0.5 0.5], 'EdgeColor', 'none','HandleVisibility','off');
    plot(xx,stat(iterd,:,1),[color{iterd},'-'],'LineWidth',linw,'DisplayName',['$d=',num2str(iterd+1),'$'])
end
xlabel('$\beta$','FontSize',18,'Interpreter','latex');
% ylabel('$\tan\angle (Q,\mathcal{K}_{d}(A,\Omega))$','FontSize',18,'Interpreter','latex')
hold off
set(gca,'xscale','log')
set(gca,'yscale','log')
axis([2^(-12),2^-1,-inf,1e12])
legend('FontSize',18,'Interpreter','latex','Location','northeast','Box','off','Orientation', 'horizontal')
set(gcf, 'Color', 'w');

export_fig('fig/EPbeta.pdf')
export_fig('fig/EPbeta.eps')

save('databeta')

function stat = genStat(data)

[m1,m2,n] = size(data);
stat = zeros(m1,m2,3);
for ii = 1:m1
    for jj = 1:m2
        outputs = sort(data(ii,jj,:));
        idx = round([0.5,0.25,0.75]*n);
        stat(ii,jj,:) = outputs(idx);

    end
end
end