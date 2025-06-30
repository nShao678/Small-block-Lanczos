close all
n = 1000;
sampleMax = 1000;
iterdMax = 3;
iteraMax = 10;
m = 60;

beta = 1e-4;
err = zeros(iterdMax,iteraMax,sampleMax);
parfor sample = 1:sampleMax
    rng(sample)
    sample
    X0 = randn(n,m);
    X0 = orth(X0);
    for iterd = 1:iterdMax
        d = iterd+1;
        b = m/d;
        for iteralpha = 1:iteraMax
            
            aa = [];
            for ii = 1:d
                alpha = ii*(2^(-iteralpha)+2*beta)/m;
                aa = [aa,linspace(alpha-beta,alpha+beta,b)];
            end
            aa = [aa,(-n+m:-1)/(n-m)-1]';
            A = @(X) X.*aa;
            X = X0(:,1:b);
            Q = Blanczos(A,X,m);
            err(iterd,iteralpha,sample) = norm(Q(m+1:n,:)/Q(1:m,:));

        end
    end
end


stat = genStat(err);


xx = 2.^(-1:-1:-iteraMax);
linw = 2;
figure
hold on
color = {'r','b','k','g'};
for iterd = 1:iterdMax
    fill([xx, fliplr(xx)], [stat(iterd,:,3), fliplr(stat(iterd,:,2))],[0.5 0.5 0.5], 'EdgeColor', 'none','HandleVisibility','off');
    plot(xx,stat(iterd,:,1),[color{iterd},'-'],'LineWidth',linw,'DisplayName',['$d=',num2str(iterd+1),'$'])
    plot(xx,xx.^(-iterd)*stat(iterd,5,1)/5/xx(5).^(-iterd),[color{iterd},'--'],'linewidth',linw,'DisplayName',['$y\propto\alpha^{',num2str(-iterd),'}$'])
end
xlabel('$\alpha$','FontSize',18,'Interpreter','latex');
hold off
set(gca,'xscale','log')
set(gca,'yscale','log')
legend('FontSize',18,'NumColumns',2,'Interpreter','latex','Location','northeast','Box','off','Orientation', 'horizontal')
axis([2^(-10),2^-1,1e2,1e20])
legend('FontSize',18,'Interpreter','latex','Location','northeast','Box','off')
set(gcf, 'Color', 'w');

export_fig('fig/EPgap.pdf')
export_fig('fig/EPgap.eps')

save('datagap')

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