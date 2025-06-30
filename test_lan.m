close all
rng(1)
n = 2000;
iterbMax = 6;
iteraMax = 4;
m = 2^(iterbMax-1);
tol = 1e-10;
err = ones(iteraMax,iterbMax,n);

tab = zeros(iteraMax,iterbMax);
X0 = randn(n,m);
X0 = orth(X0);

for iterb = 1:iterbMax
    b = 2^(iterb-1);
    d = m/b;
    for iteralpha = 1:iteraMax    
        a = 1+linspace(10^(1-iteralpha),0,m)';
        aa = [a;(-n+m:-1)'/(n-m)];
        A = @(X) X.*aa;
        X = X0(:,1:b);
        [Q,T] = Blanczos(A,X,1024);
        for iter = m/b:1024/b
            idx = 1:iter*b;
            eval = sort(eig(T(idx,idx)),'descend');
            err(iteralpha,iterb,iter) = sum(abs(eval(1:m)-a));
%             err(iteralpha,iterb,iter)
            if err(iteralpha,iterb,iter)<tol
                tab(iteralpha,iterb) = iter*b;
                break
            end
            
        end

    end
end
tabr = tab./(tab(:,1))-1;
save('introLan')

% 
% xx = 2.^(-1:-1:-iteraMax);
% linw = 2;
% figure
% hold on
% color = {'r','b','k','g'};
% for iterd = 1:iterdMax
%     fill([xx, fliplr(xx)], [stat(iterd,:,3), fliplr(stat(iterd,:,2))],[0.5 0.5 0.5], 'EdgeColor', 'none','HandleVisibility','off');
%     plot(xx,stat(iterd,:,1),[color{iterd},'-'],'LineWidth',linw,'DisplayName',['$d=',num2str(iterd+1),'$'])
%     plot(xx,xx.^(-iterd)*stat(iterd,1,1)/1.1/xx(1).^(-iterd-1),[color{iterd},'--'],'linewidth',linw,'DisplayName',['$y\propto\alpha^{',num2str(-iterd),'}$'])
% end
% xlabel('$\alpha$','FontSize',18,'Interpreter','latex');
% hold off
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% legend('FontSize',18,'NumColumns',2,'Interpreter','latex','Location','northeast','Box','off','Orientation', 'horizontal')
% axis([2^(-10),2^-1,-inf,inf])
% legend('FontSize',18,'Interpreter','latex','Location','northeast','Box','off')
% set(gcf, 'Color', 'w');
% 
% export_fig('fig/EPgap.pdf')
% export_fig('fig/EPgap.eps')

