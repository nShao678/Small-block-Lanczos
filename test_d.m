dMax = 20;
sampleMax = 1000;
b = 8;
dataNormS = zeros(sampleMax,dMax);

dataKappaOmega = zeros(1,dMax);
for d = 2:dMax
    Omega = cell(1,d);
    m = b*d;
    % A = reshape(randn(1,m),b,d);
    % A = A/max(max(abs(A)));
    A = reshape(1:m,b,d)/m;
    for sample = 1:sampleMax
        rng(sample)
        for ii = 1:d
            Omega{ii} = randn(b,b);
        end

        B = cell(1,d);
        for ii = 1:d
            B{ii} = Omega{ii}\(Omega{ii}.*(A(:,ii)));
        end

        kappaOmega = 1;
        normS = 1;
        for ii = 1:d
            [S,Omegah] = recur(A,B,Omega,ii);
            normS = max(normS,norm(inv(S)));
            for jj = 1:d
                if jj~=ii
                    kappaOmega = max(kappaOmega,cond(Omegah{jj}));
                end
            end
        end
        dataNormS(sample,d) = normS^(1/(d-1))/m;
        dataKappaOmega(sample,d) = kappaOmega^(1/(d-1));
    end
end

statKappaOmega = genStat(dataKappaOmega);
statNormS = genStat(dataNormS);

xx = 2:dMax;
linw = 2;
figure
hold on
yyaxis left
fill([xx, fliplr(xx)], [statKappaOmega(4,xx), fliplr(statKappaOmega(3,xx))],[0.5 0.5 0.5], 'EdgeColor', 'none','HandleVisibility','off');
plot(xx,statKappaOmega(1,xx),'b-','LineWidth',linw,'DisplayName','$\chi_{\Omega}$')
set(gca,'yscale','log')
yyaxis right
fill([xx, fliplr(xx)], [statNormS(4,xx), fliplr(statNormS(3,xx))],[0.5 0.5 0.5], 'EdgeColor', 'none','HandleVisibility','off');
plot(xx,statNormS(1,xx),'r-','LineWidth',linw,'DisplayName','$\chi_{S}$')
set(gca,'yscale','log')
xlabel('$d$','FontSize',18,'Interpreter','latex');
hold off
legend('FontSize',18,'Interpreter','latex','Box','off')
set(gcf, 'Color', 'w');
export_fig('fig/EPd.pdf')
export_fig('fig/EPd.eps')


function stat = genStat(data)

[n,m] = size(data);
stat = zeros(4,m);
alpha = 0.05;
for ii = 1:m
    outputs = data(:,ii);
    mean_value = mean(outputs);
    std_dev = std(outputs);
    t_value = tinv(1 - alpha/2, n - 1);
    sem = std_dev / sqrt(n);
    stat(1,ii) = mean_value;
    stat(2,ii) = std_dev;
    stat(3,ii) = mean_value - t_value * sem;
    stat(4,ii) = mean_value + t_value * sem;


end
end