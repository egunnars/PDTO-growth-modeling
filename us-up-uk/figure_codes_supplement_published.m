%% Figure S3

for k=1:count_outer
    figure;
    hold on
    scatter(error_vec_5_norm(index_set,k), error_vec_6_norm(index_set,k),50)
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    xlabel('Fitting error for Gompertz model')
    ylabel('Fitting error for logistic model')
    fplot(@(t) t)
    set(gca,'xlim',[10^(-8) 1])
    set(gca,'ylim',[10^(-8) 1])
end

for k=1:count_outer
    figure;
    hold on
    scatter(error_vec_5_norm(index_set,k), error_vec_9_norm(index_set,k),50)
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    xlabel('Fitting error for Gompertz model')
    ylabel('Fitting error for von Bertalanffy model')
    fplot(@(t) t)
    set(gca,'xlim',[10^(-8) 1])
    set(gca,'ylim',[10^(-8) 1])
end

%% Figure S4

for k=1:3
    figure;
        if k ~= 3
            z = [params_5(1:no_organoids(2*k-1),2,2*k-1);params_5(1:no_organoids(2*k),2,2*k)];
            v = log10(z(find(z>0)));
        else
            z = [params_5(1:no_organoids(2*k-1),2,2*k-1);params_5(1:no_organoids(2*k),2,2*k);params_5(1:no_organoids(2*k+1),2,2*k+1)];
            v = log10(z(find(z>0)));
        end
    mean(v)
    std(v)
    y = (v-mean(v))/std(v);
    hold on
    histogram(v,20,'Normalization','pdf')
    fplot(@(x) normpdf(x,mean(v),std(v)),'LineWidth',2)
    set(gca,'xlim',[-1.5 1.5]);
    set(gca,'ylim',[0 2.5]);
    [h,p] = kstest(y)
end

%% Figure S6
    
for ell=1:3
    figure;
    for k = 2*ell-1:2*ell
        b = find(r(1:no_organoids(k),k)==0);
        scatter(x_0_all(b,4,k),min(exp(params_5(b,2,k)./params_5(b,3,k)),10^10),50);
        set(gca,'Xscale','log')
        set(gca,'Yscale','log')
    xlabel('Observed volume on Day 5')
    ylabel('Carrying capacity')

        hold on
    end
    if ell==3
        for k = 2*ell+1
        b = 1:no_organoids(k);
        scatter(x_0_all(b,4,k),min(exp(params_5(b,2,k)./params_5(b,3,k)),10^10),50);
        set(gca,'Xscale','log')
        set(gca,'Yscale','log')
        xlabel('Observed volume on Day 5')
        ylabel('Carrying capacity')
        end
    end

    yline(10^4);
    fplot(@(t) t);
    set(gca,'xlim',[10^0 10^4]);
    set(gca,'ylim',[10^0 10^10]);
end

%% Figure S7

indices = {1:12;13:24;25:36;37:48;49:54;55:66;67:78};
n_exp = size(indices,1);

center_x = 1600;
center_y = 1600;

summary = zeros(1,3,n_exp);

ell = 1;
for count_outer = 1:n_exp
    count = 0;
    for l = cell2mat(indices(count_outer,:))
        for m = min(data_tracking{l}.TrackID):max(data_tracking{l}.TrackID)
            if ~isempty(find(data_tracking{l}.TrackID == m & data_tracking{l}.Index_t_ == ell))
                a = find(data_tracking{l}.TrackID == m & data_tracking{l}.Index_t_ == ell);
                if data_tracking{l}.Gomp2(a(1)) > 0
                    count = count+1;
                    summary(count,1,count_outer) = data_tracking{l}.Gomp1(a(1));
                    summary(count,2,count_outer) = exp(data_tracking{l}.Gomp1(a(1))/data_tracking{l}.Gomp2(a(1)));
                    summary(count,3,count_outer) = sqrt((data_tracking{l}.CenterOfMassX__m_(a(1))-center_x)^2+(data_tracking{l}.CenterOfMassY__m_(a(1))-center_y)^2)/1000;
                end
            end
        end
    end
end

for ell = 1:3
    X = [];
    y = [];
    for k = 2*ell-1:2*ell
        X = [X;[summary(1:no_organoids(k),3,k),ones(no_organoids(k),1)]];
        y = [y;summary(1:no_organoids(k),1,k)];
    end
    if ell==3
        for k=2*ell+1
            X = [X;[summary(1:no_organoids(k),3,k),ones(no_organoids(k),1)]];
            y = [y;summary(1:no_organoids(k),1,k)];
        end
    end
    figure;
    scatter(X(:,1),y,40);
    mdl = fitlm(X,y)
    x = X\y;
    hold on
    fplot(@(t) x(1)*t+x(2));
    set(gca,'Xlim',[0 1.6])
    set(gca,'Ylim',[0 6]);
end

for ell = 1:3
    X = [];
    y = [];
    for k = 2*ell-1:2*ell
        b = find(r(1:no_organoids(k),k)==0);
        X = [X;[summary(b,3,k),ones(size(b,1),1)]];
        y = [y;log10(min(summary(b,2,k),10^10))];
    end
    if ell==3
        for k=2*ell+1
            b = find(r(1:no_organoids(k),k)==0);
            X = [X;[summary(b,3,k),ones(size(b,1),1)]];
            y = [y;log10(min(summary(b,2,k),10^10))];
        end
    end
    figure;
    scatter(X(:,1),y,40);
    mdl = fitlm(X,y)
    x = X\y;
    hold on
    fplot(@(t) x(1)*t+x(2));
    set(gca,'Xlim',[0 1.6])
    set(gca,'Ylim',[0 10]);
end

%% Figure S8

perc_exp = zeros(3,1);
figure;
for k = 1:3
    if k ~= 3
        z = [r(1:no_organoids(2*k-1),2*k-1);r(1:no_organoids(2*k),2*k)];
        w = [params_5(1:no_organoids(2*k-1),2,2*k-1);params_5(1:no_organoids(2*k),2,2*k)];
        perc_exp(k) = sum(z)/size(z,1);
        w = w(find(z==1));
        w = log10(w);
        [f,xi] = ksdensity(w,-1.5:0.001:1,'Function','cdf');
        plot(xi,f);
        hold on
    else
        z = [r(1:no_organoids(2*k-1),2*k-1);r(1:no_organoids(2*k),2*k);r(1:no_organoids(2*k+1),2*k+1)];
        w = [params_5(1:no_organoids(2*k-1),2,2*k-1);params_5(1:no_organoids(2*k),2,2*k);params_5(1:no_organoids(2*k+1),2,2*k+1)];
        perc_exp(k) = sum(z)/size(z,1);
        w = w(find(z==1));
        w = log10(w);
        [f,xi] = ksdensity(w,-1.5:0.001:1,'Function','cdf');
        plot(xi,f);
    end
end
legend('UK','UP','US');
set(gca,'xlim',[-1.2 0.6])
set(gca,'ylim',[0 1]);

figure;
for k = 1:3
    if k ~= 3
        z = [r(1:no_organoids(2*k-1),2*k-1);r(1:no_organoids(2*k),2*k)];
        y = [params_5(1:no_organoids(2*k-1),2,2*k-1);params_5(1:no_organoids(2*k),2,2*k)];
        w = [params_5(1:no_organoids(2*k-1),3,2*k-1);params_5(1:no_organoids(2*k),3,2*k)];
        y = y(find(z==0));
        w = w(find(z==0));
        v = log10(exp(y./w));
        [f,xi] = ksdensity(v,0:0.001:20,'Function','cdf');
        plot(xi,f);
        hold on
    else
        z = [r(1:no_organoids(2*k-1),2*k-1);r(1:no_organoids(2*k),2*k);r(1:no_organoids(2*k+1),2*k+1)];
        y = [params_5(1:no_organoids(2*k-1),2,2*k-1);params_5(1:no_organoids(2*k),2,2*k);params_5(1:no_organoids(2*k+1),2,2*k+1)];
        w = [params_5(1:no_organoids(2*k-1),3,2*k-1);params_5(1:no_organoids(2*k),3,2*k);params_5(1:no_organoids(2*k+1),3,2*k+1)];
        y = y(find(z==0));
        w = w(find(z==0));
        v = log10(exp(y./w));
        [f,xi] = ksdensity(v,0:0.001:20,'Function','cdf');
        plot(xi,f);
    end
end
legend('UK','UP','US');
set(gca,'ylim',[0 1]);

%% Figure S9

for k=1:3
    figure;
        if k ~= 3
            z = [params_5(1:no_organoids(2*k-1),1,2*k-1);params_5(1:no_organoids(2*k),1,2*k)];
        else
            z = [params_5(1:no_organoids(2*k-1),1,2*k-1);params_5(1:no_organoids(2*k),1,2*k);params_5(1:no_organoids(2*k+1),1,2*k+1)];
        end
    histogram(z,20,'Normalization','probability')
    set(gca,'ylim',[0 0.3])
end