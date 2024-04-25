%% Figure 3(a)

for k=1:count_dataset
    zz = params_5(1:no_organoids(k),2,k);
    v = log10(zz);
    y = (v-mean(v))/std(v);
    mean(v)
    std(v)
    figure;
    if k <= 4
        histogram(v,20,'Normalization','pdf')
    else
        histogram(v,15,'Normalization','pdf')
    end
    hold on
    fplot(@(x) normpdf(x,mean(v),std(v)),'LineWidth',2)
    set(gca,'xlim',[-1.5 1.5]);
    set(gca,'ylim',[0 3]);
    [h,p] = kstest(y)
end

%% Figure 3(b)

for ell=1:3
    figure;
    for k=2*ell-1:2*ell
        b = 1:no_organoids(k); 
        scatter(params_5(b,2,k),params_5(b,3,k),40);
        set(gca,'Xscale','log')
        set(gca,'Yscale','log')
        xlabel('Initial exponential growth rate (a)')
        ylabel('Growth rate decay parameter (b)')
        a = find(r(1:no_organoids(k),k)==0);

        hold on
     
        contours = [1 2 3 4];
        for m=1:size(contours,2)
            fplot(@(t) t*log10(exp(1))/contours(m));
        end
        set(gca,'xlim',[0.05,10])
        set(gca,'ylim',[10^(-4),10])
    
        legend('Organoids','K=10','K=10^2','K=10^3','K=10^4')
    end
    if ell==3
        for k=2*ell+1
           b = 1:no_organoids(k); 
            scatter(params_5(b,2,k),params_5(b,3,k),40);
            set(gca,'Xscale','log')
            set(gca,'Yscale','log')
            xlabel('Initial exponential growth rate (a)')
            ylabel('Growth rate decay parameter (b)')
    
            hold on
         
            contours = [1 2 3 4];
            for m=1:size(contours,2)
                fplot(@(t) t*log10(exp(1))/contours(m));
            end
            set(gca,'xlim',[0.05,10])
            set(gca,'ylim',[10^(-4),10])
        
            legend('Organoids','K=10','K=10^2','K=10^3','K=10^4')
        end
    end
    if ell ~= 3
        z = [params_5(find(r(1:no_organoids(2*ell-1),2*ell-1)==0),2,2*ell-1);params_5(find(r(1:no_organoids(2*ell),2*ell)==0),2,2*ell)];
        w = [params_5(find(r(1:no_organoids(2*ell-1),2*ell-1)==0),3,2*ell-1);params_5(find(r(1:no_organoids(2*ell),2*ell)==0),3,2*ell)];
    else
        z = [params_5(find(r(1:no_organoids(2*ell-1),2*ell-1)==0),2,2*ell-1);params_5(find(r(1:no_organoids(2*ell),2*ell)==0),2,2*ell);params_5(find(r(1:no_organoids(2*ell+1),2*ell+1)==0),2,2*ell+1)];
        w = [params_5(find(r(1:no_organoids(2*ell-1),2*ell-1)==0),3,2*ell-1);params_5(find(r(1:no_organoids(2*ell),2*ell)==0),3,2*ell);params_5(find(r(1:no_organoids(2*ell+1),2*ell+1)==0),3,2*ell+1)];
    end
    corr(z,w)
end

%% Figure 5(a)

figure;
for k=1:3
if k ~= 3
    z = [params_5(1:no_organoids(2*k-1),2,2*k-1);params_5(1:no_organoids(2*k),2,2*k)];
    v = log10(z);
else
    z = [params_5(1:no_organoids(2*k-1),2,2*k-1);params_5(1:no_organoids(2*k),2,2*k);params_5(1:no_organoids(2*k+1),2,2*k+1)];
    v = log10(z);
end
y = (v-mean(v))/std(v);
kstest(y)
hold on
fplot(@(x) normpdf(x,mean(v),std(v)),'LineWidth',2)
end
legend('UK','UP','US')
set(gca,'xlim',[-1 1]);

%% Figure 5(b)

figure;
perc_gomp_class = zeros(3,5);
for ell=1:3
    count_first = 0;
    count_second = 0;
    count_third = 0;
    count_fourth = 0;
    count_fifth = 0;
    for k=2*ell-1:2*ell
        for u=1:no_organoids(k)
            if r(u,k) == 0
                if exp(params_5(u,2,k)/params_5(u,3,k)) < 10^1
                    count_first = count_first + 1;
                elseif exp(params_5(u,2,k)/params_5(u,3,k)) < 10^2
                    count_second = count_second + 1;
                elseif exp(params_5(u,2,k)/params_5(u,3,k)) < 10^3
                    count_third = count_third + 1;
                elseif exp(params_5(u,2,k)/params_5(u,3,k)) < 10^4
                    count_fourth = count_fourth + 1;
                else
                    count_fifth = count_fifth + 1;
                end
            end
        end
    end
    if ell==3
        for k=2*ell+1
            for u=1:no_organoids(k)
                if r(u,k) == 0
                    if exp(params_5(u,2,k)/params_5(u,3,k)) < 10^1
                        count_first = count_first + 1;
                    elseif exp(params_5(u,2,k)/params_5(u,3,k)) < 10^2
                        count_second = count_second + 1;
                    elseif exp(params_5(u,2,k)/params_5(u,3,k)) < 10^3
                        count_third = count_third + 1;
                    elseif exp(params_5(u,2,k)/params_5(u,3,k)) < 10^4
                        count_fourth = count_fourth + 1;
                    else
                        count_fifth = count_fifth + 1;
                    end
                end
            end
        end
    end
   perc_gomp_class(ell,:) = [count_first,count_second,count_third,count_fourth,count_fifth]/(count_first+count_second+count_third+count_fourth+count_fifth);
end
bar(perc_gomp_class,'stacked')
set(gca,'ylim',[0 1.1])
legend('UK','UP','US');

%% Figure 5(c)

figure;
for k=1:3
if k ~= 3
    z = [x_0_all(1:no_organoids(2*k-1),2,2*k-1);x_0_all(1:no_organoids(2*k),2,2*k)];
else
    z = [x_0_all(1:no_organoids(2*k-1),2,2*k-1);x_0_all(1:no_organoids(2*k),2,2*k);x_0_all(1:no_organoids(2*k+1),2,2*k+1)];
end
v = log10(z);
hold on
[f,xi] = ksdensity(v,-0.5:0.001:3.5,'Function','cdf');
plot(xi,f);
end
set(gca,'xlim',[-0.5 3.5]);
legend('UK','UP','US')
set(gca,'ylim',[0 1.05]);
xlabel("log_{10}(organoid size on Day 0)");

%% Figure 5(d)

figure;
for k=1:3
if k ~= 3
    z = [params_5(1:no_organoids(2*k-1),2,2*k-1).*exp(-params_5(1:no_organoids(2*k-1),3,2*k-1).*(params_5(1:no_organoids(2*k-1),1,2*k-1)));params_5(1:no_organoids(2*k),2,2*k).*exp(-params_5(1:no_organoids(2*k),3,2*k).*(params_5(1:no_organoids(2*k),1,2*k)))];
else
    z = [params_5(1:no_organoids(2*k-1),2,2*k-1).*exp(-params_5(1:no_organoids(2*k-1),3,2*k-1).*(params_5(1:no_organoids(2*k-1),1,2*k-1)));params_5(1:no_organoids(2*k),2,2*k).*exp(-params_5(1:no_organoids(2*k),3,2*k).*(params_5(1:no_organoids(2*k),1,2*k)));params_5(1:no_organoids(2*k+1),2,2*k+1).*exp(-params_5(1:no_organoids(2*k+1),3,2*k+1).*(params_5(1:no_organoids(2*k+1),1,2*k+1)))];
end
v = log10(z);
hold on
[f,xi] = ksdensity(v,-1.5:0.001:0.5,'Function','cdf');
plot(xi,f);
end
set(gca,'xlim',[-1.5 0.5]);
set(gca,'ylim',[0 1.05]);
legend('UK','UP','US')
xlabel("log_{10}(growth rate on Day 0)");

%% Figure 5(e)

figure;
for k=1:3
if k ~= 3
    z = [x_0_all(1:no_organoids(2*k-1),4,2*k-1);x_0_all(1:no_organoids(2*k),4,2*k)];
else
    z = [x_0_all(1:no_organoids(2*k-1),4,2*k-1);x_0_all(1:no_organoids(2*k),4,2*k);x_0_all(1:no_organoids(2*k+1),4,2*k+1)];
end
v = log10(z);
hold on
[f,xi] = ksdensity(v,-0.5:0.001:3.5,'Function','cdf');
plot(xi,f);
end
legend('UK','UP','US')
set(gca,'xlim',[-0.5 3.5]);
xlabel("log_{10}(organoid size on Day 5)");
set(gca,'ylim',[0 1.05]);

%% Figure 5(f)

figure;
for k=1:3
if k ~= 3
    z = [params_5(1:no_organoids(2*k-1),2,2*k-1).*exp(-params_5(1:no_organoids(2*k-1),3,2*k-1).*(params_5(1:no_organoids(2*k-1),1,2*k-1)+5));params_5(1:no_organoids(2*k),2,2*k).*exp(-params_5(1:no_organoids(2*k),3,2*k).*(params_5(1:no_organoids(2*k),1,2*k)+5))];
else
    z = [params_5(1:no_organoids(2*k-1),2,2*k-1).*exp(-params_5(1:no_organoids(2*k-1),3,2*k-1).*(params_5(1:no_organoids(2*k-1),1,2*k-1)+5));params_5(1:no_organoids(2*k),2,2*k).*exp(-params_5(1:no_organoids(2*k),3,2*k).*(params_5(1:no_organoids(2*k),1,2*k)+5));params_5(1:no_organoids(2*k+1),2,2*k+1).*exp(-params_5(1:no_organoids(2*k+1),3,2*k+1).*(params_5(1:no_organoids(2*k+1),1,2*k+1)+5))];
end
v = log10(z);
hold on
[f,xi] = ksdensity(v,-1.5:0.001:0.5,'Function','cdf');
plot(xi,f);
end
legend('UK','UP','US')
set(gca,'xlim',[-1.5 0.5]);
xlabel("log_{10}(growth rate on Day 5)");
set(gca,'ylim',[0 1.05]);

%% Figures 5(g) and 5(i)

iter = 100000;

params_uk = [params_5(1:no_organoids(1),:,1);params_5(1:no_organoids(2),:,2)];
params_up = [params_5(1:no_organoids(3),:,3);params_5(1:no_organoids(4),:,4)];
params_us = [params_5(1:no_organoids(5),:,5);params_5(1:no_organoids(6),:,6);params_5(1:no_organoids(7),:,7)];

no_organoids_uk = no_organoids(1)+no_organoids(2);
no_organoids_up = no_organoids(3)+no_organoids(4);
no_organoids_us = no_organoids(5)+no_organoids(6)+no_organoids(7);

F = @(x,t) exp(x(1)./x(2).*(1-exp(-x(2).*t)));

final_day = 24;
vec_range = 0:0.01:final_day;

fcn_vec_uk = zeros(iter,size(vec_range,2));
for i = 1:iter
    x = params_uk(randi(no_organoids_uk),:);
    g = @(c) (c <= 4-x(1)) + (c>4-x(1)).*F([x(2) x(3)],c-4+x(1));
    fcn_vec_uk(i,:) = g(vec_range);
end
fcn_vec_median_uk = median(fcn_vec_uk);

fcn_vec_up = zeros(iter,size(vec_range,2));
for i = 1:iter
    x = params_up(randi(no_organoids_up),:);
    g = @(c) (c <= 4-x(1)) + (c>4-x(1)).*F([x(2) x(3)],c-4+x(1));
    fcn_vec_up(i,:) = g(vec_range);
end
fcn_vec_median_up = median(fcn_vec_up);

fcn_vec_us = zeros(iter,size(vec_range,2));
for i = 1:iter
    x = params_us(randi(no_organoids_us),:);
    g = @(c) (c <= 4-x(1)) + (c>4-x(1)).*F([x(2) x(3)],c-4+x(1));
    fcn_vec_us(i,:) = g(vec_range);
end
fcn_vec_median_us = median(fcn_vec_us);

figure;
hold on
set(gca,'Yscale','log')
set(gca,'ylim',[0 10^4])
plot(vec_range-4,fcn_vec_median_uk)
plot(vec_range-4,fcn_vec_median_up)
plot(vec_range-4,fcn_vec_median_us)
set(gca,'Yscale','log')
set(gca,'ylim',[0 10^4])

fcn_vec_median_uk_cum = zeros(1,size(vec_range,2));
fcn_vec_median_up_cum = zeros(1,size(vec_range,2)); 
fcn_vec_median_us_cum = zeros(1,size(vec_range,2));
for i=1:size(vec_range,2)
    fcn_vec_median_uk_cum(i) = sum(fcn_vec_median_uk(1:i))*(vec_range(2)-vec_range(1));
    fcn_vec_median_up_cum(i) = sum(fcn_vec_median_up(1:i))*(vec_range(2)-vec_range(1));
    fcn_vec_median_us_cum(i) = sum(fcn_vec_median_us(1:i))*(vec_range(2)-vec_range(1));
end

figure;
mut_rate = 10^(-4);
plot(vec_range-4,mut_rate.*fcn_vec_median_uk_cum)
hold on
plot(vec_range-4,mut_rate.*fcn_vec_median_up_cum)
plot(vec_range-4,mut_rate.*fcn_vec_median_us_cum)
set(gca,'Yscale','log')
set(gca,'xlim',[0 20])


%% Figure 5(h)

iter = 100000;

params_uk = [params_5(1:no_organoids(1),:,1);params_5(1:no_organoids(2),:,2)];
params_up = [params_5(1:no_organoids(3),:,3);params_5(1:no_organoids(4),:,4)];
params_us = [params_5(1:no_organoids(5),:,5);params_5(1:no_organoids(6),:,6);params_5(1:no_organoids(7),:,7)];

no_organoids_uk = no_organoids(1)+no_organoids(2);
no_organoids_up = no_organoids(3)+no_organoids(4);
no_organoids_us = no_organoids(5)+no_organoids(6)+no_organoids(7);

F = @(x,t) exp(x(1)./x(2).*(1-exp(-x(2).*t)));

final_day = 24;
vec_range = 0:1:final_day;

fcn_vec_uk = zeros(iter,size(vec_range,2));
for i = 1:iter
    x = params_uk(randi(no_organoids_uk),:);
    g = @(c) (c <= 4-x(1)) + (c>4-x(1)).*F([x(2) x(3)],c-4+x(1));
    fcn_vec_uk(i,:) = g(vec_range);
end
fcn_vec_median_uk = median(fcn_vec_uk);

fcn_vec_up = zeros(iter,size(vec_range,2));
for i = 1:iter
    x = params_up(randi(no_organoids_up),:);
    g = @(c) (c <= 4-x(1)) + (c>4-x(1)).*F([x(2) x(3)],c-4+x(1));
    fcn_vec_up(i,:) = g(vec_range);
end
fcn_vec_median_up = median(fcn_vec_up);

fcn_vec_us = zeros(iter,size(vec_range,2));
for i = 1:iter
    x = params_us(randi(no_organoids_us),:);
    g = @(c) (c <= 4-x(1)) + (c>4-x(1)).*F([x(2) x(3)],c-4+x(1));
    fcn_vec_us(i,:) = g(vec_range);
end
fcn_vec_median_us = median(fcn_vec_us);

figure;
hold on
errorbar(vec_range-4,fcn_vec_median_uk,fcn_vec_median_uk-prctile(fcn_vec_uk,5),prctile(fcn_vec_uk,95)-fcn_vec_median_uk)
set(gca,'Yscale','log')
set(gca,'ylim',[0 10^8])
