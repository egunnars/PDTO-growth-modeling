tic

%% Definition of main parameters

%Fix the random seed for reproducibility.
rndseed = 1; rng(rndseed);

%Number of timepoints including Day -tau.
n_timepoints = 5;
%Number of free datapoints (initial size of each organoid is fixed).
n_datapoints = 4;
%Number of optimization runs in the least squares estimation.
iter = 1000;

%Define which wells should be grouped together.
indices = {2:2:12;50:2:60};
n_exp = size(indices,1);

%Number of parameters in the 9 models we consider,
%including the tau (starting time) parameter.
n_param = 2;
n_param_2 = 2;
n_param_3 = 2;
n_param_4 = 2;
n_param_5 = 3;
n_param_6 = 3;
n_param_7 = 3;
n_param_8 = 3;
n_param_9 = 3;

%Define bounds on parameters in the parameter estimation.
ub_powerlaw = 100;
ub_gomp = 10;
lb_growthrate = 0;
lb_decayrate = 10^(-4);
ub_cc = 10^20;
ub_tau = 7; 
lb_tau = 0;

%Define the smallest possible error.
tol = 10^(-3);

%% Definition of matrices storing results

no_organoids = zeros(1,1);

params = zeros(3,2,1);
params_2 = zeros(3,2,1);
params_3 = zeros(3,2,1);
params_4 = zeros(3,2,1);
params_5 = zeros(3,3,1);
params_6 = zeros(3,3,1);
params_7 = zeros(3,3,1);
params_8 = zeros(3,3,1);
params_9 = zeros(3,3,1);

error_vec_bic = zeros(2,6);
error_vec_2_bic = zeros(2,6);
error_vec_3_bic = zeros(2,6);
error_vec_4_bic = zeros(2,6);
error_vec_5_bic = zeros(2,6);
error_vec_6_bic = zeros(2,6);
error_vec_7_bic = zeros(2,6);
error_vec_8_bic = zeros(2,6);
error_vec_9_bic = zeros(2,6);

error_vec = zeros(2,5);
error_vec_2 = zeros(2,5);
error_vec_3 = zeros(2,5);
error_vec_4 = zeros(2,5);
error_vec_5 = zeros(2,5);
error_vec_6 = zeros(2,5);
error_vec_7 = zeros(2,5);
error_vec_8 = zeros(2,5);
error_vec_9 = zeros(2,5);

error_vec_norm = zeros(2,5);
error_vec_2_norm = zeros(2,5);
error_vec_3_norm = zeros(2,5);
error_vec_4_norm = zeros(2,5);
error_vec_5_norm = zeros(2,5);
error_vec_6_norm = zeros(2,5);
error_vec_7_norm = zeros(2,5);
error_vec_8_norm = zeros(2,5);
error_vec_9_norm = zeros(2,5);

r = zeros(2,1);
z = zeros(2,9,1);
w = zeros(2,9,1);

x_0_all = zeros(1,n_timepoints,1);
areas_all = zeros(1,3,1);
id_s_all = zeros(1,3,1);
wells_all = zeros(1,1);

max_category = 9;
categories = zeros(n_exp,max_category);
means = zeros(1,9);
medians = zeros(1,9);
means_norm = zeros(1,9);
medians_norm = zeros(1,9);
means_bic = zeros(1,9);
medians_bic = zeros(1,9);

%% Create datasets for analysis

%The tracking data is stored in a cell array called data_tracking.
%data_tracking{l} stores the data from well number l.
for count_dataset = 1:n_exp
numbers = [];
    id_s = [];
    wells = [];
    for l = cell2mat(indices(count_dataset,:))
        number = zeros(1,3);
        count = 0;
        for m = min(data_tracking{l}.TrackID):max(data_tracking{l}.TrackID)
            count = count+1;
            for ell = 0:3
                if ~isempty(find(data_tracking{l}.TrackID == m & data_tracking{l}.X1_CellTime == ell))
                    a = find(data_tracking{l}.TrackID == m & data_tracking{l}.X1_CellTime == ell);
                    number(count,ell+1) = data_tracking{l}.Y1_CellNumberOfVesiclesVesicleType_VesiclesGFP(a(1));
                end
            end
        end 
        a = find(number(:,1) == 0 & number(:,2) == 0 & number(:,3) == 0 & number(:,4) == 0);
        number(a,:) = [];
        numbers = [numbers;number];
    end

    %Only consider organoids that are growing over time.
    a = zeros(size(numbers,1),1);
    for k=1:size(numbers,1)
        if numbers(k,2) < numbers(k,1) || numbers(k,3) < numbers(k,2) || numbers(k,4) < numbers(k,3)
            a(k) = 1;
        end
    end
    numbers(find(a==1),:) = [];

    %Start each organoid with a single cell.
    x_0 = [ones(size(numbers,1),1),numbers];

    %Record the number of organoids in each dataset.
    no_organoids(count_dataset) = size(x_0,1);
    
    %Conduct growth modeling for each individual organoid.
    %A parallel for loop is used to speed up computation.
    %Taking advantage of the speed-up requires the MATLAB Parallel
    %Computing Toolbox.
    parfor k=1:no_organoids(count_dataset)
        %Store organoid cell numbers across all datasets.
        x_0_all(k,:,count_dataset) = x_0(k,:);
        %time_series contains the data for a single organoid at a time.
        time_series = x_0(k,:);

        %Suppress the displayed output of fmincon.
        options = optimoptions(@fmincon,'Display','off');

        %Exponential model fitting
        error = Inf;
            for ell=1:iter
                lb = [lb_tau,lb_growthrate];
                ub = [ub_tau,ub_gomp];
                initial_guess = (ub-lb).*rand(1,2);
                [x_temp, error_temp] = fmincon(@(x) loss_function_exponential(x_0(k,1),x(1),x(2),time_series),initial_guess,[],[],[],[],lb,ub,[],options);
                if error_temp <= error
                    x = x_temp;
                    error = error_temp;
                end
            end  
        params(k,:,count_dataset) = x;
        %Define the smallest possible error.
        if error < tol^2
            error = tol^2;
        end
        error_vec(k,count_dataset) = error;
        %Compute the normalized fitting error defined in Section 2.8.
        error_vec_norm(k,count_dataset) = sqrt(error/mean(x_0(k,:))^2);

        %Power law model fitting (gamma = 1/2)       
        error_2 = Inf;
            for ell=1:iter
                lb = [lb_tau,lb_growthrate];
                ub = [ub_tau,ub_powerlaw];
                initial_guess = (ub-lb).*rand(1,2);
                initial_guess(2) = 10^(log10(ub_powerlaw)*rand);
                [x_temp, error_temp] = fmincon(@(x) loss_function_powerlaw(x_0(k,1),x(1),x(2),1/2,time_series),initial_guess,[],[],[],[],lb,ub,[],options);
                if error_temp <= error_2
                    x = x_temp;
                    error_2 = error_temp;
                end
            end  
        params_2(k,:,count_dataset) = x;
        if error_2 < tol^2
            error_2 = tol^2;
        end
        error_vec_2(k,count_dataset) = error_2;
        error_vec_2_norm(k,count_dataset) = sqrt(error_2/mean(x_0(k,:))^2);

        %Power law model fitting (gamma = 2/3)
        error_3 = Inf;
            for ell=1:iter
                lb = [lb_tau,lb_growthrate];
                ub = [ub_tau,ub_powerlaw];
                initial_guess = (ub-lb).*rand(1,2);
                initial_guess(2) = 10^(log10(ub_powerlaw)*rand);
                [x_temp, error_temp] = fmincon(@(x) loss_function_powerlaw(x_0(k,1),x(1),x(2),2/3,time_series),initial_guess,[],[],[],[],lb,ub,[],options);
                if error_temp <= error_3
                    x = x_temp;
                    error_3 = error_temp;
                end
            end  
        params_3(k,:,count_dataset) = x;
        if error_3 < tol^2
            error_3 = tol^2;
        end
        error_vec_3(k,count_dataset) = error_3;
        error_vec_3_norm(k,count_dataset) = sqrt(error_3/mean(x_0(k,:))^2);

        %Power law model fitting (gamma = 3/4)
        error_4 = Inf;
            for ell=1:iter
                lb = [lb_tau,lb_growthrate];
                ub = [ub_tau,ub_powerlaw];
                initial_guess = (ub-lb).*rand(1,2);
                initial_guess(2) = 10^(log10(ub_powerlaw)*rand);
                [x_temp, error_temp] = fmincon(@(x) loss_function_powerlaw(x_0(k,1),x(1),x(2),3/4,time_series),initial_guess,[],[],[],[],lb,ub,[],options);
                if error_temp <= error_4
                    x = x_temp;
                    error_4 = error_temp;
                end
            end  
        params_4(k,:,count_dataset) = x;
        if error_4 < tol^2
            error_4 = tol^2;
        end
        error_vec_4(k,count_dataset) = error_4;
        error_vec_4_norm(k,count_dataset) = sqrt(error_4/mean(x_0(k,:))^2);

        options = optimoptions(@fmincon,'Display','off');
        
        %Gompertz model fitting     
        error_5 = Inf;
        for ell=1:iter
            lb = [lb_tau,lb_growthrate,lb_decayrate];
            ub = [ub_tau,ub_gomp,ub_gomp];
            initial_guess = (ub-lb).*rand(1,3);
            [x_temp, error_temp] = fmincon(@(x) loss_function_Gompertz(x_0(k,1),x(1),x(2),x(3),time_series),initial_guess,[],[],[],[],lb,ub,[],options);
            if error_temp <= error_5
                x = x_temp;
                error_5 = error_temp;
            end
        end
        params_5(k,:,count_dataset) = x;
        if error_5 < tol^2
            error_5 = tol^2;
        end
        error_vec_5(k,count_dataset) = error_5;
        error_vec_5_norm(k,count_dataset) = sqrt(error_5/mean(x_0(k,:))^2);

        options = optimoptions(@fmincon,'Display','off');

        %Logistic model fitting
        error_6 = Inf;
        for ell=1:iter
            lb = [lb_tau,lb_growthrate,0];
            ub = [ub_tau,ub_gomp,ub_cc];
            initial_guess = (ub-lb).*rand(1,3);
            initial_guess(3) = 10^(log10(ub_cc)*rand);
            [x_temp, error_temp] = fmincon(@(x) loss_function_logistic(x_0(k,1),x(1),x(2),x(3),time_series),initial_guess,[],[],[],[],lb,ub,[],options);
            if error_temp <= error_6
                x = x_temp;
                error_6 = error_temp;
            end
        end
        params_6(k,:,count_dataset) = x;
        if error_6 < tol^2
            error_6 = tol^2;
        end
        error_vec_6(k,count_dataset) = error_6;
        error_vec_6_norm(k,count_dataset) = sqrt(error_6/mean(x_0(k,:))^2);
    
        %von Bertalanffy model fitting (gamma = 1/2)
        error_7 = Inf;
        for ell=1:iter
            lb = [lb_tau,lb_growthrate,lb_decayrate];
            ub = [ub_tau,ub_powerlaw,ub_powerlaw];
            initial_guess = (ub-lb).*rand(1,3);
            initial_guess(2) = 10^(log10(ub_powerlaw)*rand);
            A = [];
            b = [];
            [x_temp, error_temp] = fmincon(@(x) loss_function_vonBertalanffy(x_0(k,1),x(1),x(2),x(3),1/2,time_series),initial_guess,[],[],[],[],lb,ub,[],options);
            if error_temp <= error_7
                x = x_temp;
                error_7 = error_temp;
            end
        end
        params_7(k,:,count_dataset) = x;
        if error_7 < tol^2
            error_7 = tol^2;
        end
        error_vec_7(k,count_dataset) = error_7; 
        error_vec_7_norm(k,count_dataset) = sqrt(error_7/mean(x_0(k,:))^2); 

        %von Bertalanffy model fitting (gamma = 2/3)
        error_8 = Inf;
        for ell=1:iter
            lb = [lb_tau,lb_growthrate,lb_decayrate];
            ub = [ub_tau,ub_powerlaw,ub_powerlaw];
            initial_guess = (ub-lb).*rand(1,3);
            initial_guess(2) = 10^(log10(ub_powerlaw)*rand);
            A = [];
            b = [];
            [x_temp, error_temp] = fmincon(@(x) loss_function_vonBertalanffy(x_0(k,1),x(1),x(2),x(3),2/3,time_series),initial_guess,[],[],[],[],lb,ub,[],options);
            if error_temp <= error_8
                x = x_temp;
                error_8 = error_temp;
            end
        end
        params_8(k,:,count_dataset) = x;
        if error_8 < tol^2
            error_8 = tol^2;
        end
        error_vec_8(k,count_dataset) = error_8;  
        error_vec_8_norm(k,count_dataset) = sqrt(error_8/mean(x_0(k,:))^2);

        %von Bertalanffy model fitting (gamma = 3/4)
        error_9 = Inf;
        for ell=1:iter
            lb = [lb_tau,lb_growthrate,lb_decayrate];
            ub = [ub_tau,ub_powerlaw,ub_powerlaw];
            initial_guess = (ub-lb).*rand(1,3);
            initial_guess(2) = 10^(log10(ub_powerlaw)*rand);
            A = [];
            b = [];
            [x_temp, error_temp] = fmincon(@(x) loss_function_vonBertalanffy(x_0(k,1),x(1),x(2),x(3),3/4,time_series),initial_guess,A,b,[],[],lb,ub,[],options);
            if error_temp <= error_9
                x = x_temp;
                error_9 = error_temp;
            end
        end
        params_9(k,:,count_dataset) = x;
        if error_9 < tol^2
            error_9 = tol^2;
        end
        error_vec_9(k,count_dataset) = error_9;  
        error_vec_9_norm(k,count_dataset) = sqrt(error_9/mean(x_0(k,:))^2);

        %Compute the BIC for each individual organoid and each model.
        error_vec_bic(k,count_dataset) = n_datapoints*(1+log(2*pi))+n_datapoints*log(error_vec(k,count_dataset)/n_datapoints)+(n_param+1)*log(n_datapoints);
        error_vec_2_bic(k,count_dataset) = n_datapoints*(1+log(2*pi))+n_datapoints*log(error_vec_2(k,count_dataset)/n_datapoints)+(n_param_2+1)*log(n_datapoints);
        error_vec_3_bic(k,count_dataset) = n_datapoints*(1+log(2*pi))+n_datapoints*log(error_vec_3(k,count_dataset)/n_datapoints)+(n_param_3+1)*log(n_datapoints);
        error_vec_4_bic(k,count_dataset) = n_datapoints*(1+log(2*pi))+n_datapoints*log(error_vec_4(k,count_dataset)/n_datapoints)+(n_param_4+1)*log(n_datapoints);
        error_vec_5_bic(k,count_dataset) = n_datapoints*(1+log(2*pi))+n_datapoints*log(error_vec_5(k,count_dataset)/n_datapoints)+(n_param_5+1)*log(n_datapoints);
        error_vec_6_bic(k,count_dataset) = n_datapoints*(1+log(2*pi))+n_datapoints*log(error_vec_6(k,count_dataset)/n_datapoints)+(n_param_6+1)*log(n_datapoints);
        error_vec_7_bic(k,count_dataset) = n_datapoints*(1+log(2*pi))+n_datapoints*log(error_vec_7(k,count_dataset)/n_datapoints)+(n_param_7+1)*log(n_datapoints);
        error_vec_8_bic(k,count_dataset) = n_datapoints*(1+log(2*pi))+n_datapoints*log(error_vec_8(k,count_dataset)/n_datapoints)+(n_param_8+1)*log(n_datapoints);
        error_vec_9_bic(k,count_dataset) = n_datapoints*(1+log(2*pi))+n_datapoints*log(error_vec_9(k,count_dataset)/n_datapoints)+(n_param_9+1)*log(n_datapoints);

        %Gather modeling results for all growth models.
        z(k,:,count_dataset) = [error_vec_bic(k,count_dataset), error_vec_2_bic(k,count_dataset), error_vec_3_bic(k,count_dataset), error_vec_4_bic(k,count_dataset), error_vec_5_bic(k,count_dataset), error_vec_6_bic(k,count_dataset), error_vec_7_bic(k,count_dataset), error_vec_8_bic(k,count_dataset), error_vec_9_bic(k,count_dataset)];
        w(k,:,count_dataset) = [error_vec(k,count_dataset), error_vec_2(k,count_dataset), error_vec_3(k,count_dataset), error_vec_4(k,count_dataset), error_vec_5(k,count_dataset), error_vec_6(k,count_dataset), error_vec_7(k,count_dataset), error_vec_8(k,count_dataset), error_vec_9(k,count_dataset)];
        w_norm(k,:,count_dataset) = [error_vec_norm(k,count_dataset), error_vec_2_norm(k,count_dataset), error_vec_3_norm(k,count_dataset), error_vec_4_norm(k,count_dataset), error_vec_5_norm(k,count_dataset), error_vec_6_norm(k,count_dataset), error_vec_7_norm(k,count_dataset), error_vec_8_norm(k,count_dataset), error_vec_9_norm(k,count_dataset)];
    end
end

%Identify the organoids referred to as "exponential" in the paper.
for count_dataset = 1:n_exp
    for k=1:no_organoids(count_dataset)
        r(k,count_dataset) = abs(params_5(k,3,count_dataset)-lb_decayrate)<lb_decayrate*10^(-2);
    end
end

%Collect results on model fit quality.
for k=1:count_dataset
    index_set = 1:no_organoids(k);

    for ell = 1:max_category
        categories(k,ell) = size(find(z(1:no_organoids(k),ell,k) == min(z(1:no_organoids(k),1:max_category,k)')'),1);
    end

    means(k,:) = mean(w(index_set,:,k));
    means_norm(k,:) = mean(w_norm(index_set,:,k));
    means_bic(k,:) = mean(z(index_set,:,k));
    medians(k,:) = median(w(index_set,:,k));
    medians_norm(k,:) = median(w_norm(index_set,:,k));
    medians_bic(k,:) = median(z(index_set,:,k));
end

toc

%Define the sum of least squares to be minimized for each growth model.
function sum = loss_function_exponential(v_0,t_0,b,a)
    sum = 0;
    times = [-t_0,0,1,3,6];
    for u=1:5
        sum = sum + (a(u)-v_0*exp(b*(t_0+times(u))))^2;
    end
end

function sum = loss_function_powerlaw(v_0,t_0,b,gamma,a)
    sum = 0;
    times = [-t_0,0,1,3,6];
    for u=1:5
        sum = sum + (a(u)-(v_0^(1-gamma)+(1-gamma)*b*(t_0+times(u)))^(1/(1-gamma)))^2;
    end
end

function sum = loss_function_Gompertz(v_0,t_0,b,c,a)
    sum = 0;
    times = [-t_0,0,1,3,6];
    for u=1:5
        sum = sum + (a(u)-v_0*exp(b/c*(1-exp(-c*(t_0+times(u))))))^2;
    end
end

function sum = loss_function_logistic(v_0,t_0,b,k,a)
    sum = 0;
    times = [-t_0,0,1,3,6];
    for u=1:5
        sum = sum + (a(u)-v_0*k/(v_0+(k-v_0)*exp(-b*(t_0+times(u)))))^2;
    end
end

function sum = loss_function_vonBertalanffy(v_0,t_0,b,c,gamma,a)
    sum = 0;
    times = [-t_0,0,1,3,6];
    for u=1:5
        sum = sum + (a(u)-(b/c+(v_0^(1-gamma)-b/c)*exp(-c*(1-gamma)*(t_0+times(u))))^(1/(1-gamma)))^2;
    end
end