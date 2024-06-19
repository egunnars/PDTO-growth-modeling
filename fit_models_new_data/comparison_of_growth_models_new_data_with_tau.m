%This code fits the 9 growth models considered in the paper on a new
%user-provided dataset. Each model includes the delay (starting time)
%parameter tau as explained in the accompanying paper.

%It is assumed that the data is in a matrix called "data", where the number
%of rows is the number of organoids in the dataset and the number of
%columns is the number of time points data was collected at.

%% Definition of main parameters
%The user needs to redefine the variables in this section according to their
%own dataset before running the code.

%Artifical data is defined below to show the structure of the data matrix.
%The initial organoid size, which is assumed to be identical for all
%organoids, is not included in the data matrix.
%If the user provides their own dataset, the following three lines should
%be commented out before running the code.
data = [150,200,250;
    170,240,310;
    120,140,160];

%Time points at which measurements were collected (row vector). The
%number of time points should correspond to the number of columns in the
%"data" matrix provided by the user.
%In the accompanying paper, the data was collected on Days 0, 3 and 5.
timepoints = [0,3,5];

%Define the initial size of each organoid. The initial size will be added
%to the time series data for each organoid later in the code.
initial_size = 100;

%The user needs to define an upper bound on the delay parameter tau.
%The upper bound should be the number of time units from the start of the
%experiment until the first measurement is taken. In the accompanying
%paper, the organoids were seeded on Day -4 and data was first collected on
%Day 0, which means there were were 4 days from the start of the experiment
%until the first measurement was taken.
ub_tau = 4; 
lb_tau = 0;

%Here upper and lower bounds are defined on other parameters involved in
%the growth modeling.
ub_powerlaw = 100;
ub_gomp = 10;
lb_growthrate = 0;
lb_decayrate = 10^(-4);
ub_cc = 10^20;

%Number of optimization runs in the least squares estimation for each
%growth model.
iter = 1000;

%Smallest possible error. See Section 2.8 of the paper.
tol = 10^(-3);

%Number of fitted datapoints for calculation of BIC.
n_datapoints = size(timepoints,2);

%% Definition of other parameters and matrices

%Number of parameters in the 9 growth models we consider,
%including the tau (starting time) parameter.
%The models considered are shown in Table 2 of the paper.
n_param = 2;
n_param_2 = 2;
n_param_3 = 2;
n_param_4 = 2;
n_param_5 = 3;
n_param_6 = 3;
n_param_7 = 3;
n_param_8 = 3;
n_param_9 = 3;

%Definition of matrices storing results.
params = zeros(3,2);
params_2 = zeros(3,2);
params_3 = zeros(3,2);
params_4 = zeros(3,2);
params_5 = zeros(3,3);
params_6 = zeros(3,3);
params_7 = zeros(3,3);
params_8 = zeros(3,3);
params_9 = zeros(3,3);

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
z = zeros(2,9);
w = zeros(2,9);
w_norm = zeros(2,9);

max_category = 9;
categories = zeros(1,max_category);
means = zeros(1,9);
medians = zeros(1,9);
means_norm = zeros(1,9);
medians_norm = zeros(1,9);
means_bic = zeros(1,9);
medians_bic = zeros(1,9);

%% Growth modeling

%Here, the initial organoid size is added to the time series for each
%organoid.
data = [initial_size*ones(size(data,1),1),data];
x_0 = data;
no_organoids = size(x_0,1);
    
%Conduct growth modeling for each individual organoid.
%A parallel for loop is used to speed up computation.
%Taking advantage of the speed-up requires the MATLAB Parallel
%Computing Toolbox.
parfor k=1:no_organoids
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
            [x_temp, error_temp] = fmincon(@(x) loss_function_exponential(x_0(k,1),x(1),x(2),time_series,timepoints),initial_guess,[],[],[],[],lb,ub,[],options);
            if error_temp <= error
                x = x_temp;
                error = error_temp;
            end
        end  
    params(k,:) = x;
    %Define the smallest possible error.
    if error < tol^2
        error = tol^2;
    end
    error_vec(k) = error;
    %Compute the normalized fitting error defined in Section 2.8.
    error_vec_norm(k) = sqrt(error/mean(x_0(k,:))^2);

    %Power law model fitting (gamma = 1/2)       
    error_2 = Inf;
        for ell=1:iter
            lb = [lb_tau,lb_growthrate];
            ub = [ub_tau,ub_powerlaw];
            initial_guess = (ub-lb).*rand(1,2);
            initial_guess(2) = 10^(log10(ub_powerlaw)*rand);
            [x_temp, error_temp] = fmincon(@(x) loss_function_powerlaw(x_0(k,1),x(1),x(2),1/2,time_series,timepoints),initial_guess,[],[],[],[],lb,ub,[],options);
            if error_temp <= error_2
                x = x_temp;
                error_2 = error_temp;
            end
        end  
    params_2(k,:) = x;
    if error_2 < tol^2
        error_2 = tol^2;
    end
    error_vec_2(k) = error_2;
    error_vec_2_norm(k) = sqrt(error_2/mean(x_0(k,:))^2);

    %Power law model fitting (gamma = 2/3)
    error_3 = Inf;
        for ell=1:iter
            lb = [lb_tau,lb_growthrate];
            ub = [ub_tau,ub_powerlaw];
            initial_guess = (ub-lb).*rand(1,2);
            initial_guess(2) = 10^(log10(ub_powerlaw)*rand);
            [x_temp, error_temp] = fmincon(@(x) loss_function_powerlaw(x_0(k,1),x(1),x(2),2/3,time_series,timepoints),initial_guess,[],[],[],[],lb,ub,[],options);
            if error_temp <= error_3
                x = x_temp;
                error_3 = error_temp;
            end
        end  
    params_3(k,:) = x;
    if error_3 < tol^2
        error_3 = tol^2;
    end
    error_vec_3(k) = error_3;
    error_vec_3_norm(k) = sqrt(error_3/mean(x_0(k,:))^2);

    %Power law model fitting (gamma = 3/4)
    error_4 = Inf;
        for ell=1:iter
            lb = [lb_tau,lb_growthrate];
            ub = [ub_tau,ub_powerlaw];
            initial_guess = (ub-lb).*rand(1,2);
            initial_guess(2) = 10^(log10(ub_powerlaw)*rand);
            [x_temp, error_temp] = fmincon(@(x) loss_function_powerlaw(x_0(k,1),x(1),x(2),3/4,time_series,timepoints),initial_guess,[],[],[],[],lb,ub,[],options);
            if error_temp <= error_4
                x = x_temp;
                error_4 = error_temp;
            end
        end  
    params_4(k,:) = x;
    if error_4 < tol^2
        error_4 = tol^2;
    end
    error_vec_4(k) = error_4;
    error_vec_4_norm(k) = sqrt(error_4/mean(x_0(k,:))^2);

    options = optimoptions(@fmincon,'Display','off');
    
    %Gompertz model fitting     
    error_5 = Inf;
    for ell=1:iter
        lb = [lb_tau,lb_growthrate,lb_decayrate];
        ub = [ub_tau,ub_gomp,ub_gomp];
        initial_guess = (ub-lb).*rand(1,3);
        [x_temp, error_temp] = fmincon(@(x) loss_function_Gompertz(x_0(k,1),x(1),x(2),x(3),time_series,timepoints),initial_guess,[],[],[],[],lb,ub,[],options);
        if error_temp <= error_5
            x = x_temp;
            error_5 = error_temp;
        end
    end
    params_5(k,:) = x;
    if error_5 < tol^2
        error_5 = tol^2;
    end
    error_vec_5(k) = error_5;
    error_vec_5_norm(k) = sqrt(error_5/mean(x_0(k,:))^2);

    options = optimoptions(@fmincon,'Display','off');

    %Logistic model fitting
    error_6 = Inf;
    for ell=1:iter
        lb = [lb_tau,lb_growthrate,0];
        ub = [ub_tau,ub_gomp,ub_cc];
        initial_guess = (ub-lb).*rand(1,3);
        initial_guess(3) = 10^(log10(ub_cc)*rand);
        [x_temp, error_temp] = fmincon(@(x) loss_function_logistic(x_0(k,1),x(1),x(2),x(3),time_series,timepoints),initial_guess,[],[],[],[],lb,ub,[],options);
        if error_temp <= error_6
            x = x_temp;
            error_6 = error_temp;
        end
    end
    params_6(k,:) = x;
    if error_6 < tol^2
        error_6 = tol^2;
    end
    error_vec_6(k) = error_6;
    error_vec_6_norm(k) = sqrt(error_6/mean(x_0(k,:))^2);

    %von Bertalanffy model fitting (gamma = 1/2)
    error_7 = Inf;
    for ell=1:iter
        lb = [lb_tau,lb_growthrate,lb_decayrate];
        ub = [ub_tau,ub_powerlaw,ub_powerlaw];
        initial_guess = (ub-lb).*rand(1,3);
        initial_guess(2) = 10^(log10(ub_powerlaw)*rand);
        A = [];
        b = [];
        [x_temp, error_temp] = fmincon(@(x) loss_function_vonBertalanffy(x_0(k,1),x(1),x(2),x(3),1/2,time_series,timepoints),initial_guess,[],[],[],[],lb,ub,[],options);
        if error_temp <= error_7
            x = x_temp;
            error_7 = error_temp;
        end
    end
    params_7(k,:) = x;
    if error_7 < tol^2
        error_7 = tol^2;
    end
    error_vec_7(k) = error_7; 
    error_vec_7_norm(k) = sqrt(error_7/mean(x_0(k,:))^2); 

    %von Bertalanffy model fitting (gamma = 2/3)
    error_8 = Inf;
    for ell=1:iter
        lb = [lb_tau,lb_growthrate,lb_decayrate];
        ub = [ub_tau,ub_powerlaw,ub_powerlaw];
        initial_guess = (ub-lb).*rand(1,3);
        initial_guess(2) = 10^(log10(ub_powerlaw)*rand);
        A = [];
        b = [];
        [x_temp, error_temp] = fmincon(@(x) loss_function_vonBertalanffy(x_0(k,1),x(1),x(2),x(3),2/3,time_series,timepoints),initial_guess,[],[],[],[],lb,ub,[],options);
        if error_temp <= error_8
            x = x_temp;
            error_8 = error_temp;
        end
    end
    params_8(k,:) = x;
    if error_8 < tol^2
        error_8 = tol^2;
    end
    error_vec_8(k) = error_8;  
    error_vec_8_norm(k) = sqrt(error_8/mean(x_0(k,:))^2);

    %von Bertalanffy model fitting (gamma = 3/4)
    error_9 = Inf;
    for ell=1:iter
        lb = [lb_tau,lb_growthrate,lb_decayrate];
        ub = [ub_tau,ub_powerlaw,ub_powerlaw];
        initial_guess = (ub-lb).*rand(1,3);
        initial_guess(2) = 10^(log10(ub_powerlaw)*rand);
        A = [];
        b = [];
        [x_temp, error_temp] = fmincon(@(x) loss_function_vonBertalanffy(x_0(k,1),x(1),x(2),x(3),3/4,time_series,timepoints),initial_guess,A,b,[],[],lb,ub,[],options);
        if error_temp <= error_9
            x = x_temp;
            error_9 = error_temp;
        end
    end
    params_9(k,:) = x;
    if error_9 < tol^2
        error_9 = tol^2;
    end
    error_vec_9(k) = error_9;  
    error_vec_9_norm(k) = sqrt(error_9/mean(x_0(k,:))^2);

    %Compute the BIC for each individual organoid and each model.
    error_vec_bic(k) = n_datapoints*(1+log(2*pi))+n_datapoints*log(error_vec(k)/n_datapoints)+(n_param+1)*log(n_datapoints);
    error_vec_2_bic(k) = n_datapoints*(1+log(2*pi))+n_datapoints*log(error_vec_2(k)/n_datapoints)+(n_param_2+1)*log(n_datapoints);
    error_vec_3_bic(k) = n_datapoints*(1+log(2*pi))+n_datapoints*log(error_vec_3(k)/n_datapoints)+(n_param_3+1)*log(n_datapoints);
    error_vec_4_bic(k) = n_datapoints*(1+log(2*pi))+n_datapoints*log(error_vec_4(k)/n_datapoints)+(n_param_4+1)*log(n_datapoints);
    error_vec_5_bic(k) = n_datapoints*(1+log(2*pi))+n_datapoints*log(error_vec_5(k)/n_datapoints)+(n_param_5+1)*log(n_datapoints);
    error_vec_6_bic(k) = n_datapoints*(1+log(2*pi))+n_datapoints*log(error_vec_6(k)/n_datapoints)+(n_param_6+1)*log(n_datapoints);
    error_vec_7_bic(k) = n_datapoints*(1+log(2*pi))+n_datapoints*log(error_vec_7(k)/n_datapoints)+(n_param_7+1)*log(n_datapoints);
    error_vec_8_bic(k) = n_datapoints*(1+log(2*pi))+n_datapoints*log(error_vec_8(k)/n_datapoints)+(n_param_8+1)*log(n_datapoints);
    error_vec_9_bic(k) = n_datapoints*(1+log(2*pi))+n_datapoints*log(error_vec_9(k)/n_datapoints)+(n_param_9+1)*log(n_datapoints);

    %Gather modeling results for all growth models.
    z(k,:) = [error_vec_bic(k), error_vec_2_bic(k), error_vec_3_bic(k), error_vec_4_bic(k), error_vec_5_bic(k), error_vec_6_bic(k), error_vec_7_bic(k), error_vec_8_bic(k), error_vec_9_bic(k)];
    w(k,:) = [error_vec(k), error_vec_2(k), error_vec_3(k), error_vec_4(k), error_vec_5(k), error_vec_6(k), error_vec_7(k), error_vec_8(k), error_vec_9(k)];
    w_norm(k,:) = [error_vec_norm(k), error_vec_2_norm(k), error_vec_3_norm(k), error_vec_4_norm(k), error_vec_5_norm(k), error_vec_6_norm(k), error_vec_7_norm(k), error_vec_8_norm(k), error_vec_9_norm(k)];
end

%Identify the organoids referred to as "exponential" in the paper.
for k=1:no_organoids
    r(k) = abs(params_5(k,3)-lb_decayrate)<lb_decayrate*10^(-2);
end

%Collect results on model fit quality.
index_set = 1:no_organoids;

for ell = 1:max_category
    categories(k,ell) = size(find(z(1:no_organoids,ell) == min(z(1:no_organoids,1:max_category)')'),1);
end

%means_bic gives the matrix shown in Table 2 of the paper.
%means_norm gives the matrix shown in Table 3.
means(1,:) = mean(w(index_set,:));
means_norm(1,:) = mean(w_norm(index_set,:));
means_bic(1,:) = mean(z(index_set,:));
medians(1,:) = median(w(index_set,:));
medians_norm(1,:) = median(w_norm(index_set,:));
medians_bic(1,:) = median(z(index_set,:));

%Define the sum of least squares to be minimized for each growth model.
function sum = loss_function_exponential(v_0,t_0,b,a,timepoints)
    sum = 0;
    times = [timepoints(1)-t_0,timepoints];
    for u=1:size(times,2)
        sum = sum + (a(u)-v_0*exp(b*(t_0-timepoints(1)+times(u))))^2;
    end
end

function sum = loss_function_powerlaw(v_0,t_0,b,gamma,a,timepoints)
    sum = 0;
    times = [timepoints(1)-t_0,timepoints];
    for u=1:size(times,2)
        sum = sum + (a(u)-(v_0^(1-gamma)+(1-gamma)*b*(t_0-timepoints(1)+times(u)))^(1/(1-gamma)))^2;
    end
end

function sum = loss_function_Gompertz(v_0,t_0,b,c,a,timepoints)
    sum = 0;
    times = [timepoints(1)-t_0,timepoints];
    for u=1:size(times,2)
        sum = sum + (a(u)-v_0*exp(b/c*(1-exp(-c*(t_0-timepoints(1)+times(u))))))^2;
    end
end

function sum = loss_function_logistic(v_0,t_0,b,k,a,timepoints)
    sum = 0;
    times = [timepoints(1)-t_0,timepoints];
    for u=1:size(times,2)
        sum = sum + (a(u)-v_0*k/(v_0+(k-v_0)*exp(-b*(t_0-timepoints(1)+times(u)))))^2;
    end
end

function sum = loss_function_vonBertalanffy(v_0,t_0,b,c,gamma,a,timepoints)
    sum = 0;
    times = [timepoints(1)-t_0,timepoints];
    for u=1:size(times,2)
        sum = sum + (a(u)-(b/c+(v_0^(1-gamma)-b/c)*exp(-c*(1-gamma)*(t_0-timepoints(1)+times(u))))^(1/(1-gamma)))^2;
    end
end