%%%Thrust Magnetic Bearing System Identification
%%%RLS with Forgetting Factor Lambda = [0.1  0.975  0.95]

%Pertama-tama set forgetting factor 
%Forgetting Factor (Lambda) akan terbagi menjadi 3 yaitu : 
list_lambda=[1 0.975 0.95]; 
N=length(thrust_in);
Qt=thrust_out;  %Input TMB
Qs=thrust_in;   %Output TMB
%Matrix konvolusi dengan 16-tap filter
%Asumsi data sebelum adalah rata-rata masukan TMB 
Regressor=...
    [Qs';...
    [0 (Qs(1:N-1))'];...
    [0 0 (Qs(1:N-2))'];...
    [0 0 0 (Qs(1:N-3))'];...
    [0 0 0 0 (Qs(1:N-4))'];...
    [0 0 0 0 0 (Qs(1:N-5))'];...
    [0 0 0 0 0 0 (Qs(1:N-6))'];...
    [0 0 0 0 0 0 0 (Qs(1:N-7))'];...
    [0 0 0 0 0 0 0 0 (Qs(1:N-8))'];...
    [0 0 0 0 0 0 0 0 0 (Qs(1:N-9))'];...
    [0 0 0 0 0 0 0 0 0 0 (Qs(1:N-10))'];...
    [0 0 0 0 0 0 0 0 0 0 0 (Qs(1:N-11))'];...
    [0 0 0 0 0 0 0 0 0 0 0 0 (Qs(1:N-12))'];...
    [0 0 0 0 0 0 0 0 0 0 0 0 0 (Qs(1:N-13))'];...
    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 (Qs(1:N-14))'];...
    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 (Qs(1:N-15))'];]';
%Solusi Initial dengan menggunakan metode Least-Squares 
%Solusi Initial menggunakan 4 data pertama 
Init_Regressor=Regressor(1:16,:); 
Init_y=Qs(1:16); 
all_theta_est=[];

%Recursive Least Squares tiap Lambda 
for count_lambda = 1:length(list_lambda) 
    lambda=list_lambda(count_lambda); 
    power_lambda=(1/lambda).^[1:16]; 
    Lambda = diag(power_lambda); 
    back=(Init_Regressor.')*Lambda*Init_y; 
    front=(Init_Regressor.')*Lambda*Init_Regressor; 
    init_theta_est=inv(front)*back 
    Pk = inv((Init_Regressor.')*Init_Regressor); 
    theta_k_est = init_theta_est;   
                        
    num_curr_obs = length(Init_y); 
    collect_theta_est_lambda=init_theta_est; 

    for count=1:(length(Qs)-length(Init_y)) 
        
        output_k=Qs(num_curr_obs+count); 
        Regrow_k=Regressor(num_curr_obs+count,:); 
        Error_k=output_k-Regrow_k*theta_k_est; 

         
        num = Pk*((Regrow_k.')*Regrow_k)*Pk; 
        denom = lambda+ Regrow_k*Pk*(Regrow_k.'); 
        oldPk=Pk; 
        Pk = (oldPk-(num/denom))/lambda; 

        old_theta=theta_k_est; 
        theta_k_est = old_theta +Pk*Regrow_k'* Error_k; 
    
        collect_theta_est_lambda=[collect_theta_est_lambda theta_k_est]; 
    end 
    final_RLS_theta_est=theta_k_est
    all_theta_est=[all_theta_est;collect_theta_est_lambda];
    end 