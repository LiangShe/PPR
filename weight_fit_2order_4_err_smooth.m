function [weight, CRF_par, r_pred] = weight_fit_2order_4_err_smooth (stim, resp, ind_train, ind_validation, weight_ini, learn_rate, n_iter, lamda)	
% polyfit only training data
% 
% the criteria for validation data during gradient descent is MSE
% 
% add Regularization term for the filters to get smooth one.
% lamda: weight of smooth term, tradeoff between fitting the responses accurately and 
% estimating smooth filters.
% 

%% initialize parameters
showfig=0;

learn_rate_new=learn_rate;

weight=weight_ini;

stim_train=stim(:,ind_train);

gLo=[0,0,1,0,0;    0,2,-8,2,0;    1,-8,20,-8,1;    0,2,-8,2,0;...
    0,0,1,0,0;]; % Gradient Operator of Laplacian Operator 
res2=length(weight);
res=sqrt(res2);
grad_smooth=zeros(res,res);

if showfig
figure(1);set(1,'position',[71 115 204 414]); colormap(gray)
h11=subplot(2,1,1,'Parent',1);
h12=subplot(2,1,2,'Parent',1);

figure(2);set(2,'position',[300 115 204 414]);
h21=subplot(2,1,1,'Parent',2);
h22=subplot(2,1,2,'Parent',2);
end


mse_validation_min = inf;
mse_validation=zeros(n_iter,1);

if showfig
    CC_validation_max = -1;
    CC_validation=zeros(n_iter,1);
    CC_train=zeros(n_iter,1);
    mse_train=zeros(n_iter,1);
end
%%
for ii=1:n_iter
    
    r_lin=stim'*weight;
    CRF_par=polyfit(r_lin(ind_train), resp(ind_train), 2); 
	r_pred = polyval(CRF_par,r_lin); 
    
	err=resp-r_pred;
    err_train=err(ind_train);
    
    if showfig
        err_validation=err(ind_validation);
        
        mse_train(ii)=mean(err_train.^2);
        mse_validation(ii)=mean(err_validation.^2);
        
        tmp=corrcoef(resp(ind_train),r_pred(ind_train));
        CC_train(ii)=tmp(1,2);
        
        tmp=corrcoef(resp(ind_validation),r_pred(ind_validation));
        CC_validation(ii)=tmp(1,2);
    end

    err_validation=err(ind_validation);
    mse_validation(ii)=mean(err_validation.^2);
  
	if mse_validation(ii)<=mse_validation_min
        
        mse_validation_min=mse_validation(ii);
        weight_old=weight;
		CRF_par_old=CRF_par;
		r_pred_old=r_pred;
		learn_rate_new=learn_rate_new*1.05;
        
        % update gradient
        ActFunc_deriv_train=2*CRF_par(1)*r_lin(ind_train)+CRF_par(2);
        grad=stim_train*(err_train.*ActFunc_deriv_train)/length(err_train);
        
        % grad_smooth_valid=conv2(reshape(weight,res,res),gLo,'valid');  % ignore edge
        % grad_smooth(3:end-2,3:end-2)= grad_smooth_valid * lamda;
        grad_smooth=conv2(reshape(weight,res,res),gLo,'same')* lamda;    % include edge
        
        weight=weight+learn_rate_new*(grad - grad_smooth(:));
        weight=weight/std(weight);
        
	else
		weight=weight_old;
		CRF_par=CRF_par_old;
		r_pred=r_pred_old;
		learn_rate_new=learn_rate_new*0.7;
        
        weight=weight+learn_rate_new*(grad - grad_smooth(:));
        weight=weight/std(weight);
        
	end

    if learn_rate_new<10^(-4)
        break
    end
    
if showfig

    imagesc(reshape(weight,res,res ),'Parent',h11)
    axis(h11,'off')
    axis(h11,'equal')
    plot(resp,'Parent',h12)
    hold(h12,'on')
    plot(r_pred,'r','Parent',h12)
    hold(h12,'off')

    plot(mse_train(1:ii),'Parent',h21)
    hold(h21,'on')
    plot(mse_validation(1:ii),'r','Parent',h21)
    hold(h21,'off')
    title(h21,'error')
    
    plot(CC_train(1:ii),'Parent',h22)
    hold(h22,'on')
    plot(CC_validation(1:ii),'r','Parent',h22)
    hold(h22,'off')
    title(h22,'CC')
    pause(0.1)
end

end

