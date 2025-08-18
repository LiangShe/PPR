function models = compute_validation_and_test_ccs( result_path, stim, resp, nseg, MaxSubunit )
% For PPR4 result
if nargin<4
    nseg = 20;
end
if nargin<5
    MaxSubunit=9;
end

count=0;
models.PPR4=[];

[y, x, z]=size(stim);
stim=reshape(stim,x*y,z);
stim=stim-mean(stim(:));
stim=stim/std(stim(:));

result_files = dir(fullfile (result_path, "PPR4_result*.mat"));

for i=1:length(result_files)
    result_file = fullfile(result_path, result_files(i).name);
    load(result_file,'cell_id','data_bw','lamda','Part','N');
    
    for k=1:MaxSubunit  % for each model
        nfilter=size(data_bw(k).filters,2);

        % prediction of validation data
        ind_validation=Part(N,:);
        zt = sum(ind_validation);
        r_pred=zeros(zt,1);
        for ifilter=1:nfilter
            r_lin = stim(:,ind_validation)'*data_bw(k).filters(:,ifilter);
            r_pred_i = polyval(data_bw(k).CRF_par(:,ifilter),r_lin);
            r_pred=r_pred+r_pred_i;
        end

        ccs=zeros(nseg,1);
        respt=resp(ind_validation);
        for j=1:nseg
            seg_ind= j:nseg:zt; % interleaved segmentation
            ccs(j)=corr(r_pred(seg_ind), respt(seg_ind));
        end

        % prediction of test data
        ind_test=Part(end,:);
        r_pred=zeros(sum(ind_test),1);
        for ifilter=1:nfilter
            r_lin = stim(:,ind_test)'*data_bw(k).filters(:,ifilter);
            r_pred_i = polyval(data_bw(k).CRF_par(:,ifilter),r_lin);
            r_pred=r_pred+r_pred_i;
        end
        test_cc = corr( r_pred, resp(ind_test) );
        
        count=count+1;
        models.PPR4(count).filters = data_bw(k).filters; 
        models.PPR4(count).n = nfilter; 
        models.PPR4(count).CRF_par = data_bw(k).CRF_par;
        models.PPR4(count).ccs = ccs;
        models.PPR4(count).lamda = lamda;
        models.PPR4(count).test_cc = test_cc;
        models.PPR4(count).Nth_part = N;
        models.PPR4(count).cell_id = cell_id;
                
    end
end
                
end