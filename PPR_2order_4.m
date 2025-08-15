function PPR_2order_4(stim, resp, lamda, cell_id, outputpath)
% PPR_2order_4(stim, resp, lamda, outputpath)
%

%%%%%%% PPR parameters %%%%%%%%%%%%%
n_filter_forward=20;
n_iter_forward=1000;
n_iter_backward=1000;
LR_f=0.01;
LR_b=0.01;
Npart=5; %  dividing data into N part for training validation and test data
polynomial_order='2nd order';

%
analysis_patch_proportion=[0 0 1 1]; %[top left height width]

Nrepeat=1;

% data partitioning
% seed = sum(100*clock);
seed = 21416; % use fixed seed for reproducibility

len=length(resp);
Part = DataPartitioning(len, Npart, seed);

[yr, xr, f]=size(stim);
analysis_patch=zeros(1,4);
analysis_patch(1)=round(analysis_patch_proportion(1)*(yr-1)+1);
analysis_patch(2)=round(analysis_patch_proportion(2)*(xr-1)+1);
analysis_patch(3)=round(analysis_patch_proportion(3)*(yr-1)+1);
analysis_patch(4)=round(analysis_patch_proportion(4)*(xr-1)+1);

if any(analysis_patch_proportion~=[0 0 1 1])
    stim=stim(analysis_patch(1):analysis_patch(1)+analysis_patch(3)-1,...
        analysis_patch(2):analysis_patch(2)+analysis_patch(4)-1,:);
end

[yr, xr, f]=size(stim);
stim=stim-mean(stim(:));
stim=stim/std(stim(:));


%%%%%%%%%%%%% PPR %%%%%%%%%%%%%%%%%


for i_repeat=1:Nrepeat
    fndate = datestr(now,'yyyy.mm.dd.HH.MM.SS');
    for N=1:Npart-1 % one part was reserved as test data
        fprintf('computing %d/%d repeats, %d/%d parts ...\n',i_repeat,Nrepeat,N,Npart-1); 
        stim=reshape(stim,yr*xr,f);

        ind_validation=Part(N,:);
        ind_train=false(1,f);
        for i=[1:N-1 N+1:Npart-1]
            ind_train=ind_train | Part(i,:);
        end
        nvalidation=sum(ind_validation);



        %% %%%%%%%%%%%%% forward step %%%%%%%%%%%%%%%%

        filters=zeros(yr*xr,n_filter_forward);
        CRF_par=zeros(3,n_filter_forward);   % 2nd order polynomial parameter
        amp=zeros(n_filter_forward,1);
        r_pred=zeros(f,n_filter_forward);

        % h=figure('Position',[528 111 557 419]);
        % colormap(gray)
        % suptitle(num2str(N))


        r_residue=resp;
        for ii=1:n_filter_forward

        	weight_ini=randn(yr*xr,1);
            weight_ini=weight_ini/std(weight_ini);
        	[filters(:,ii), CRF_par(:,ii), r_pred(:,ii)]  =  ...
                weight_fit_2order_4_err_smooth (stim, r_residue, ind_train, ind_validation, weight_ini, LR_f, n_iter_forward, lamda);

            r_residue=r_residue-r_pred(:,ii);
            amp(ii)=std(r_pred(:,ii));

            % 	subplot(ceil(n_filter_forward/5),5, ii,'Parent',h)
            % 	imagesc(reshape(filters(:,ii),yr,xr))
            %     axis equal
            %     axis off

        end


        %% %%%%%%%%%%%%% backward step %%%%%%%%%%%%%%%%

        data_bw=struct([]); % backward data storage

        n_elim = n_filter_forward-1;
        for ii=1:n_elim

            % sort
        	[~, seq]=sort(amp, 'descend');
        	filters=filters(:,seq);
            CRF_par=CRF_par(:,seq);
            amp=amp(seq);
            r_pred=r_pred(:,seq);

            n_filter=length(amp);

            % cumulative validation cc
            r_pred_cum_validation=zeros(nvalidation,1);
            validation_cc_cum=zeros(n_filter,1);
            for jj=1:n_filter
                r_pred_cum_validation=r_pred_cum_validation+r_pred(ind_validation,jj);
                tmp=corrcoef(r_pred_cum_validation, resp(ind_validation));
                validation_cc_cum(jj)=tmp(1,2);
            end

            data_bw(n_filter).filters=filters;
            data_bw(n_filter).CRF_par=CRF_par;
            data_bw(n_filter).amp=amp;
            data_bw(n_filter).validation_cc_cum=validation_cc_cum;


            % eliminate last one
            filters=filters(:,1:end-1);
        	CRF_par=CRF_par(:,1:end-1);
        	amp=amp(1:end-1);
        	r_pred=r_pred(:,1:end-1);

            n_filter=n_filter-1;

            for jj=1:n_filter

                ind_rest=true(n_filter,1);
                ind_rest(jj)=false;
                r_residual=resp-sum( r_pred(:,ind_rest), 2 );

                [filters(:,jj), CRF_par(:,jj), r_pred(:,jj)] = ...
                    weight_fit_2order_4_err_smooth (stim, r_residual, ind_train, ind_validation, filters(:,jj), LR_b, n_iter_backward, lamda);

                amp(jj)=std(r_pred(:,jj));

            end

            %     figure('Position', [57 875 1478 97]);colormap(gray)
            %     for  kk=1:length(amp)
            %         subplot	(1, n_filter, kk)
            %         imagesc( reshape(filters(:,kk),yr,xr) )
            %         axis equal
            %         axis off
            %     end

        end

        % validation cc for last filter  (n_filter=1)
        validation_cc_cum=zeros(n_filter,1);
        tmp=corrcoef(r_pred(ind_validation,jj), resp(ind_validation));
        validation_cc_cum(jj)=tmp(1,2);

        data_bw(n_filter).filters=filters;
        data_bw(n_filter).CRF_par=CRF_par;
        data_bw(n_filter).amp=amp;
        data_bw(n_filter).validation_cc_cum=validation_cc_cum;


        output_name=['PPR4_result_2nd_err_' cell_id '[' num2str(N) '-' num2str(Npart) '].' fndate '.mat'];
        save( fullfile(outputpath, output_name), 'cell_id', ...
            'Part','N','Npart', ...
            'polynomial_order', 'analysis_patch', ...
            'data_bw', ...
            'n_filter_forward','n_iter_forward','n_iter_backward','LR_f','LR_b', ...
            'lamda' );


    end
end
