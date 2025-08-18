function [ model_best ] = find_the_best_model ( models, stim, resp )
% FIND_BEST_MODEL
% Find the model with significant highest cc with least subunits.
% 1. find the model with maximum mean of cc
% 2. the best model should has ccs not significant different with the
%     maximum one and having least number of subunits.
%
%  ignore NaN in ccs 

% verify models are from the same cell
cell_id_all = {models.PPR4.cell_id};
cell_id = cell_id_all{1};
assert(all(strcmp(cell_id, cell_id_all)));

Tp=0.01; % p threshold

n=length(models.PPR4);

ccsall={models.PPR4.ccs};
ccsall=cell2mat(ccsall);
mccsall=nanmean(ccsall);
[~, max_ind]=max(mccsall);

max_ccs=ccsall(:,max_ind);
% fprintf('%e\n',models.PPR4(max_ind).lamda);
% max_ccs=max_ccs(~isnan(max_ccs));
max_n=models.PPR4(max_ind).n; % number of subunits
max_lamda=models.PPR4(max_ind).lamda ;

model_best=models.PPR4(max_ind);
model_best.index=max_ind;

for i=1:n
    current_ccs=models.PPR4(i).ccs;
%     current_ccs=current_ccs(~isnan(current_ccs));
    p=signrank(max_ccs,  current_ccs);
    
    if p>Tp
        
        if model_best.n > models.PPR4(i).n
            model_best=models.PPR4(i);
            model_best.index=i;
            
        elseif model_best.n == models.PPR4(i).n
            ccbest=nanmean(model_best.ccs);
            cc = nanmean(models.PPR4(i).ccs);
            if cc > ccbest
                model_best=models.PPR4(i);
                model_best.index=i;
            end

        end
        
    end
    
end

% estimate nonlinear function
model_best.nonlin = struct('x',[],'y',[],'e',[],'ns',[],'rlin_std',[]);
[yres, xres, frames]=size(stim);
model_best.filters = model_best.filters - mean(model_best.filters);
nfilters = size(model_best.filters, 2);
for j=1:nfilters
    rlin=model_best.filters(:,j)'*reshape(stim,xres*yres,frames);
    rlin = rlin/std(rlin);
    model_best.nonlin(j) = Compute_nonlin_func(rlin, resp);
end

model_best.cell_id = cell_id;
model_best.valid_cc = nanmean(model_best.ccs);

model_best.stimpos = [];
model_best.max_lamda = max_lamda;

end

