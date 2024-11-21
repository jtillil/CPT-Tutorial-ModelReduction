%% test file for redmodel_population_error

model = simulate_and_reduce_numqss('in_vitro_PTtest_highTF');
I = model.I;


%% PT test low TF
redmodel_PTlow=struct;
error_PTlow=struct;
%save('results/sampling_redmodel_PTlow.mat','-struct','redmodel_PTlow','-v7.3')
%save('results/sampling_error_PTlow.mat','-struct','error_PTlow','-v7.3')

parfor seed=1:100
    for sample_size=[20,50,100,200]%[500,1000]%[5,10,20,50,100,200]%,200]%[20,50,100,200,500,1000]
        for relerrTOL=[0.1]%[0.01,0.1]%[0.01,0.05,0.1,0.2]
            for sampling=["naive","lhs","lhs_eog"]%
                 fprintf(['working on reduction for ',char(sampling),'_',int2str(sample_size),'_',int2str(relerrTOL*100),'_',int2str(seed)]);
                 redmodel_PTlow.([char(sampling),'_',int2str(sample_size),'_',int2str(relerrTOL*100),'_',int2str(seed)])=...
                     simulate_and_reduce_numqss("in_vitro_PTtest_lowTF",'sample_size', sample_size,...
                     'share',0.95,'force_lhs',true,'reduce',true,'relerrTOL',relerrTOL,'sampling',sampling,...
                     'output_time',(0:0.5:360)/3600,'seed',seed,'only_random',true,'noqss',true);
                 error_PTlow.([char(sampling),'_',int2str(sample_size),'_',int2str(relerrTOL*100),'_',int2str(seed)])=...
                     redmodel_population_error(redmodel_PTlow.([char(sampling),'_',int2str(sample_size),'_',int2str(relerrTOL*100),'_',int2str(seed)]).I,...
                     'in_vitro_PTtest_lowTF',(0:0.5:360)/3600);
            end
            dyn_struct=structfun(@(x) I.nmstate(x.I.dyn),redmodel_PTlow,'UniformOutput',false);
            mySave('results/sampling_dyn_PTlow.mat',dyn_struct)
            env_struct=structfun(@(x) I.nmstate(x.I.env),redmodel_PTlow,'UniformOutput',false);
            mySave('results/sampling_env_PTlow.mat',env_struct)
            mySave('results/sampling_error_PTlow.mat',error_PTlow)
            redmodel_PTlow=struct;
            error_PTlow=struct;
            dyn_struct=struct;
            env_struct=struct;
        end
    end
end

function mySave(filenm, variable)
    save(filenm,'-struct','variable','-append');
end

%% for evaluation:
%structfun(@(x) quantile(x,0.95),error_PTlow)
%
%structfun(@(x) x.Y_refs'*3600,redmodel_PTlow,'UniformOutput',false)