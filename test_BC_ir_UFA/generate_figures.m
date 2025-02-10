%% warfarin article figures

%% Figure 4 loglog-plot
load('results/loglogplot.mat','ref_model')
factor_product = permute(prod(permute(ref_model.X_refs(:,:,[I.II,I.VII,I.X]),[3,2,1])),[3,2,1]);
figure;
scatter(factor_product(:,[2,6,12,31])./factor_product(:,1),ref_model.Y_refs(:,[2,6,12,31])./ref_model.Y_refs(:,1),'ob')
hold on;
plot([1,10^-4],[1,(10^-4)^-0.1975])
xlabel('normalised II*VII*X')
ylabel('INR')
set(gca, 'yScale', 'log')
set(gca, 'xScale', 'log')

%% Figure 5a approximation in population
load('results/population_approximation.mat','ref_model','red_model')
factor_product_red = permute(prod(permute(red_model.X_reds(:,:,[I.II,I.VII,I.X]),[3,2,1])),[3,2,1]);
figure;
scatter(ref_model.Y_refs./ref_model.Y_refs(:,1),(factor_product_red./factor_product_red(:,1)).^-0.1975,'ob','MarkerEdgeAlpha',.1)


%% Figure 5b different CYP2C9 genotypes
load('results/cyp_simulations.mat','model11','model12','model22','model13','model23','model33','model33_doseopt')
figure(6); hold on
plot(model11.t_ref,model11.Y_ref)
plot(model12.t_ref,model12.Y_ref)
plot(model22.t_ref,model22.Y_ref)
plot(model13.t_ref,model13.Y_ref)
plot(model23.t_ref,model23.Y_ref)
plot(model33.t_ref,model33.Y_ref)
plot(model33.t_ref,model33_doseopt.Y_ref)
plot(model11.t_ref,model11.Y_red,'--')
plot(model12.t_ref,model12.Y_red,'--')
plot(model22.t_ref,model22.Y_red,'--')
plot(model13.t_ref,model13.Y_red,'--')
plot(model23.t_ref,model23.Y_red,'--')
plot(model33.t_ref,model33.Y_red,'--')
plot(model33.t_ref,model33_doseopt.Y_red,'--')
makeplotbold(6)
legend("*1/*1","*1/*2","*2/*2","*1/*3","*2/*3","*3/*3","*3/*3 reduced dose")
xlabel('time [h]')
ylabel('INR')
title('QSP model simulation vs. reduced model simulation for different CYP2C9 genotypes')

%% Figure 5c different VKORC1 genotypes
load('results/vko_simulations.mat','modelv11','modelv12','modelv22')
figure(7); hold on
col=lines(7);
plot(modelv11.t_ref,modelv11.Y_ref,'Color',col(1,:))
plot(modelv12.t_ref,modelv12.Y_ref,'Color',col(2,:))
plot(modelv22.t_ref,modelv22.Y_ref,'Color',col(3,:))
plot(modelv11.t_ref,modelv11.Y_red,'--','Color',col(1,:))
plot(modelv11.t_ref,modelv12.Y_red,'--','Color',col(2,:))
plot(modelv12.t_ref,modelv22.Y_red,'--','Color',col(3,:))
%makeplotbold(7)
legend("GG","GA","AA")
xlabel('time [h]')
ylabel('INR')
title('QSP model simulation vs. reduced model simulation for different VKORC1 genotypes')


%% Figure Supplement alternative VKORC1 genotype parametrisation
load('results/vko_simulations_alternative.mat','model11','model12','model22')
figure(9)
plot(model11.t_ref,model11.Y_ref)
hold on
plot(model11.t_ref,model12.Y_ref)
plot(model11.t_ref,model22.Y_ref)
legend('GG','GA','AA')
xlabel('time [h]')
ylabel('INR')
