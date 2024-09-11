clear
rng(1234)
m = 4*1e2;
n = 10; 

s = 4*n; % sketching parameter
uqr = 'double';
rndseed = 20242024;
b = zeros(m,1);

cond_no_list = 10.^(0:16);
cond_no_list_length = length(cond_no_list); 

cond_preconditioned_augmented = zeros(cond_no_list_length,1);
cond_preconditioned_augmented_single = zeros(cond_no_list_length,1);
boundARinv_cond = zeros(cond_no_list_length,1);
boundARinv_norms = zeros(cond_no_list_length,1);
boundARinv_cond_single = zeros(cond_no_list_length,1);
boundARinv_norms_single = zeros(cond_no_list_length,1);

% generate A in double with different condition numbers 1e1 ro 1e16
for ind = 1:cond_no_list_length
    cond_no = cond_no_list(ind);
    A = gallery('randsvd',[m,n],cond_no,3);

    % generate R with sketching in double and single
    R = sketchgauss(A,b,s,'double',uqr,rndseed);
    Rs = sketchgauss(A,b,s,'single',uqr,rndseed);

    % compute the preconditioned matrices and their svd
    AR = A/R;
    ARs = A/Rs;

    svAR = svd(AR);
    svARs = svd(ARs);

    % compute the condition numbers
    cond_preconditioned_augmented(ind) = cond([eye(m) AR; AR' zeros(n)]);
    cond_preconditioned_augmented_single(ind) = cond([eye(m) ARs; ARs' zeros(n)]);

    % compute the bounds
    boundARinv_cond(ind) = 2*svAR(1)/svAR(end);
    boundARinv_cond_single(ind) = 2*svARs(1)/svARs(end);

    if svAR(end) >= sqrt(2)
        boundARinv_norms(ind) = 1 + svAR(1);
    else
        boundARinv_norms(ind) = 2*(1 + svAR(1))/(sqrt(1 + 4*svAR(end)^2) - 1);
    end

    if svARs(end) >= sqrt(2)
        boundARinv_norms_single(ind) = 1 + svARs(1);
    else
        boundARinv_norms_single(ind) = 2*(1 + svARs(1))/(sqrt(1 + 4*svARs(end)^2) - 1);
    end
     

end


%% plot the norms and the bounds
figure; semilogy(cond_preconditioned_augmented,'k','LineWidth',10); hold on

semilogy(boundARinv_cond,'b--','LineWidth',10); hold on
semilogy(boundARinv_norms,'r:','LineWidth',10); hold on

legend('$\kappa(M_L^{-1} \tilde{A} M_R^{-1})$','$2 \kappa(A \hat{R}^{-1})$', 'gen. bound', 'Interpreter', 'LaTeX')
legend('Location','northwest')
xticks(1:2:17)
xticklabels({'1e0','1e2','1e4','1e6','1e8','1e10','1e12','1e14','1e16'})
xlabel('\kappa(A)')
xlim([0.5,17.5])
fontsize(gcf, 50, 'points')

% sketching in single
figure; semilogy(cond_preconditioned_augmented_single,'k','LineWidth',10); hold on

semilogy(boundARinv_cond_single,'b--','LineWidth',10); hold on
semilogy(boundARinv_norms_single,'r:','LineWidth',10); hold on

xticks(1:2:17)
xticklabels({'1e0','1e2','1e4','1e6','1e8','1e10','1e12','1e14','1e16'})
xlabel('\kappa(A)')
ylim([1e-1,1e17])
yticks([1e0,1e4,1e8,1e12,1e16])
xlim([0.5,17.5])
fontsize(gcf, 50, 'points')



