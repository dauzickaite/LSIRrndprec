clear
rng(1234)
m = 4*1e2;
n = 10; 

cond_no_list = 10.^(0:16);
cond_no_list_length = length(cond_no_list); 

normARsinv = zeros(cond_no_list_length,1);
normARhinv = zeros(cond_no_list_length,1);


for ind = 1:cond_no_list_length
    % generate A in double with different condition numbers
    cond_no = cond_no_list(ind);
    A = gallery('randsvd',[m,n],cond_no,3);

    % cast A to single As and half Ah
    As = single(A);
    Ah = half(A);

    % compute R factors of As and Ah in double
    [~,Rs] = qr(double(As),0);
    [~,Rh] = qr(double(Ah),0);

     % compute norms of pinv(A/Rs) and pinv(A/Rh)
    svARs = svd(A/Rs);
    svARh = svd(A/Rh);

    normARsinv(ind) = 1/min(svARs);
    normARhinv(ind) = 1/min(svARh);

end


%% plot the norms and the bounds
figure; semilogy(normARsinv,'ko','LineWidth',8,'MarkerSize',20); hold on
semilogy(normARhinv,'rx','LineWidth',8,'MarkerSize',20); hold on

semilogy((1+ sqrt(n)*2^(-24)*cond_no_list),'k--','LineWidth',8); hold on
semilogy((1+ sqrt(n)*2^(-11)*cond_no_list),'r:','LineWidth',8); hold on


legend('||  (A R^{-1}_s)^+ ||_2','|| (A R^{-1}_h)^+ ||_2', '1+n^{1/2}u_s \kappa(A)', '1+n^{1/2} u_h \kappa(A)')
legend('Location','northwest')
xticks(1:2:17)
xticklabels({'1e0','1e2','1e4','1e6','1e8','1e10','1e12','1e14','1e16'})
xlabel('\kappa(A)')
ylim([1e-1,1e14])
yticks([1e0,1e4,1e8,1e12])
xlim([0.5,17.5])

fontsize(gcf, 40, 'points')

