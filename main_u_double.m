clear
addpath '../AdvanpixMCT-4.8.5.14607'
addpath '../chop-master'

m = 1e3;
n = 1e2;
s=4*n;

ur = 'quad';
iritmax = 30;
irtol = 2*eps('double');
lsqritmax = 2*n;
lsqrtol = 1e-12;
gmresitmax = 50;
gmrestol = 1e-12;

kappa_list = [1e2,1e4,1e6,1e8,1e10,1e12,1e14,1e15];
kappa_no = length(kappa_list);

rng(112358)
b = rand(m,1);
b = b/norm(b);

% fgmres precisions
precisions_fgmres = struct('u','double','uA','double','upr','double');
precisions_fgmres_increased = struct('u','double','uA','quad','upr','quad');
% precision for computing QR of a sketched matrix
uqr = 'double';

cond_aug = zeros(kappa_no,1);
xtrue_norm = zeros(kappa_no,1);
rtrue_norm = zeros(kappa_no,1);

variants_no = 3;

results(1:variants_no) = struct('R', zeros(n,n), 'xlsqr', zeros(n,1), ...
    'rlsqr',zeros(m,1),'cond_aug_pr',zeros(kappa_no,1),...
    'norm_ARinv', zeros(kappa_no,1),'norm_pinv_ARinv', zeros(kappa_no,1),...
    'lsqrit',zeros(kappa_no,1), 'x0err', zeros(kappa_no,1), ...
    'r0err', zeros(kappa_no,1), 'lsirit', zeros(kappa_no,1),...
    'ir_convergence',false,'gmresit', zeros(kappa_no,1),...
    'x_error_ir', zeros(kappa_no,iritmax),'r_error_ir', zeros(kappa_no,iritmax),...
    'increased_fgmres_prec', false(kappa_no,1), 'lsqrflag',zeros(kappa_no,1));

results(1).us = 'half';
results(2).us = 'single';
results(3).us = 'double';

rndseed = 98765;
for ind = 1:kappa_no
    kappa = kappa_list(ind);

    % generate the data and true x and r
    rng(13213455)
    A = gallery('randsvd',[m,n],kappa,3);

    cond_aug(ind) = cond([eye(m) A; A' zeros(n)]);

    xtrue = mp(A,64)\mp(b,64);
    xtruen = norm(xtrue);
    rtrue = mp(b,64) - mp(A,64)*mp(xtrue,64);
    rtruen = norm(mp(rtrue,64));

    xtrue_norm(ind) = xtruen;
    rtrue_norm(ind) = rtruen;

    % generate the preconditioner and solve via LSQR
    results(1) = ...
        sketchandsolveLSQR(A,b,s,uqr,rndseed,lsqrtol,lsqritmax,...
        xtrue,rtrue,xtruen,rtruen,results(1),ind);


    results(2) = ...
        sketchandsolveLSQR(A,b,s,uqr,rndseed,lsqrtol,lsqritmax,...
        xtrue,rtrue,xtruen,rtruen,results(2),ind);
    

    results(3) = ...
        sketchandsolveLSQR(A,b,s,uqr,rndseed,lsqrtol,lsqritmax,...
        xtrue,rtrue,xtruen,rtruen,results(3),ind);

    % run IR; if no converge, increase uA and upr in FGMRES and run IR
    % again
        if kappa <= 1e4        
            results(1) = ...
                lsirfgmres(A,b,xtrue,rtrue,xtruen,rtruen,...
                iritmax,irtol,gmrestol,gmresitmax,results(1),ind,precisions_fgmres,ur);
        end

        if kappa <= 1e6
            results(2) = ...
                lsirfgmres(A,b,xtrue,rtrue,xtruen,rtruen,...
                iritmax,irtol,gmrestol,gmresitmax,results(2),ind,precisions_fgmres,ur);
        end

        if kappa <= 1e10
            results(3) = ...
            lsirfgmres(A,b,xtrue,rtrue,xtruen,rtruen,...
            iritmax,irtol,gmrestol,gmresitmax,results(3),ind,precisions_fgmres,ur);
        else
            results(3).increased_fgmres_prec(ind) = true;
            gmresitmax = 80;    
            results(3) = ...
             lsirfgmres(A,b,xtrue,rtrue,xtruen,rtruen,...
             iritmax,irtol,gmrestol,gmresitmax,results(3),ind,...
             precisions_fgmres_increased,ur);
        end

end

