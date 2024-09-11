function results_str = ...
    sketchandsolveLSQR(A,b,s,uqr,rndseed,lsqrtol,lsqritmax,xtrue,rtrue,...
    xtruen,rtruen,results_str,ind)

    us = results_str.us;
    
    % generate the preconditioner and sketch-and-solve solution
    [R,x0] = sketchgauss(A,b,s,us,uqr,rndseed);
    
    results_str.R = R;
   
     % compute norms and condition number
     ARinv = double(A/R);
     [m,n] = size(A);
     results_str.cond_aug_pr(ind) = cond([eye(m) ARinv; ARinv' zeros(n)]);
     results_str.norm_ARinv(ind) = norm(ARinv);
     results_str.norm_pinv_ARinv(ind) = 1/svds(ARinv,1,'smallest');

     % initial solve via LSQR
    [x,~,~,lsqrit] = lsqr(A,b,lsqrtol,lsqritmax,R,[],x0);

    results_str.lsqrit(ind) = lsqrit;

    r = b - A*x;

    results_str.xlsqr = x;
    results_str.rlsqr = r;

     % error in x and r after the initial solve
     results_str.x0err(ind) = norm(xtrue - x)/xtruen;
     results_str.r0err(ind) = norm(rtrue - r)/rtruen;
