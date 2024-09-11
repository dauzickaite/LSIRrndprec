function results_str = lsirfgmres(A,b,xtrue,rtrue,xtruen,rtruen,...
    ir_it_max,irtol,gmrestol,gmresitmax,results_str,indkappa,...
    precisions_fgmres,ur)
    
    % fgmres precisions
    u = precisions_fgmres.u;
    uA = precisions_fgmres.uA;
    upr = precisions_fgmres.upr;
    
    [m,n] = size(A);
    Aa = [eye(m) A; A' zeros(n)];

    R = results_str.R;
    Pr = @(x, convprc) Pblock_R(x,convprc,R,m,u,upr);
    Pl = @(x, convprc) Pblock_R(x,convprc,R',m,u,upr);
    
    % check if we need to run IR
    converged = false;
    if results_str.x0err(indkappa) <= irtol && results_str.r0err(indkappa) <= irtol
        converged = true;
    end
    ind = 0;
    gmres_it_total = 0;
    
    x_relerror = zeros(ir_it_max,1);
    r_relerror = zeros(ir_it_max,1);

    x = results_str.xlsqr;
    r = results_str.rlsqr;
    
    while ~converged && ind < ir_it_max
        ind = ind+1;
    
        % compute the residuals for augmented system
        switch ur
            case 'double'                
                f = double(b) - double(A)*double(x) - double(r);
                g = -double(A)'*double(r);
            case 'quad'
                f = mp(b,34) - mp(A,34)*mp(x,34) - mp(r,34);
                g = -mp(A,34)'*mp(r,34);
        end
        
        switch u
            case 'single'
                % solve
                [z,~,~,git] = mpgmres(Aa,single([f;g]),zeros(m+n,1),gmrestol,1,...
                    gmresitmax,Pr,Pl,u,uA);
                        
                % update
                x = single(x) + single(z(m+1:end));
                r = single(r) + single(z(1:m));

            case 'double'
                % solve
                [z,~,~,git] = mpgmres(Aa,double([f;g]),zeros(m+n,1),gmrestol,1,...
                    gmresitmax,Pr,Pl,u,uA);
                                            
                % update
                x = double(x) + double(z(m+1:end));
                r = double(r) + double(z(1:m));
        end
        
        gmres_it_total = gmres_it_total + git;

        % compute the error 
        r_relerror(ind) = norm(mp(r,64) - mp(rtrue,64))/rtruen;
        x_relerror(ind) = norm(mp(x,64) - mp(xtrue,64))/xtruen;
        
        % check for convergence
        if x_relerror(ind) <= irtol && r_relerror(ind) <= irtol
            converged = true;
        end
        
    end
    
    x_relerror(ind+1:end) = nan;
    r_relerror(ind+1:end) = nan;

    results_str.lsirit(indkappa) = ind;
    results_str.gmresit(indkappa) = gmres_it_total;
    results_str.ir_convergence(indkappa) = converged;
    results_str.x_error_ir(indkappa,:) = x_relerror;
    results_str.r_error_ir(indkappa,:) = r_relerror;

