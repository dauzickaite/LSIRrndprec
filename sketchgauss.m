function [R,x] = sketchgauss(A,b,s,us,uqr,rndseed)
% half precision simulated in matrx-matrix product via chop
rng(rndseed)
[m,n] = size(A);

% precision us
switch us
    case {'half','single'}
        A = single(A);
        b = single(b);
        G = (1/sqrt(s))*single(randn(s,m));
    case 'double'
        A = double(A);
        b = double(b);
        G = (1/sqrt(s))*randn(s,m);
end

% skecth A and b in us
switch us 
    case 'half' 
        % scaling
        G = (1/sqrt(s))*randn(s,m);
        opt.format = 'h';
        chop([],opt)
        Y = zeros(s,n);
        for i = 1:m
            Y = chop(Y + chop(G(:,i)*A(i,:)));
        end
        Gb = zeros(s,1);
        b = chop(b);
        for i = 1:m
            Gb = chop(Gb+chop(b(i,1)*G(:,i)));
        end
    case {'single','double'}
        Y = G*A;
        Gb = G*b;
end

% compute QR and x in uqr
switch uqr
    case {'half','single'}
        Y = single(Y);
        Gb = single(Gb);
    case 'double'
        Y = double(Y);
        Gb = double(Gb);
end


switch uqr
    case 'half' 
        addpath '../Multi_precision_NLA_kernels-master'
        [Q,R] = house_qr_lp(Y,0); 
        QTGb = zeros(s,1);
        Gb = chop(Gb);
        for i = 1:m
            QTGb = chop(QTGb+chop(Gb(i,1)*Q(i,:)'));
        end
        x = trisol(R,QTGb);
    case {'single','double'} 
        [Q,R] = qr(Y,"econ");
        x = R\(Q'*Gb);
end


