function nested = main_functions_3D_b_4
    % Return handles to subfunctions for 3D 3D B-spline registration
    nested.BSplineTransformation3D = @BSplineTransformation3D;
    nested.ImageRegistration3D = @ImageRegistration3D;
    nested.applyTransformation3D = @applyTransformation3D;
    
    %% ----------------------------
    function T = BSplineTransformation3D(x, beta, k, z)
        % x: 3 x N points (columns)
        % beta: 3 * prod(z) x 1, control point displacements stored as
        %       [cp1_x; cp1_y; cp1_z; cp2_x; cp2_y; cp2_z; ...]
        % k: control point spacing [kx,ky,kz]
        % z: number of control points [z1,z2,z3]
    
        N = size(x,2);
        T = zeros(3,N);
    
        for idx = 1:N
            pt = x(:,idx)';        % [i j l]
            def = zeros(1,3);
    
            u = (pt - 1) ./ k;    % map voxel -> control-point coords (zero-based)
            base = floor(u);
            f = u - base;
    
            Bx = computeBSplineWeights(f(1));
            By = computeBSplineWeights(f(2));
            Bz = computeBSplineWeights(f(3));
    
            % iterate over 4x4x4 neighborhood
            for ii = 0:3
                for jj = 0:3
                    for kk = 0:3
                        cx = base(1) + ii; cy = base(2) + jj; cz = base(3) + kk;
                        if cx < 0 || cy < 0 || cz < 0 || cx >= z(1) || cy >= z(2) || cz >= z(3)
                            continue;
                        end
    
                        % convert to 1-based subscript for sub2ind
                        sub_i = cx + 1; sub_j = cy + 1; sub_k = cz + 1;
                        cpIndex = sub2ind(z, sub_i, sub_j, sub_k);     % 1..prod(z)
                        idx_beta = (cpIndex-1)*3;                       % start index in beta (0-based block)
                        disp_vec = beta(idx_beta + (1:3));
    
                        w = Bx(ii+1) * By(jj+1) * Bz(kk+1);
                        def = def + w * disp_vec(:)';
                    end
                end
            end
            T(:,idx) = def';
        end
    end

    %% ----------------------------
    function [beta, Image, Dssd, U, V, W] = ImageRegistration3D(R, T, k)
        % 3D Image Registration using Gauss-Newton + line search

        % Validate input
        if length(k) ~= 3
            error('k must have 3 elements');
        end
        [n1,n2,n3] = size(R);
        [t1,t2,t3] = size(T);
        if any([n1 n2 n3] ~= [t1 t2 t3])
            error('R and T must have the same size');
        end

        % Control points
        z = ceil([n1 n2 n3]./k) + 3;
        lG = prod(z);

        % Initialize
        beta = zeros(3*lG,1);
        lambda = 20; maxIter = 15; p = 0.5;
        alphaSteps = [1 0.5 0.25 0.1];
        Dssd = zeros(maxIter,1);

        % Compute image gradients
        [dxT, dyT, dzT] = imgradientxyz(T,'central');

        % Iterative optimization
        for iter = 1:maxIter
            [J,f] = Gradient3D(R, T, dxT, dyT, dzT, beta, k, z);
            JJ = J'*J + lambda*speye(size(J,2));
            s = -JJ \ (J'*f);
            beta = lineSearch3D(R,T,beta,s,k,z,alphaSteps,p);
            Dssd(iter) = DSSD3D(R,T,beta,k,z);
        end

        % Apply final transformation
        [Image,U,V,W] = applyTransformation3D(R,T,beta,k,z);
        Dssd = Dssd(1:iter);
    end

    %% ----------------------------
    function out = computeBSplineWeights(u)
        % Cubic B-spline weights for fractional position u
        out = [((1-u)^3)/6, (3*u^3 - 6*u^2 + 4)/6, ...
               (-3*u^3 + 3*u^2 + 3*u +1)/6, (u^3)/6];
    end

    %% ----------------------------
    function idxs = getControlPointIndices(pt, k, z)
        % Returns linear indices (1-based) of control points in 4x4x4 neighborhood
        u = (pt - 1) ./ k;
        base = floor(u);
        idxs = [];
        for ii = 0:3
            for jj = 0:3
                for kk = 0:3
                    cx = base(1) + ii; cy = base(2) + jj; cz = base(3) + kk;
                    if all([cx cy cz] >= 0) && cx < z(1) && cy < z(2) && cz < z(3)
                        sub_i = cx + 1; sub_j = cy + 1; sub_k = cz + 1;
                        idx = sub2ind(z, sub_i, sub_j, sub_k);
                        idxs = [idxs idx];
                    end
                end
            end
        end
    end

    %% ----------------------------
    function w = computeWeight(voxel, cp_idx, k, z)
        % cp_idx is linear index (1-based)
        [cp_i, cp_j, cp_k] = ind2sub(z, cp_idx);
    
        % distance in units of spacing
        u = (voxel(1) - (cp_i-1)*k(1)) / k(1);
        v = (voxel(2) - (cp_j-1)*k(2)) / k(2);
        w_ = (voxel(3) - (cp_k-1)*k(3)) / k(3);
    
        % do NOT clamp here; cubicBSpline expects distances (can be up to 2)
        w = cubicBSpline(u) * cubicBSpline(v) * cubicBSpline(w_);
    end
    
    function B = cubicBSpline(t)
        t = abs(t);
        if t < 1
            B = (3/2)*t^3 - (5/2)*t^2 + 1;
        elseif t < 2
            B = (-1/2)*t^3 + (5/2)*t^2 - 4*t + 2;
        else
            B = 0;
        end
    end


    %% ----------------------------
    function beta_new = lineSearch3D(R,T,beta,s,k,z,alphaSteps,p)
        % Line search to update beta
        alpha = 1.0; best_beta = beta;
        min_cost = DSSD3D(R,T,beta,k,z);
        found = false;
        for i=1:length(alphaSteps)
            step = alphaSteps(i);
            beta_test = beta + step*s;
            cost = DSSD3D(R,T,beta_test,k,z);
            if cost<min_cost
                min_cost=cost; best_beta=beta_test; found=true; break;
            end
        end
        if ~found
            for i=1:5
                beta_test = beta + alpha*s;
                cost = DSSD3D(R,T,beta_test,k,z);
                if cost<min_cost, best_beta=beta_test; break; end
                alpha=alpha*p;
            end
        end
        beta_new = best_beta;
    end

    %% ----------------------------
    function [Image,U,V,W] = applyTransformation3D(R,T,beta,k,z)
        [m,n,p] = size(R);
        Image = zeros(m,n,p);
        U = zeros(m,n,p); V = zeros(m,n,p); W = zeros(m,n,p);
        for i=1:m
            for j=1:n
                for l=1:p
                    disp_vec = BSplineTransformation3D([i;j;l], beta, k, z);
                    U(i,j,l)=disp_vec(1); V(i,j,l)=disp_vec(2); W(i,j,l)=disp_vec(3);
                    Image(i,j,l) = TrilinearApp(T,[i-disp_vec(1), j-disp_vec(2), l-disp_vec(3)]);
                end
            end
        end
    end

    %% ----------------------------
    function cost = DSSD3D(R,T,beta,k,z)
        [transformed_T,~,~,~] = applyTransformation3D(R,T,beta,k,z);
        cost = sum((R(:)-transformed_T(:)).^2);
    end

    %% ----------------------------
    function val = TrilinearApp(T, x)
        [m,n,p] = size(T);
        x1 = floor(x(1)); x2 = floor(x(2)); x3 = floor(x(3));
        x1h = ceil(x(1)); x2h = ceil(x(2)); x3h = ceil(x(3));
        if x(1)<1||x(1)>m||x(2)<1||x(2)>n||x(3)<1||x(3)>p, val=0; return; end
        a000 = safeGet3D(T,x1,x2,x3);
        a100 = safeGet3D(T,x1h,x2,x3);
        a010 = safeGet3D(T,x1,x2h,x3);
        a001 = safeGet3D(T,x1,x2,x3h);
        a110 = safeGet3D(T,x1h,x2h,x3);
        a101 = safeGet3D(T,x1h,x2,x3h);
        a011 = safeGet3D(T,x1,x2h,x3h);
        a111 = safeGet3D(T,x1h,x2h,x3h);
        xd = x(1)-x1; yd=x(2)-x2; zd=x(3)-x3;
        c00=a000*(1-xd)+a100*xd; c01=a001*(1-xd)+a101*xd;
        c10=a010*(1-xd)+a110*xd; c11=a011*(1-xd)+a111*xd;
        c0=c00*(1-yd)+c10*yd; c1=c01*(1-yd)+c11*yd;
        val = c0*(1-zd)+c1*zd;
    end

    function val = safeGet3D(T,x,y,z)
        [m,n,p] = size(T);
        if x<1||x>m||y<1||y>n||z<1||z>p, val=0; else val=T(x,y,z); end
    end

    %% ----------------------------
    function [J,f] = Gradient3D(R,T,dxT,dyT,dzT,beta,k,z)
        [m,n,p] = size(R);
        numCP = prod(z);
        J = sparse(m*n*p, 3*numCP);
        f = zeros(m*n*p,1);
        index = 1;
        for i=1:m
            for j=1:n
                for kIndex=1:p
                    disp_vec = BSplineTransformation3D([i;j;kIndex],beta,k,z);
                    pt = [i-disp_vec(1), j-disp_vec(2), kIndex-disp_vec(3)];
                    intensityDiff = TrilinearApp(T,pt) - R(i,j,kIndex);
                    f(index) = intensityDiff;
                    gradX = TrilinearApp(dxT,pt);
                    gradY = TrilinearApp(dyT,pt);
                    gradZ = TrilinearApp(dzT,pt);
                    cpIdxs = getControlPointIndices([i,j,kIndex], k, z);
                    for c=1:length(cpIdxs)
                        cp = cpIdxs(c);
                        w = computeWeight([i,j,kIndex],cp,k,z);
                        J(index,cp) = w*gradX;
                        J(index,cp+numCP) = w*gradY;
                        J(index,cp+2*numCP) = w*gradZ;
                    end
                    index = index + 1;
                end
            end
        end
    end
end
