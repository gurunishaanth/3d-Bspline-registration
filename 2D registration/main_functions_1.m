function nested = main_functions_1()
    % Return handles to subfunctions
    nested.BSplineTransformation = @BSplineTransformation;
    nested.ImageRegistration = @ImageRegistration;

    function [out] = BSplineTransformation(x, beta, k, z)
        % Optimized B-Spline Transformation
        k1 = k(1);
        k2 = k(2);
        z2 = z(2);

        % Calculate grid indices
        y1 = floor((x(1) - 1) / k1) * z2 * 2;
        y2 = floor((x(2) - 1) / k2) * 2 + 1;

        % Calculate u values
        u = (x(1) - 1) / k1 - floor((x(1) - 1) / k1);
        u2 = (x(2) - 1) / k2 - floor((x(2) - 1) / k2);

        % Precompute B-spline weights
        zz = computeBSplineWeights(u);
        zz2 = computeBSplineWeights(u2);

        % Precompute index offsets
        tBase = y2 + y1;
        offsets = [0, 2, 4, 6];
        zOffsets = 2 * z2 * (0:3);

        % Compute output
        out = zeros(2, 1);  % Initialize output as 2D vector
        for i = 1:4
            for j = 1:4
                idx = tBase + offsets(i) + offsets(j) + zOffsets((i - 1) + 1);
                maxIdx = numel(beta);
                if idx > maxIdx
                    idx = maxIdx;  % Limit index to valid range
                end
                if idx + 1 > maxIdx
                    out = out + beta(idx) * zz(i) * zz2(j);
                else
                    out = out + beta(idx:idx + 1) * zz(i) * zz2(j);
                end

            end
        end
    end

    function [beta, Image, Dssd] = ImageRegistration(R, T, k)
        % Optimized Image Registration

        % Validate input dimensions
        if length(k) ~= 2
            error('k must have 2 elements.');
        end

        % Image dimensions
        [n1, n2] = size(R);
        [t1, t2] = size(T);
        if n1 ~= t1 || n2 ~= t2
            error('R and T must have the same dimensions.');
        end

        % Control points calculation
        z = ceil([n1, n2] ./ k) + 3;
        lG = prod(z);

        % Initialize beta and parameters
        beta = zeros(2 * lG, 1);
        lambda = 2;
        maxIter = 5;
        p = 0.3;
        l = 0.7;
        alphaSteps = [1, 0.5, 0.1];
        Dssd = zeros(maxIter, 1);

        % Compute gradients of T
        [dxT, dyT] = imgradientxy(T, 'central');

        % Iterative optimization
        for iter = 1:maxIter
            iter
            [J, f] = Gradient(R, T, dxT, dyT, beta, k, z);
            JJ = J' * J + lambda * speye(size(J, 2));
            s = -JJ \ (J' * f);
            beta = lineSearch(R, T, beta, s, k, z, alphaSteps, p, l);
            Dssd(iter) = DSSD(R, T, beta, k, z);
        end

        % Apply transformation
        Image = applyTransformation(R, T, beta, k, z);
        Dssd = Dssd(1:iter);
    end

    function [out] = computeBSplineWeights(u)
        % Precompute B-spline weights
        out = [((1 - u)^3) / 6, ...
               (3 * u^3 - 6 * u^2 + 4) / 6, ...
               (-3 * u^3 + 3 * u^2 + 3 * u + 1) / 6, ...
               (u^3) / 6];
    end

    function beta = lineSearch(R, T, beta, s, k, z, alphaSteps, p, l)
        % Perform line search
        for alpha = alphaSteps
            betaNew = beta + alpha * s;
            if SimilarityMeasure(R, T, betaNew, k, z) < SimilarityMeasure(R, T, beta, k, z)
                beta = betaNew;
                break;
            end
        end
    end

    function Image = applyTransformation(R, T, beta, k, z)
        % Apply B-spline transformation
        [m, n] = size(R);
        Image = zeros(m, n);
        for i = 1:m
            for j = 1:n
                new_u = BSplineTransformation([i, j], beta, k, z);
                Image(i, j) = BilinearApp(T, [i - new_u(1), j - new_u(2)]);
            end
        end
    end

    function [out] = BilinearApp(T, x)
        % Optimized bilinear interpolation
        [m, n] = size(T);
        y11 = floor(x(1));
        y12 = floor(x(2));
        y21 = ceil(x(1));
        y22 = ceil(x(2));

        % Boundary checks
        if x(1) < 1 || x(1) > m || x(2) < 1 || x(2) > n
            out = 0;
            return;
        end

        % Retrieve neighboring values
        a11 = safeGet(T, y11, y12);
        a12 = safeGet(T, y11, y22);
        a21 = safeGet(T, y21, y12);
        a22 = safeGet(T, y21, y22);

        % Linear interpolation
        z1 = (y21 - x(1)) * a11 + (x(1) - y11) * a21;
        z2 = (y21 - x(1)) * a12 + (x(1) - y11) * a22;
        out = (y22 - x(2)) * z1 + (x(2) - y12) * z2;
    end

    function val = safeGet(T, x, y)
        % Safely retrieve value from T with boundary checks
        [m, n] = size(T);
        if x < 1 || x > m || y < 1 || y > n
            val = 0;
        else
            val = T(x, y);
        end
    end
    function [A, F] = Gradient(R, T, dxT, dyT, beta, k, z)
        % Optimized Gradient and Residual Computation for Image Registration
        k1 = k(1);
        k2 = k(2);
        z2 = z(2);
        [m, n] = size(T);
        [lG, ~] = size(beta);
    
        % Estimate maximum possible nonzero elements
        maxNonZeros = m * n * 32;
        I = zeros(maxNonZeros, 1);
        J = zeros(maxNonZeros, 1);
        V = zeros(maxNonZeros, 1);
        F = zeros(n * m, 1);
        Index = 0;
    
        for i = 1:m
            floor_i = floor((i - 1) / k1);
            t1 = floor_i * z2 * 2;
            u = (i - 1) / k1 - floor_i;
            
            % Precompute B-spline basis functions for x
            zz = [((1 - u)^3) / 6, (3 * u^3 - 6 * u^2 + 4) / 6, ...
                  (-3 * u^3 + 3 * u^2 + 3 * u + 1) / 6, (u^3) / 6];
    
            for j = 1:n
                floor_j = floor((j - 1) / k2);
                t2 = floor_j * 2 + 1;
                u2 = (j - 1) / k2 - floor_j;
                
                % Precompute B-spline basis functions for y
                zz2 = [((1 - u2)^3) / 6, (3 * u2^3 - 6 * u2^2 + 4) / 6, ...
                       (-3 * u2^3 + 3 * u2^2 + 3 * u2 + 1) / 6, (u2^3) / 6];
    
                % Compute transformation indices in a vectorized manner
                tt_base = t2 + t1;
                tt_offsets = [0, 2, 4, 6, 2 * z2, 2 + 2 * z2, 4 + 2 * z2, 6 + 2 * z2, ...
                              4 * z2, 2 + 4 * z2, 4 + 4 * z2, 6 + 4 * z2, ...
                              6 * z2, 2 + 6 * z2, 4 + 6 * z2, 6 + 6 * z2];
                tt = tt_base + tt_offsets;
    
                % Compute displacement using dot product (faster than explicit summation)
                % Extract x-displacement and y-displacement separately
                new_u_x = sum(sum(reshape(beta(tt), 4, 4) .* (zz' * zz2))); 
                new_u_y = sum(sum(reshape(beta(tt + 1), 4, 4) .* (zz' * zz2))); 
                
                % Store them in a 2-element vector
                new_u = [new_u_x; new_u_y]; 
    
                % Compute new pixel position
                ii = i - new_u(1);
                jj = j - new_u(2);
    
                % Compute gradients using bilinear interpolation
                a1 = BilinearApp(dyT, [ii; jj]); % Gradient in y
                a2 = BilinearApp(dxT, [ii; jj]); % Gradient in x
    
                % Store nonzero Jacobian entries efficiently
                for o = 0:3
                    for w = 0:3
                        coeff = zz(o + 1) * zz2(w + 1);
                        newa1 = a1 * coeff;
                        newa2 = a2 * coeff;
    
                        if abs(newa1) > 1e-3
                            Index = Index + 1;
                            I(Index) = (i - 1) * n + j;
                            J(Index) = t1 + t2 + 2 * z2 * o + 2 * w;
                            V(Index) = -newa1;
                        end
                        if abs(newa2) > 1e-3
                            Index = Index + 1;
                            I(Index) = (i - 1) * n + j;
                            J(Index) = t1 + t2 + 2 * z2 * o + 2 * w + 1;
                            V(Index) = -newa2;
                        end
                    end
                end
    
                % Compute residual using bilinear interpolation
                F((i - 1) * n + j) = BilinearApp(T, [i - new_u(1), j - new_u(2)]) - R(i, j);
            end
        end
    
        % Finalize sparse matrix construction
        I = I(1:Index);
        J = J(1:Index);
        V = V(1:Index);
        A = sparse(I, J, V, n * m, lG);
    end
    function [out] = DSSD(R,T,beta,k,z)
        out = 0;
        [m,n] = size(R);
        %outw = (norm(F(R,T,beta,k,z),2))^2;
        %outw = outw/2;
        for i = 1:m
        for j = 1:n
        new_u = BSplineTransformation([i,j],beta,k,z);
        out = out + (BilinearApp(T,[i-new_u(1),j-new_u(2)])-R(i,j))^2;
        end
        end
        out = out/(m*n);
        
    end
    function [out] = SimilarityMeasure(R, T, beta, k, z)
        out = 0;
        k1 = k(1);
        k2 = k(2);
        z2 = z(2);
        [m, n] = size(R);
    
        % Define B-spline basis function coefficients
        zz = @(u) [((1 - u)^3) / 6, (3*u^3 - 6*u^2 + 4) / 6, ...
                   (-3*u^3 + 3*u^2 + 3*u + 1) / 6, (u^3) / 6];
    
        % Loop over the image grid
        for i = 1:m
            y1 = floor((i - 1) / k1) * z2 * 2;
            u = (i - 1) / k1 - floor((i - 1) / k1);
            zz_i = zz(u);  % Basis function values for i
    
            for j = 1:n
                y2 = floor((j - 1) / k2) * 2 + 1;
                u2 = (j - 1) / k2 - floor((j - 1) / k2);
                zz_j = zz(u2);  % Basis function values for j
    
                % Compute control point indices
                t_base = y1 + y2;
                t = t_base + [0, 2, 4, 6, ...
                              2*z2, 2+2*z2, 4+2*z2, 6+2*z2, ...
                              4*z2, 2+4*z2, 4+4*z2, 6+4*z2, ...
                              6*z2, 2+6*z2, 4+6*z2, 6+6*z2];
    
                % Ensure beta values are properly reshaped to match 4×4 weighting
                beta_x = reshape(beta(t), [4, 4]);      % Extract x-displacement control points
                beta_y = reshape(beta(t + 1), [4, 4]);  % Extract y-displacement control points
    
                % Compute displacement with element-wise multiplication
                new_u_x = sum(sum(beta_x .* (zz_i' .* zz_j)));  % Use element-wise multiplication
                new_u_y = sum(sum(beta_y .* (zz_i' .* zz_j)));  
    
                new_u = [new_u_x; new_u_y];
    
                % Compute similarity measure
                out = out + (BilinearApp(T, [i - new_u(1), j - new_u(2)]) - R(i, j))^2;
            end
        end
    
        out = out / 2;
    end


end
