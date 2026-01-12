%functions
function nested = main_functions_2d()
    % Return handles to subfunctions
    nested.BSplineTransformation = @BSplineTransformation;
    nested.ImageRegistration = @ImageRegistration;
    function[out]=BSplineTransformation(x,beta,k,z)
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
        zz = [((1 - u)^3) / 6, ...
               (3 * u^3 - 6 * u^2 + 4) / 6, ...
               (-3 * u^3 + 3 * u^2 + 3 * u + 1) / 6, ...
               (u^3) / 6];
        zz2 = [((1 - u2)^3) / 6, ...
               (3 * u2^3 - 6 * u2^2 + 4) / 6, ...
               (-3 * u2^3 + 3 * u2^2 + 3 * u2 + 1) / 6, ...
               (u2^3) / 6];

        % Precompute index offsets
        tBase = y2 + y1;
        offsets = [0, 2, 4, 6];
        zOffsets = 2 * z2 * (0:3);

        % Compute output
        out = zeros(2, 1);  % Initialize output as 2D vector
        for i = 1:4
            for j = 1:4
                idx = tBase  + offsets(j) + zOffsets(i);
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
    % R = template, T = shifted image, k = [k1, k2, ...] step size
    % Here 2D model, so k = [k1, k2]
    
    % Test if the dimensions of k are correct
    test = length(k);
    if test ~= 2
    error('k does not have enough arguments');
    end
    clear test
    
    % n = [n1, n2] Number of pixels in x- and y-direction
    [n(1), n(2)] = size(R);
    
    % Test if R and T have the same dimensions
    [n2(1), n2(2)] = size(T);
    if n(1) ~= n2(1) || n(2) ~= n2(2)
    error('R and T do not have the same dimensions');
    end
    clear n2
    
    % z = [z1, z2] Number of control points in x- and y-direction
    % Since the image is surrounded by control points lying outside, +3
    z = [0; 0];
    for i = 1:2
    z(i) = ceil(n(i) / k(i)) + 3;
    end
    
    % lG = number of control points
    lG = z(1) * z(2);
    
    % beta = vector of all control point positions
    % For each control point, x- and y-coordinates are stored sequentially
    % Initial vector = 0
    beta = zeros(2 * lG, 1);
    
    % Later to be made into a while loop
    lambda = 2; % Regularization parameter
    p = 0.3;
    l = 0.7;
    test = 99999;
    new = 0;
    Iteration = 1;
    Dssd = zeros(1, 100);
    Dssd(Iteration) = DSSD(R, T, beta, k, z);
    [dxT, dyT] = imgradientxy(T, 'central');
    next = SimilarityMeasure(R, T, beta, k, z);
    Steps = [1, 0.5, 0];
    v = 2;
    
    for Iteration = 1:20 % Iterative optimization loop
    Iteration
    [J, f] = Gradient(R, T, dxT, dyT, beta, k, z);
    
    JJ = J' * J;
    JJ = Add(JJ, lambda);
    s = (JJ) \ ((-J)' * f);
    alpha = Steps(1);
    step = 1;
    test = next;
    next = SimilarityMeasure(R, T, beta + alpha * s, k, z);
    
    while (next > 0.995 * test && alpha > 0)
    step = step + 1;
    alpha = Steps(step);
    
    % Due to rounding errors
    if alpha < 0.01
        next = test;
    else
        next = SimilarityMeasure(R, T, beta + alpha * s, k, z);
    end
    end
    
    next
    
    if test == next
    r = 0;
    else
    r = (test - next) / (-s' * alpha * J' * f);
    end
    
    beta = beta + alpha * s;
    
    if r < p % Reduce lambda
    if test - new < 0
        lambda = lambda * 10;
    else
        lambda = lambda / 10;
    end
    else
    if r > l % Increase lambda
        lambda = lambda / 10;
    end
    end
    
    Iteration = Iteration + 1;
    Dssd(Iteration) = DSSD(R, T, beta, k, z);
    end
    
    Image = zeros(n(1), n(2));
    
    for i = 1:n(1)
        for j = 1:n(2)
            new_u = BSplineTransformation([i, j], beta, k, z);
            Image(i, j) = BilinearApp(T, [i - new_u(1), j - new_u(2)]);
        end
    end
    
    Dssd = Dssd(1:Iteration - 1);
    end
    
    function [A] = Add(A, lambda)
    [m, n] = size(A);
    
    for i = 1:m
    A(i, i) = A(i, i) + lambda;
    end
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
    
    
    function [out] = SimilarityMeasure(R,T,beta,k,z)
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


    function [out] =  BilinearApp(T,x)
    
    
    [m,n] = size(T);
    y11 = floor(x(1));
    y12 = floor(x(2));
    y21 = ceil(x(1));
    y22= ceil(x(2));
    
    if x(1)>= 1 && x(1)<= m && x(2)>= 1 && x(2)<= n
    a11 = T(y11,y12); 
    a12 = T(y11,y22);
    a21 = T(y21,y12);
    a22 = T(y21,y22);
    else
    if y11 < 1 || y11 > m
    a11 = 0;
    a12 = 0;
    else
    if y12 < 1 || y12 > n
       a11 = 0;
    else
       a11 = T(y11,y12);
    end 
    if y22< 1 || y22> n
       a12 = 0;
    else
       a12 = T(y11,y22);
    end
    end
    
    if y21 < 1 || y21 > m
    a21 = 0;
    a22 = 0;
    else
    if y12 < 1 || y12 > n
        a21 = 0;
    else
        a21 = T(y21,y12);
    end
    if y22< 1 || y22> n
        a22 = 0;
    else
        a22 = T(y21,y22);
    end
    end
    end
    
    
    %Since the step size is always one, there is no need to consider the difference quotient, only the distance.
    if y21 ~= y11
    z1 = (y21-x(1))*a11 + (x(1)-y11)*a21;
    z2 = (y21-x(1))*a12 + (x(1)-y11)*a22;
    else
    z1 = a21;
    z2 = a22;
    end
    
    if y22~= y12
    out = (y22-x(2))*z1 + (x(2)-y12)*z2;
    else
    out = z2;
    end
    
    end

    function [A, F] = Gradient(R, T, dxT, dyT, beta, k, z)
    % GradientAndResidual computes the Jacobian matrix (A) and residual vector (F)
    % for image registration using B-spline transformations.
    % 
    % Inputs:
    %   R - Reference image (2D matrix)
    %   T - Transformed image (2D matrix)
    %   dxT, dyT - Partial derivatives of T with respect to x and y
    %   beta - Control points for B-spline transformation
    %   k - Step sizes in x and y directions
    %   z - Number of control points in x and y directions
    % Outputs:
    %   A - Sparse Jacobian matrix of size (m*n, lG)
    %   F - Residual vector of size (m*n, 1)
    
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
end
