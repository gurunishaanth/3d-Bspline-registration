%functions
function nested = main_functions()
    % Return handles to subfunctions
    nested.BSplineTransformation = @BSplineTransformation;
    nested.ImageRegistration = @ImageRegistration;
    function[out]=BSplineTransformation(x,beta,k,z)
    %x,k,g are 2D; beta is 2*lG x 1
    %k-Grid spacing
    %x-pixel coordinates
    %z-no of ctrl pts
    k1 = k(1);
    k2 = k(2);
    z2 = z(2);
    %Since there are negative control points, everything is shifted by +3.
    %beta stores control points row-wise, i.e., the j-th row.
    %Starts at index (j-1)*z(2)*2+1.
    y1 = floor((x(1)-1)/k(1))*z(2)*2;
    %Since beta stores x- and y-coordinates consecutively, corresponds to x2.
    %+1 because MATLAB starts counting from 1.
    y2 = floor((x(2)-1)/k(2))*2 +1;
    %Computes the simplified cubic spline.
    %i = i-th interval.
    %k = Grid spacing.
    %x = Value.
    %Order is reversed.
    u = (x(1)-1)/k1 - floor((x(1)-1)/k1);
    u2 = (x(2)-1)/k2 - floor((x(2)-1)/k2);
    
    zz1 = ((1-u)^3)/6;
    zz2 = (3*u^3 -6*u^2 +4)/6;
    zz3 = (-3*u^3 + 3*u^2 +3*u+1)/6;
    zz4 = (u^3)/6;
    
    zz21 = ((1-u2)^3)/6;
    zz22 = (3*u2^3 -6*u2^2 +4)/6;
    zz23 = (-3*u2^3 + 3*u2^2 +3*u2+1)/6;
    zz24 = (u2^3)/6;
    
    t   = y2 +y1;
    t1  = y2 +y1+2;
    t2  = y2 +y1+4;
    t3  = y2 +y1+6;
    t4  = y2 +y1+2*z2;
    t5  = y2 +y1+2+2*z2;
    t6  = y2 +y1+4+2*z2;
    t7  = y2 +y1+6+2*z2;
    t8  = y2 +y1+4*z2;
    t9  = y2 +y1+2+4*z2;
    t10 = y2 +y1+4+4*z2;
    t11 = y2 +y1+6+4*z2;
    t12 = y2 +y1+6*z2;
    t13 = y2 +y1+2+6*z2;
    t14 = y2 +y1+4+6*z2;
    t15 = y2 +y1+6+6*z2;
    
    
    out = beta(t:t+1)*zz1*zz21;
    out = out + beta(t1:t1+1)*zz1*zz22;
    out = out + beta(t2:t2+1)*zz1*zz23;
    out = out + beta(t3:t3+1)*zz1*zz24;
    out = out + beta(t4:t4+1)*zz2*zz21;
    out = out + beta(t5:t5+1)*zz2*zz22;
    out = out + beta(t6:t6+1)*zz2*zz23;
    out = out + beta(t7:t7+1)*zz2*zz24;
    out = out + beta(t8:t8+1)*zz3*zz21;
    out = out + beta(t9:t9+1)*zz3*zz22;
    out = out + beta(t10:t10+1)*zz3*zz23;
    out = out + beta(t11:t11+1)*zz3*zz24;
    out = out + beta(t12:t12+1)*zz4*zz21;
    out = out + beta(t13:t13+1)*zz4*zz22;
    out = out + beta(t14:t14+1)*zz4*zz23;
    out = out + beta(t15:t15+1)*zz4*zz24;
    
    % for i = 0:3              
    %     for j = 0:3
    %         t = y(2)+i*z(2)*2  +y(1)+j*2;
    %         out = out + beta(t:t+1)*zz(i+1)*zz2(j+1);
    %     end
    % end
    
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
    
    for Iteration = 1:5 % Iterative optimization loop
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
    [m,n] = size(R);
    
    
    %outw = (norm(F(R,T,beta,k,z),2))^2;
    %outw = outw/2;
    
    for i = 1:m
    y1 = floor((i-1)/k2)*z2*2;
    u = (i-1)/k1 - floor((i-1)/k1);
    zz1 = ((1-u)^3)/6;
    zz2 = (3*u^3 -6*u^2 +4)/6;
    zz3 = (-3*u^3 + 3*u^2 +3*u+1)/6;
    zz4 = (u^3)/6;
    for j = 1:n
    %       new_u = BSplineTransformation([i,j],beta,k,z);
    y2 = floor((j-1)/k2)*2 +1;
    u2 = (j-1)/k2 - floor((j-1)/k2);
    
    zz21 = ((1-u2)^3)/6;
    zz22 = (3*u2^3 -6*u2^2 +4)/6;
    zz23 = (-3*u2^3 + 3*u2^2 +3*u2+1)/6;
    zz24 = (u2^3)/6;
    
    t   = y2 +y1;
    t1  = y2 +y1+2;
    t2  = y2 +y1+4;
    t3  = y2 +y1+6;
    t4  = y2 +y1+2*z2;
    t5  = y2 +y1+2+2*z2;
    t6  = y2 +y1+4+2*z2;
    t7  = y2 +y1+6+2*z2;
    t8  = y2 +y1+4*z2;
    t9  = y2 +y1+2+4*z2;
    t10 = y2 +y1+4+4*z2;
    t11 = y2 +y1+6+4*z2;
    t12 = y2 +y1+6*z2;
    t13 = y2 +y1+2+6*z2;
    t14 = y2 +y1+4+6*z2;
    t15 = y2 +y1+6+6*z2;
    
    
    new_u = beta(t:t+1)*zz1*zz21;
    new_u = new_u + beta(t1:t1+1)*zz1*zz22;
    new_u = new_u + beta(t2:t2+1)*zz1*zz23;
    new_u = new_u + beta(t3:t3+1)*zz1*zz24;
    new_u = new_u + beta(t4:t4+1)*zz2*zz21;
    new_u = new_u + beta(t5:t5+1)*zz2*zz22;
    new_u = new_u + beta(t6:t6+1)*zz2*zz23;
    new_u = new_u + beta(t7:t7+1)*zz2*zz24;
    new_u = new_u + beta(t8:t8+1)*zz3*zz21;
    new_u = new_u + beta(t9:t9+1)*zz3*zz22;
    new_u = new_u + beta(t10:t10+1)*zz3*zz23;
    new_u = new_u + beta(t11:t11+1)*zz3*zz24;
    new_u = new_u + beta(t12:t12+1)*zz4*zz21;
    new_u = new_u + beta(t13:t13+1)*zz4*zz22;
    new_u = new_u + beta(t14:t14+1)*zz4*zz23;
    new_u = new_u + beta(t15:t15+1)*zz4*zz24;
    out = out + (BilinearApp(T,[i-new_u(1),j-new_u(2)])-R(i,j))^2;
    end
    end
    out = out/2;
    end
    
    function [out] =  BilinearApp(T,x)
    
    % a11---- a12
    %   |           |
    %   |           |
    % a21---- a22
    
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
    
    k1 = k(1);
    k2 = k(2);
    z2 = z(2);
    [m, n] = size(T);      % Image dimensions
    [lG, ~] = size(beta);  % Number of control points
    I = zeros(n * m * 32, 1);
    J = zeros(n * m * 32, 1);
    V = zeros(n * m * 32, 1);
    F = zeros(n * m, 1);
    Index = 0;
    
    for i = 1:m
    t1 = floor((i - 1) / k1) * z2 * 2;  % Index offset for control points
    u = (i - 1) / k1 - floor((i - 1) / k1); % Relative position in the grid
    
    % B-spline basis functions for x
    zz(1) = ((1 - u)^3) / 6;
    zz(2) = (3 * u^3 - 6 * u^2 + 4) / 6;
    zz(3) = (-3 * u^3 + 3 * u^2 + 3 * u + 1) / 6;
    zz(4) = (u^3) / 6;
    
    for j = 1:n
    t2 = floor((j - 1) / k2) * 2 + 1;  % Index offset for control points
    u2 = (j - 1) / k2 - floor((j - 1) / k2); % Relative position in the grid
    
    % B-spline basis functions for y
    zz2(1) = ((1 - u2)^3) / 6;
    zz2(2) = (3 * u2^3 - 6 * u2^2 + 4) / 6;
    zz2(3) = (-3 * u2^3 + 3 * u2^2 + 3 * u2 + 1) / 6;
    zz2(4) = (u2^3) / 6;
    
    % B-spline transformation indices
    tt   = t2 + t1;
    tt1  = t2 + t1 + 2;
    tt2  = t2 + t1 + 4;
    tt3  = t2 + t1 + 6;
    tt4  = t2 + t1 + 2 * z2;
    tt5  = t2 + t1 + 2 + 2 * z2;
    tt6  = t2 + t1 + 4 + 2 * z2;
    tt7  = t2 + t1 + 6 + 2 * z2;
    tt8  = t2 + t1 + 4 * z2;
    tt9  = t2 + t1 + 2 + 4 * z2;
    tt10 = t2 + t1 + 4 + 4 * z2;
    tt11 = t2 + t1 + 6 + 4 * z2;
    tt12 = t2 + t1 + 6 * z2;
    tt13 = t2 + t1 + 2 + 6 * z2;
    tt14 = t2 + t1 + 4 + 6 * z2;
    tt15 = t2 + t1 + 6 + 6 * z2;
    
    % Compute displacement due to control points
    new_u = beta(tt:tt + 1) * zz(1) * zz2(1) ...
          + beta(tt1:tt1 + 1) * zz(1) * zz2(2) ...
          + beta(tt2:tt2 + 1) * zz(1) * zz2(3) ...
          + beta(tt3:tt3 + 1) * zz(1) * zz2(4) ...
          + beta(tt4:tt4 + 1) * zz(2) * zz2(1) ...
          + beta(tt5:tt5 + 1) * zz(2) * zz2(2) ...
          + beta(tt6:tt6 + 1) * zz(2) * zz2(3) ...
          + beta(tt7:tt7 + 1) * zz(2) * zz2(4) ...
          + beta(tt8:tt8 + 1) * zz(3) * zz2(1) ...
          + beta(tt9:tt9 + 1) * zz(3) * zz2(2) ...
          + beta(tt10:tt10 + 1) * zz(3) * zz2(3) ...
          + beta(tt11:tt11 + 1) * zz(3) * zz2(4) ...
          + beta(tt12:tt12 + 1) * zz(4) * zz2(1) ...
          + beta(tt13:tt13 + 1) * zz(4) * zz2(2) ...
          + beta(tt14:tt14 + 1) * zz(4) * zz2(3) ...
          + beta(tt15:tt15 + 1) * zz(4) * zz2(4);
    
    % Residual and gradients
    ii = i - new_u(1);
    jj = j - new_u(2);
    a1 = BilinearApp(dyT, [ii; jj]); % Gradient in y
    a2 = BilinearApp(dxT, [ii; jj]); % Gradient in x
    
    for o = 0:3
        for w = 0:3
            % Contribution from each basis function
            newa = a1 * zz(o + 1) * zz2(w + 1);
            if abs(newa) > 1e-3
                Index = Index + 1;
                I(Index) = (i - 1) * n + j;
                J(Index) = t1 + t2 + 2 * z2 * o + 2 * w;
                V(Index) = -newa;
            end
    
            newa = a2 * zz(o + 1) * zz2(w + 1);
            if abs(newa) > 1e-3
                Index = Index + 1;
                I(Index) = (i - 1) * n + j;
                J(Index) = t1 + t2 + 2 * z2 * o + 2 * w + 1;
                V(Index) = -newa;
            end
        end
    end
    
    % Compute residual
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
