function [center, radii, chi2] = computeSphereFit(xyz)
    %% least square fitting 3D point cloud to sphere shell
    % 
    % INPUT:
    %   xyz:
    %       N by 3 matrix, as x,y,z coordinates for each localization
    %
    %
    % OUTPUT:
    %   center:
    %       1 by 3 vector, as x,y,z coordinate of the fitted sphere center
    %
    %   radii:
    %       scalar, as fitted sphere radius
    %
    %   chi2:
    %       sum of residual error
    %
    % code adopted from: 
    %   Yury (2024). Ellipsoid fit 
    %   https://www.mathworks.com/matlabcentral/fileexchange/24693-ellipsoid-fit
    %
    %
    % <Ziqiang.Huang@embl.de>
    % 05.03.2024

    if size(xyz, 1) < 4
        return;
    end
    x = xyz(:,1); y = xyz(:,2); z = xyz(:,3);
    D = [ 2*x, 2*y, 2*z, 1 + 0*x ];
    % solve the normal system of equations
    d2 = x .* x + y .* y + z .* z; % the RHS of the llsq problem (y's)
    u = ( D' * D ) \ ( D' * d2 );  % solution to the normal equations
    v = [ -1 -1 -1 0 0 0 u( 1 : 4 )' ];
    v = v';
    % form the algebraic form of the sphere
    A = [ v(1) v(4) v(5) v(7); ...
          v(4) v(2) v(6) v(8); ...
          v(5) v(6) v(3) v(9); ...
          v(7) v(8) v(9) v(10) ];
    % find the center of the sphere
    center = -A( 1:3, 1:3 ) \ v( 7:9 );
    % form the corresponding translation matrix
    T = eye( 4 );
    T( 4, 1:3 ) = center';
    % translate to the center
    R = T * A * T';
    % solve the eigenproblem
    [ evecs, evals ] = eig( R( 1:3, 1:3 ) / -R( 4, 4 ) );
    radii = sqrt( 1 ./ diag( abs( evals ) ) );
    sgns = sign( diag( evals ) );
    radii = radii .* sgns;
    % calculate difference of the fitted points from the actual data normalized by the conic radii
    d = [ x - center(1), y - center(2), z - center(3) ]; % shift data to origin
    d = d * evecs; % rotate to cardinal axes of the conic;
    d = [ d(:,1) / radii(1), d(:,2) / radii(2), d(:,3) / radii(3) ]; % normalize to the conic radii
    chi2 = sum( abs( 1 - sum( d.^2 .* repmat( sgns', size( d, 1 ), 1 ), 2 ) ) );

end