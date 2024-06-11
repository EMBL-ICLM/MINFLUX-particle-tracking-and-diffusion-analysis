function sphericity = computeSphericity (xyz)
    %% compute sphericity of 3D points
    % 
    % INPUT:
    %   xyz:
    %       N by 3 matrix, as x,y,z coordinates for each localization
    % 
    %
    % OUTPUT:
    %   sphericity:
    %       a value in range [0, 1] representing how sphere-like the input
    %       3D point clouds are.
    %       sphericity = V^2 / A^3 * factor;
    %       factor = Vsphere ^2 / Asphere^3
    %              = (4/3*pi*R^3)^2 / (4*pi*R^2)^3
    %              = 1 / 36 / pi;
    %
    % <Ziqiang.Huang@embl.de
    % 05.03.2024
    
    % get convex hull of the point cloud
    [K, volume] = convhulln(xyz);
    
    % triangulation of surface area
    area= ...
        sum(sqrt(sum(( ...
        [xyz(K(:,1),2).*xyz(K(:,2),3) - xyz(K(:,1),3).*xyz(K(:,2),2) ...
        xyz(K(:,1),3).*xyz(K(:,2),1) - xyz(K(:,1),1).*xyz(K(:,2),3)  ...
        xyz(K(:,1),1).*xyz(K(:,2),2) - xyz(K(:,1),2).*xyz(K(:,2),1)] + ...
        [xyz(K(:,2),2).*xyz(K(:,3),3) - xyz(K(:,2),3).*xyz(K(:,3),2) ...
        xyz(K(:,2),3).*xyz(K(:,3),1) - xyz(K(:,2),1).*xyz(K(:,3),3)  ...
        xyz(K(:,2),1).*xyz(K(:,3),2) - xyz(K(:,2),2).*xyz(K(:,3),1)] + ...
        [xyz(K(:,3),2).*xyz(K(:,1),3) - xyz(K(:,3),3).*xyz(K(:,1),2) ...
        xyz(K(:,3),3).*xyz(K(:,1),1) - xyz(K(:,3),1).*xyz(K(:,1),3)  ...
        xyz(K(:,3),1).*xyz(K(:,1),2) - xyz(K(:,3),2).*xyz(K(:,1),1)]).^2,2))) ...
        /2;
    
    sphericity = (pi)^(1/3)*(6*volume)^(2/3)/area;

end