%% Process data file

function [xq,yq,vq] = process(q12_list, q23_list, k_data)
    x = q12_list;
    y = q23_list;
    %v = kox_list.*prefactor.*prefactor;
    v = k_data;
    %v = dos_max;

    [xq,yq] = meshgrid(1:.01:5, 1:.01:5);

    % Duplicate values
    [theta,uid]  = unique([x(:), y(:)], 'rows');
    % Non-unique elements
    nuid = 1:1:length(x);
    nuid(uid) = [];
    % Check non-unique elements
    % q12_list(nuid)
    % q23_list(nuid)

    % Avoid duplicate values
    ux = theta(:,1);
    uy = theta(:,2);
    uv = v(uid);
    vq = griddata(ux,uy,uv,xq,yq, 'cubic');
end
