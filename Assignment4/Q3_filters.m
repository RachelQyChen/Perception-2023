%Time filter
function [f1, f2] = Q3_filters(x, t, deltaT, tau)
    
    [x_size, y_size, t_length] = size(x);
    y = zeros(x_size, y_size, t_length, 7);
    f1 = zeros(x_size, y_size, t_length);
    f2 = zeros(x_size, y_size, t_length);

    for tt = 1:t_length - 1
        for type = 1:7
            if type == 1
                deltaY = (deltaT / tau) * (-y(:, :, tt, type) + x(:, :, tt));
                y(:, :, tt + 1, type) = y(:, :, tt, type) + deltaY;
            else
                deltaY = (deltaT / tau) * (-y(:, :, tt, type) + y(:, :, tt, type - 1));
                y(:, :, tt + 1, type) = y(:, :, tt, type) + deltaY;
            end
        end
        % filter 1
        f1(:, :, tt) = y(:, :, tt, 3) - y(:, :, tt, 5); 
         % filter 2
        f2(:, :, tt) = y(:, :, tt, 5) - y(:, :, tt, 7);
    end
end
