function Mc = matlab2c(Mm);

% function Mc = matlab2c(Mm);

Mc = Mm*[1 0 0 1; 0 1 0 1; 0 0 1 1; 0 0 0 1]; % Convert from 1-based to 0-based indexing

return
