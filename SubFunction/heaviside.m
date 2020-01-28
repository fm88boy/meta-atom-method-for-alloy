% Heavisde½×Ô¾º¯Êı
function Y = heaviside(X)
Y = zeros(size(X));
Y(X > 0) = 1;
end
