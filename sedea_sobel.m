function y = sedea_sobel(x)

%sobel coefs
h = [1 2 1;...
    0 0 0:...
    -1 -2 -1];
y = abs(conv2(x, h));
end