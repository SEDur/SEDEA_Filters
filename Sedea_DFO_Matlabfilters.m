function output = Sedea_DFO_Matlabfilters(input, acoefs, bcoefs)

% myinput = zeros(size(input) + size(acoefs));
% myinput(size(acoefs):end) = input;

xdel = zeros(1,length(acoefs));

ydel = zeros(1,length(bcoefs));

output = zeros(size(input));

for i = 1 : 1 : size(input)
    xdel(1) = input(i);
    xsamps = sum(xdel .* acoefs);
    ysamps = sum(ydel(2:end) .* -bcoefs(2:end));
    output(i) = xsamps + ysamps;
    xdel(2:end) = xdel(1:length(xdel)-1);
    ydel(2:end) = ydel(1:length(ydel)-1);
    ydel(1) = output(i);
end


end