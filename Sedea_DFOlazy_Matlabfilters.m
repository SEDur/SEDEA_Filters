function output = Sedea_DFO_Matlabfilters(input, bcoefs, acoefs)

xdel = zeros(1, size(input)+ size(acoefs));
ydel = zeros(1, size(bcoefs));

xdel(1, size(acoefs)+1 :end) = input;

output = zeros(size(input));

for i = 1 : 1 : size(input)
    xsamps = sum(xdel(i-size(acoefs) : i) .* acoefs);
    ysamps = sum(ydel .* bcoefs);
    output(i) = (input(i) * acoefs(1)) + () + () + () + ()
end


end