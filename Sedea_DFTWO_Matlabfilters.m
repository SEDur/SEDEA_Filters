function output = Sedea_DFTWO_Matlabfilters(input, bcoefs, acoefs)

% myinput = zeros(size(input) + size(acoefs));
% myinput(size(acoefs):end) = input;

dels = zeros(1,length(acoefs)-1);

output = zeros(size(input));

    for i = 1 : 1 : size(input)
        xsamps = sum(dels .* -acoefs(2:end));
        xsamps = xsamps + input(i);
        ysamps = sum(dels .* bcoefs(2:end));
        output(i) = (xsamps * bcoefs(1)) + ysamps;
        dels(2:end) = dels(1:length(dels)-1);
        dels(1) = xsamps;
    end
end