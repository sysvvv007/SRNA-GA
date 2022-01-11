function x = transform2to10(Population)
Bitlength = size(Population,2);
x = Population(Bitlength);
for i = 1:Bitlength-1
    x = x+Population(Bitlength-i)*power(2,i);
end
end