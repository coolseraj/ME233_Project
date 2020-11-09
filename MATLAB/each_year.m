function [loc] = each_year(und, num_p)
%this function finds the number of each undergrad group in a specific
%loaction
for i = 1:length(und)
    ratio(i) = und(i)/sum(und);
end
for i = 1:length(num_p)
for j = 1:length(und)
    loc(i,j) = round(num_p(i)*ratio(j));
end
end
end

