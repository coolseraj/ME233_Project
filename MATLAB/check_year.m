function [year_L, year_M] = check_year(L, M, und)
%checks the year of L and M particles
if L <= und(1)
    year_L = 1; 
elseif und(1) < L <= und(2)
    year_L = 2;
elseif und(2) < L <= und(3)
    year_L = 3;
else
    year_L = 4;
end

if M <= und(1)
    year_M = 1; %we can assign another probability for infection rate
elseif und(1) < M <= und(2)
    year_M = 2;
elseif und(2) < M <= und(3)
    year_M = 3;
else
    year_M = 4;
end
end

