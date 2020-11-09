function [X, PV] = BC(XF,XR,X,PV)
%checks the BC and returns the value of X
while X<XF || X>XR
    if X<XF
        X = 2*XF - X; %specular reflection eq.11.7
        PV = -PV;
    end
    if X>XR
        X = 2*XR - X;
        PV = -PV;
    end
end
end

