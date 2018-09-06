function [L, Ldot] = domaingrowth(t,params)

num = params.Linfty*exp(params.a*(t-params.ts));
den = (params.Linfty/params.L0+exp(params.a*(t-params.ts))-1);

L = num / den + params.k0
%L = 1;
if nargout > 1
    Ldot = num*params.a/den - params.a*num*exp(params.a*(t-params.ts))/(den^2)
    %[L, Ldot]
    %Ldot = 0;
end

return