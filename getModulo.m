function [modI, modIminus, modIplus] = getModulo(i, n)

modI = mod(i,n);
if modI == 0
    modI = n; 
end

modIminus = modI - 1;
modIplus = modI + 1;

if modIminus == 0
    modIminus = n;
end

if modIplus > n
    modIplus = 1;
end

end