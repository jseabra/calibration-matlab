% -----------------------------------------
function R = pitch(ang, off)
% rotation around X

a = pi *ang / 180;
ca = cos(a);
sa = sin(a);

R = [1  , 0     , 0         ,off(1) ;
    0  ,ca , -sa   ,off(2) ;
    0  ,sa , ca    ,off(3) ;
    0  , 0     ,    0      , 1 ];

end
