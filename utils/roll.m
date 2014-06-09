% -----------------------------------------
% Rotation around Y
function R = roll(ang, off)

a = pi *ang / 180;
ca = cos(a);
sa = sin(a);

R = [ca , 0 , sa , off(1) ;
    0      , 1 , 0      , off(2) ;
    -sa, 0 , ca , off(3) ;
    0      , 0 , 0      , 1 ];

end