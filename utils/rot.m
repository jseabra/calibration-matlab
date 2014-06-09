function R = rot(ang , off)
% ROT Rotation around Z
% -----------------------------------------

a = pi *ang / 180;
ca = cos(a);
sa = sin(a);

R = [ca , -sa , 0 , off(1) ;
    sa , ca  , 0 , off(2) ;
    0      , 0       , 1 , off(3) ;
    0      , 0       , 0 , 1 ];
end