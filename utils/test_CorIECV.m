function test_CorIECV

T1 = [

         0.813033420493566        -0.274470729561351         0.513972335651686                         0;
         0.102549274887616        -0.801830235774619        -0.588850169956354                         0;
         0.573114563922918         0.530786295614043        -0.623768004952609                         0;
                         0                         0                         0                         1];

T2 = [  0.825485886868745         0.348455628663437        -0.442126789303014                         0
        -0.119859777247723        -0.656145555952667        -0.744743180670404                         0
        -0.551549348996668         0.669366629177415        -0.499873082138257                         0
                         0                         0                         0                         1];
     

coord = Cor_IECV(T1);
rz = coord(4);
ry = coord(5);
rx = coord(6);

% IEC rotation:
T1_ = rot(rz,[0,0,0])*pitch(rx,[0,0,0])*roll(ry,[0,0,0])
T1-T1_

coord = Cor_IECV(T2);
rz = coord(4);
ry = coord(5);
rx = coord(6);

% IEC rotation:
T2_ = rot(rz,[0,0,0])*pitch(rx,[0,0,0])*roll(ry,[0,0,0])
T2-T2_


[x,y,z] = decompose_rotation(T2)
R = compose_rotation(x,y,z)

end

function [x,y,z] = decompose_rotation(R)
	x = atan2(R(3,2), R(3,3));
	y = atan2(-R(3,1), sqrt(R(3,2)*R(3,2) + R(3,3)*R(3,3)));
	z = atan2(R(2,1), R(1,1));
    
end

function R = compose_rotation(x, y, z)
	X = eye(3,3);
	Y = eye(3,3);
	Z = eye(3,3);

    X(2,2) = cos(x);
    X(2,3) = -sin(x);
    X(3,2) = sin(x);
    X(3,3) = cos(x);

    Y(1,1) = cos(y);
    Y(1,3) = sin(y);
    Y(3,1) = -sin(y);
    Y(3,3) = cos(y);

    Z(1,1) = cos(z);
    Z(1,2) = -sin(z);
    Z(2,1) = sin(z);
    Z(2,2) = cos(z);

	R = Z*Y*X;
    
end