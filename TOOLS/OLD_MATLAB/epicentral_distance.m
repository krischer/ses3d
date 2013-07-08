function e=epicentral_distance(lat_1,lon_1,lat_2,lon_2)

	lat_1=lat_1*pi/180;
	lat_2=lat_2*pi/180;
	lon_1=lon_1*pi/180;
	lon_2=lon_2*pi/180;

	v_1(1)=cos(lon_1)*cos(lat_1);
	v_1(2)=sin(lon_1)*cos(lat_1);
	v_1(3)=sin(lat_1);
	
    v_2(1)=cos(lon_2)*cos(lat_2);
	v_2(2)=sin(lon_2)*cos(lat_2);
	v_2(3)=sin(lat_2);
    
    d=v_1(1)*v_2(1)+v_1(2)*v_2(2)+v_1(3)*v_2(3);
    
    e=acos(d)*180/pi;