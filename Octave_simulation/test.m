function test
    segment = [0,0,1;0,2,3;0,3,1];
    theta = turningAngle (segment)
endfunction

function theta = turningAngle (segment)
## output format: theta is the angle formed by the segment in degree
    dif = diff(segment,1,2)
    (dif(1)*dif(4)+dif(2)*dif(5))/(sqrt((dif(1))^2+(dif(2))^2)*sqrt((dif(4))^2+(dif(5))^2))
    theta = acosd((dif(1)*dif(4)+dif(2)*dif(5))/(sqrt((dif(1))^2+(dif(2))^2)*sqrt((dif(4))^2+(dif(5))^2)));
endfunction