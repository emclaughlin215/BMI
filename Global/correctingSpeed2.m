function newSpeed = correctingSpeed2(Param, x, v)
    longi = [cos(Param.prefdir),sin(Param.prefdir)];
    ortho = [-sin(Param.prefdir), cos(Param.prefdir)];
    
    k = 0.003;
    
    error = abs(x*ortho');
    attractor = norm(x)*longi;
    correction = k*error*(attractor-x);%/norm(attractor-x);
    
    newSpeed = v+correction;
   
end