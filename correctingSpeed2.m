function newSpeed = correctingSpeed2(Param, v, trial)
    longi = [cos(Param.prefdir),sin(Param.prefdir)];
    ortho = [-sin(Param.prefdir),cos(Param.prefdir)];
    x = Param.decodedPos;
    
    k = 0.005;
    Magic = [99.97;97.86;98.81;97.78;94.03;92.67;98.14;99.76];
    magic = Magic(Param.idx,1);
    
    error = abs(x*ortho');
    attractor = magic*longi+trial.startHandPos;
    correction = k*error*(attractor-x);%/norm(attractor-x);
    
    newSpeed = v+correction;
   
end