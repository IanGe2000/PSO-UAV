function test
    particle = rand(1,10)*10 - 5 + linspace(0,30,10);
    particle(1) = 0;
    particle(end) = 30;
    particle = particle - (particle < 0) .* particle - (particle > 30) .* particle + (particle > 30)*30;
    plot(particle)
    hold on;
    adjustedparticle = particleAdjust (particle);
    plot(adjustedparticle)
    hold off;
endfunction

function adjustedparticle = particleAdjust (particle)
    [v, i] = min(diff(particle)>=0);
    if v == 0
        [w, j] = max((particle(i:end) - particle(i))>0);
        if w == 1
            particle(i+1) = particle(i) + (particle(i+j-1) - particle(i))/(j-1)*rand();
            adjustedparticle = particleAdjust(particle);
        else
            particle(i:end) = particle(i);
            adjustedparticle = particle;
        endif
    else
        adjustedparticle = particle;
    endif
endfunction