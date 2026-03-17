clc; clear;

//============================================================
// LATTICE POSITION GENERATOR
//============================================================
function [a1,a2,a3] = lattice_pos(index, L)
    n       = ceil(index^(1/3));
    spacing = L / n;
    count   = 1;
    for i = 0:n-1
        for j = 0:n-1
            for k = 0:n-1
                if count == index then
                    a1 = (i+0.5)*spacing;
                    a2 = (j+0.5)*spacing;
                    a3 = (k+0.5)*spacing;
                    return
                end
                count = count + 1;
            end
        end
    end
endfunction

//============================================================
// LJ ENERGY FOR SINGLE PARTICLE j
//============================================================
function en = particle_energy(pos, j, npart, L, rc)
    en = 0;  rc2 = rc*rc;
    for k = 1:npart
        if k <> j then
            xr = pos(j,1)-pos(k,1);  xr = xr - L*round(xr/L);
            yr = pos(j,2)-pos(k,2);  yr = yr - L*round(yr/L);
            zr = pos(j,3)-pos(k,3);  zr = zr - L*round(zr/L);
            r2 = xr^2 + yr^2 + zr^2;
            if r2 < rc2 & r2 > 0.64 then
                r2i = 1/r2;  r6i = r2i^3;
                en  = en + 4*r6i*(r6i - 1);
            end
        end
    end
endfunction

//============================================================
// TOTAL LJ ENERGY (all pairs)
//============================================================
function en = total_energy(pos, npart, L, rc)
    en = 0;  rc2 = rc*rc;
    for i = 1:npart-1
        for j = i+1:npart
            xr = pos(i,1)-pos(j,1);  xr = xr - L*round(xr/L);
            yr = pos(i,2)-pos(j,2);  yr = yr - L*round(yr/L);
            zr = pos(i,3)-pos(j,3);  zr = zr - L*round(zr/L);
            r2 = xr^2 + yr^2 + zr^2;
            if r2 < rc2 & r2 > 0.64 then
                r2i = 1/r2;  r6i = r2i^3;
                en  = en + 4*r6i*(r6i - 1);
            end
        end
    end
endfunction

//============================================================
// SIMULATION PARAMETERS
//============================================================
npart = 100;
L     = 7.5;
temp  = 1.0;
rc    = 2.5;
delta = 0.2;
steps = 3000;
equil = 300;

pos = zeros(npart, 3);
for i = 1:npart
    [x,y,z]  = lattice_pos(i, L);
    pos(i,1) = x;  pos(i,2) = y;  pos(i,3) = z;
end

en = total_energy(pos, npart, L, rc);
printf('Initial E/N = %.4f\n', en/npart);

fID  = mopen('ENE_MC100.dat',  'wt');
fTrj = mopen('Traj_MC100.xyz', 'wt');
accept = 0;
PE     = zeros(1, steps);

//============================================================
// MONTE CARLO LOOP
//============================================================
for step = 1:steps

    for sweep = 1:npart
        j = int(rand()*npart) + 1;
        if j > npart then j = npart; end

        en_old = particle_energy(pos, j, npart, L, rc);
        xold = pos(j,1);  yold = pos(j,2);  zold = pos(j,3);

        pos(j,1) = modulo(pos(j,1) + (rand()-0.5)*delta, L);
        pos(j,2) = modulo(pos(j,2) + (rand()-0.5)*delta, L);
        pos(j,3) = modulo(pos(j,3) + (rand()-0.5)*delta, L);

        en_new = particle_energy(pos, j, npart, L, rc);
        dE     = en_new - en_old;

        if dE < 0 | rand() < exp(-dE/temp) then
            en     = en + dE;
            accept = accept + 1;
        else
            pos(j,1) = xold;
            pos(j,2) = yold;
            pos(j,3) = zold;
        end
    end

    PE(step) = en / npart;
    mfprintf(fID, '%d\t%.6f\n', step, en/npart);

    mfprintf(fTrj, '%d\n\n', npart);
    for i = 1:npart
        mfprintf(fTrj, 'Ar\t%.4f\t%.4f\t%.4f\n', ...
                 pos(i,1), pos(i,2), pos(i,3));
    end

    if modulo(step,300)==0 then
        printf('Step %4d | E/N = %8.5f | Acc = %.3f\n', ...
               step, en/npart, accept/(step*npart));
    end

end
mclose(fID);  mclose(fTrj);

//============================================================
// PLOT
//============================================================
scf(1); clf();
plot(1:steps, PE, 'b-');
xlabel('Monte Carlo Steps');
ylabel('Potential Energy per Particle (epsilon)');
title('NVT MC Simulation -- LJ System  N=100  T*=1.0');
xgrid();
avg_E = mean(PE(equil:steps));
printf('Equilibrium E/N  = %.6f epsilon\n', avg_E);
printf('Acceptance ratio = %.3f\n', accept/(steps*npart));
