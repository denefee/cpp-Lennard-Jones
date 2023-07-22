#include <iostream>
#include <cmath>

struct wector {
    public:
        long double x;
        long double y;
        long double z;
        void multiply(wector W, long double n) {
            W.x = W.x * n;
            W.y = W.y * n;
            W.z = W.z * n;            
        }
        void multiplyW(wector W1, wector W2) {
            W1.x = W1.x * W2.x;
            W1.y = W1.y * W2.y;
            W1.z = W1.z * W2.z;         
        }
}

class Particle {


};

/*class Particle:
    """Particle class"""
    def __init__(self, c = np.zeros(3), v = np.zeros(3), a = np.zeros(3), lc = np.zeros(3), way = np.zeros(3)):
        self.c = c # coordinate
        self.v = v # velocity
        self.a = a # acceleration
        self.lc = lc # last coordinate
        self.way = way # the movement of a particle from the beginning of time
        
    def to_border(c):
        # returns the particle to the borders of the box
        for i in np.arange(3):
            while ((c[i] >= Leng)or(c[i] < 0)):
                c[i] %= Leng
             
    def vec_to_virtual_copy(partc, part1c):
        # returns a vector directed to a virtual copy of particle "part1"
        vect_r = part1c - partc
        for i in np.arange(3):
            if (vect_r[i] > half):
                vect_r[i] -= Leng
            if (vect_r[i] < -half):
                vect_r[i] += Leng
        return vect_r
       
    def first_move(self):
        # moves the particle for the first time 
        self.lc = np.zeros(3) + self.c
        delta_r = dt*(self.v) + 0.5*(self.a)*dt**2
        self.way += delta_r
        self.c += delta_r
        Particle.to_border(self.c)
        self.v += dt*(self.a)
    
    def move(self):
        # moves the particle using the Verlet scheme
        delta_r = Particle.vec_to_virtual_copy(self.lc, self.c) + self.a*dt**2
        self.lc = np.zeros(3) + self.c
        self.way += delta_r
        self.c += delta_r
        Particle.to_border(self.c)
        self.v += self.a*dt

import math
import numpy as np
import time

from particleclass import Particle

# глобальные переменные
N = int(512)  # количество частиц
Vmax = float(2.0)  # максимальная скорость частицы
dt = float(0.001)  # тик
Leng = int(16)  # длина коробки
half = Leng/2  # половина длины коробки


# opens files with datas
# impt = open('imp.txt', 'w')
kint = open('kin.txt', 'w')
pott = open('pot.txt', 'w')
mect = open('mec.txt', 'w')
wayt = open('way.txt', 'w')
# coord = open('coord.txt', 'w')
maxwtx = open('maxwx.txt', 'w')
maxwty = open('maxwy.txt', 'w')
maxwtz = open('maxwz.txt', 'w')


def cell_gen(particles):
    # cell generation
    n = 0
    particle_is_even = True
    edge = math.ceil(np.cbrt(N))
    dl = Leng/edge
    dl_half = dl/2
    for i in np.arange(edge):
        for j in np.arange(edge):
            for k in np.arange(edge):
                c = dl_half + np.array([i, j, k])*dl
                if particle_is_even:
                    v = np.random.uniform(-Vmax, Vmax, (3))
                    if (n == N-1):
                        v = np.zeros(3)
                        particles.append(Particle(c, v))
                        return 0
                    particles.append(Particle(c, v))
                    particle_is_even = False
                    n += 1
                else:
                    particles.append(Particle(c, -v))
                    if (n == N-1):
                        return 0
                    particle_is_even = True
                    n += 1


def null_axel(particles):
    # nullifies all accelerations
    for i in np.arange(N):
        particles[i].a = np.zeros(3)


def axel(part, part1):
    # calculates the forces of interaction between these particles
    # and changes their accelerations
    vect_r = Particle.vec_to_virtual_copy(part.c, part1.c)
    abs_r = np.linalg.norm(vect_r)
    ac = 24*(2*np.power(abs_r, -14) - np.power(abs_r, -8))*vect_r
    part.a -= ac
    part1.a += ac


def calc_axel(particles):
    # calculates the accelerations of all particles and changes them
    null_axel(particles)
    for i in np.arange(N-1):
        for j in np.arange(i+1, N):
            axel(particles[i], particles[j])


def first_move(particles):
    # moves all particles for the first time
    calc_axel(particles)
    for i in np.arange(N):
        Particle.first_move(particles[i])


def move(particles):
    # moves all particles
    calc_axel(particles)
    for i in np.arange(N):
        Particle.move(particles[i])


def potentwo(part, part1):
    # calculates the potential energy of the interaction of two particles
    vect_r = Particle.vec_to_virtual_copy(part.c, part1.c)
    abs_r = np.linalg.norm(vect_r)
    u = 4*(np.power(abs_r, -12) - np.power(abs_r, -6))
    return u


def impulse(particles):
    # calculates the total momentum of the system
    summ = np.zeros(3)
    for i in np.arange(N):
        summ += particles[i].v
    impt.write(np.array2string(summ) + '\n')


def poten_eng(particles):
    # calculates the potential energy of the interaction of all particles
    pot = 0.0
    for i in np.arange(N-1):
        for j in np.arange(i+1, N):
            pot += potentwo(particles[i], particles[j])
    pott.write(str(pot) + '\n')
    return pot


def kinetic_eng(particles):
    # calculates the total kinetic energy of the system
    kin = 0.0
    for i in np.arange(N):
        kin += (np.linalg.norm(particles[i].v)**2)/2
    kint.write(str(kin) + '\n')
    return kin


def energy(particles):
    # calculates the total mechanical energy of the system
    pot = poten_eng(particles)
    kin = kinetic_eng(particles)
    summ = pot + kin
    mect.write(str(summ) + '\n')


def maxwell(particles):
    listx = np.zeros(N)
    listy = np.zeros(N)
    listz = np.zeros(N)
    for i in np.arange(N):
        listx[i] = particles[i].v[0]
        listy[i] = particles[i].v[1]
        listz[i] = particles[i].v[2]
    listx = np.sort(listx)
    listy = np.sort(listy)
    listz = np.sort(listz)
    for i in np.arange(N):
        maxwtx.write(str(listx[i]) + '\n')
        maxwty.write(str(listy[i]) + '\n')
        maxwtz.write(str(listz[i]) + '\n')


def average_way(particles):
    summ = 0.0
    for i in np.arange(N):
        summ += (np.linalg.norm(particles[i].way))**2
    summ = summ/N
    wayt.write(str(summ) + '\n')

    
def display_coordinates(particles):
    coord.write(str(N) + '\n')
    coord.write('Lattice="8.0 0.0 0.0 0.0 8.0 0.0 0.0 0.0 8.0" Properties=S:1:pos:R:3' + '\n')
    for i in np.arange(N):
        for j in np.arange(3):
            coord.write(str(particles[i].c[j]) + ' ')
        coord.write('\n')


def timego(particles, tick):
    "starts the simulation"
    # display_coordinates(particles)
    print(0, '%')
    first_move(particles)
    # display_coordinates(particles)
    # impulse(particles)
    energy(particles)
    average_way(particles)
    for i in np.arange(1, tick):
        move(particles)
        # display_coordinates(particles)
        # impulse(particles) commented out because momentum is maintained
        energy(particles) # commented out because energy is maintained
        average_way(particles)
        if i % (tick//1000) == 0:
            print(i*100/tick, '%')
    maxwell(particles)
    print(100, '%')


def main():
    t = int(25000)  # ticks
    start = time.time()  # точка отсчета времени
    particles = [] # particle array
    cell_gen(particles)  # генерация сеткой  
    timego(particles, t)
    end = time.time() - start  # собственно время работы программы
    print(end)  # вывод времени


if __name__ == "__main__":
    main()


# close files with data
# impt.close()
kint.close()
mect.close()
pott.close()
wayt.close()
# coord.close()
maxwtx.close()
maxwty.close()
maxwtz.close()*/