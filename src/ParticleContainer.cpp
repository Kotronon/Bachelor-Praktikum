//
// Created by maraconda on 02.11.23.
//

#include <iostream>
#include <vector>
#include <functional>
#include "spdlog/spdlog.h"
#include "ParticleContainer.h"
#include "utils/ArrayUtils.h"

/**
 * creates a new empty ParticleContainer
 */
ParticleContainer::ParticleContainer() {
   spdlog::info("Empty Particle Container generated!");
}

/**
 * creates a new ParticleContainer as a copy from other ParticleContainer
 * @param other
 */
ParticleContainer::ParticleContainer(const ParticleContainer &other) {
    containedParticles = other.containedParticles;
    spdlog::info("Particle Container generated by copy!");
}

ParticleContainer::~ParticleContainer() { spdlog::info("Particle Container destructed!"); }

/**
 * adds an existing Particle
 * @param particle
 */
void ParticleContainer::addParticle(const Particle &particle) {
    const Particle& new_particle(particle);
    containedParticles.emplace_back(new_particle);
    spdlog::info("Added copy of particle to container!");
}

/**
 * add a new empty particle to the container
 * @param type_arg
 */
void ParticleContainer::addParticle(int type_arg) {
    Particle new_particle(type_arg);
    containedParticles.emplace_back(new_particle);
    spdlog::info("Added newly created particle to container!");
}


/**
 * adds new Particle
 * @param x_arg
 * @param v_arg
 * @param m_arg
 * @param type_arg
 */
void ParticleContainer::addParticle(std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg, int type_arg, double sig, double eps) {
    Particle new_particle(x_arg, v_arg, m_arg, sig, eps, type_arg);
    containedParticles.emplace_back(new_particle);
    spdlog::info("Added newly created particle to container!");
}

/**
 * add anther ParticleContainer to the current container
 * @param container
 */
void ParticleContainer::addParticleContainer(ParticleContainer &container) {
    for (auto &p: container) {
        containedParticles.emplace_back(p);
    }
    spdlog::info("Added particles contained in other container to this container!");
}

/**
 * iterate through the Particles pairwise to calculate the force of each particle
 * @param forceCalculation
 */
void ParticleContainer::applyForcePairwise(const std::function<void(Particle*, Particle*)>& forceCalculation, double Ggrav){
    auto first = begin();
    auto last = end();
    for (; first != last; ++first) {
        for(auto next = std::next(first); next != last; ++next)
            forceCalculation(&(*first), &(*next));
        std::array<double, 3> gravity = {0, first->getM(), 0};
        first->setF(first->getF() + gravity);
    }
}

/**
 * returns the vector of all PArticles in the ParticleContainer with Pointer at first element
 * @return
 */
std::vector<Particle>::iterator ParticleContainer::begin() {
    return containedParticles.begin();
}

/**
 * returns the vector of all PArticles in the ParticleContainer with Pointer at last element
 * @return
 */
std::vector<Particle>::iterator ParticleContainer::end() {
    return containedParticles.end();
}

/**
 * returns number of Particles in ParticleContainer
 * @return
 */
int ParticleContainer::size() {
    return containedParticles.size();
}

/**
 * iterate through the Particles and set the neighbours
 * @param p Particle whose neighbours are supposed to be set
 * */

void ParticleContainer::SetAllNeighbours( Particle &p) {
    double x = p.getX()[0];
    double y = p.getX()[1];

    if(containedParticles.size() >1) {
        for (Particle particle: containedParticles) {
            if (!particle.operator==(p)) {
                double xp = particle.getX()[0];
                double yp = particle.getX()[1];

                if (xp == x + 1) {
                    if (yp == y + 1) {
                        p.setNeighbourDiagonalRightUp(&particle);
                        particle.setNeighbourDiagonalLeftDown(&p);
                    }
                    if (yp == y) {
                        p.setNeighbourRight(&particle);
                        particle.setNeighbourLeft(&p);
                    }

                    if(yp == y -1){
                        p.setNeighbourDiagonalRightDown(&particle);
                        particle.setNeighbourDiagonalLeftUp(&p);
                    }
                }

                if(xp == x){
                    if (yp == y + 1) {
                       p.setNeighbourUp(&particle);
                       particle.setNeighbourDown(&p);
                    }


                    if(yp == y -1){
                        p.setNeighbourDown(&particle);
                        particle.setNeighbourUp(&p);
                    }
                }

                if(xp == x-1){
                    if (yp == y + 1) {
                        p.setNeighbourDiagonalLeftUp(&particle);
                        particle.setNeighbourDiagonalRightDown(&p);
                    }
                    if (yp == y) {
                        p.setNeighbourLeft(&particle);
                        particle.setNeighbourRight(&p);
                    }

                    if(yp == y -1){
                        p.setNeighbourDiagonalLeftDown(&particle);
                        particle.setNeighbourDiagonalRightUp(&p);
                    }
                }

            }
        }
    }

   /* for (int i = 0; i < size(); ++i) {
        Particle p = containedParticles[i];

        int x = p.getX()[0];
        int y = p.getX()[1];
        int z = p.getX()[2];

        if(x ==0 && y == 0){
            break;
            return;

        }
        if(x ==0){

        }

    }*/

}
