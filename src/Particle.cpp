/*
 * Particle.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "Particle.h"

#include <iostream>
#include "utils/ArrayUtils.h"
#include "spdlog/spdlog.h"

/**
 * creates an empty particle
 * @param type_arg
 */
Particle::Particle(int type_arg) {
  type = type_arg;
  //spdlog::info("Particle generated!");
  f = {0., 0., 0.};
  old_f = {0., 0., 0.};
}

/**
 * makes an copy of another particle
 * @param other
 */
Particle::Particle(const Particle &other) {
  x = other.x;
  v = other.v;
  f = other.f;
  old_f = other.old_f;
  m = other.m;
  type = other.type;
  sig = other.sig;
  eps = other.eps;
  //spdlog::info("Particle generated by copy!");
}

/**
 * creates a new Particle
 * @param x_arg
 * @param v_arg
 * @param m_arg
 * @param sig
 * @param eps
 * @param type_arg
 */
Particle::Particle(std::array<double, 3> x_arg, std::array<double, 3> v_arg,
                   double m_arg, double sigma, double epsilon, int type_arg) {
    x = x_arg;
    v = v_arg;
    m = m_arg;
    type = type_arg;
    f = {0., 0., 0.};
    old_f = {0., 0., 0.};
    sig = sigma;
    eps = epsilon;
    //spdlog::info("Particle generated!");
}

Particle::~Particle() {
    //spdlog::info("Particle destructed!");
}

/**
 * returns x coordinates of particle
 * @return
 */
const std::array<double, 3> &Particle::getX() const { return x; }

/**
 * returns velocity of particle
 * @return
 */
const std::array<double, 3> &Particle::getV() const { return v; }

/**
 * returns force of particle
 * @return
 */
const std::array<double, 3> &Particle::getF() const { return f; }

/**
 * returns old force of particle
 * @return
 */
const std::array<double, 3> &Particle::getOldF() const { return old_f; }

/**
 * returns mass of particle
 * @return
 */
double Particle::getM() const { return m; }

/**
 * returns type of particle
 * @return
 */
int Particle::getType() const { return type; }

/**
 * returns sigma of particle
 * @return
 */
 double Particle::getSig() const {return sig;}

 /**returns epsilon of particle
  * @return
  */
  double Particle::getEps() const {return eps;}

/**
 * returns particle as string
 * @return
 */
std::string Particle::toString() const {
  std::stringstream stream;
  stream << "Particle: X:" << x << " v: " << v << " f: " << f
         << " old_f: " << old_f << " type: " << type;
  return stream.str();
}

/**
 * comparision with another particle
 * @param other
 * @return
 */
bool Particle::operator==(Particle &other) {
  return (x == other.x) && (v == other.v) && (f == other.f) &&
         (type == other.type) && (m == other.m) && (old_f == other.old_f);
}

/**
 * adds Particle to an given stream as string
 * @param stream
 * @param p
 * @return
 */
std::ostream &operator<<(std::ostream &stream, Particle &p) {
  stream << p.toString();
  return stream;
}

/**
 * sets x coordinates of Particle
 * @param x_arg
 */
void Particle::setX(std::array<double, 3> x_arg) {
   x = x_arg;
 }

 /**
  * sets velocity of Particle
  * @param v_arg
  */
void Particle::setV(std::array<double, 3> v_arg) {
  v = v_arg; 
}

/**
 * sets force of Particle
 * @param f_arg
 */
void Particle::setF(std::array<double, 3> f_arg){
  f = f_arg;
}

/**
 * sets old force of Particle
 * @param old_f_arg
 */
void Particle::setOldF(std::array<double, 3> old_f_arg) {
  old_f = old_f_arg;
}

