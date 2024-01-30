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

void Particle::setNeighbours(Particle down, Particle up, Particle right, Particle left, Particle diagonal_r_down,
                             Particle diagonal_r_up, Particle diagonal_l_down, Particle diagonal_l_up) {

    neighbour_up = &up;
    neighbour_down = &down;
    neighbour_right = &right;
    neighbour_left = &left;
    neighbour_diagonal_left_down = &diagonal_l_down;
    neighbour_diagonal_left_up = &diagonal_l_up;
    neighbour_diagonal_right_down = &diagonal_r_down;
    neighbour_diagonal_right_up = &diagonal_r_up;


}

Particle *Particle::getNeighbourRight() const {
    return neighbour_right;
}

Particle *Particle::getNeighbourLeft() const {
    return neighbour_left;
}

Particle *Particle::getNeighbourUp() const {
    return neighbour_up;
}

Particle *Particle::getNeighbourDown() const {
    return neighbour_down;
}

Particle *Particle::getNeighbourDiagonalRightDown() const {
    return neighbour_diagonal_right_down;
}

Particle *Particle::getNeighbourDiagonalLeftDown() const {
    return neighbour_diagonal_left_down;
}

Particle *Particle::getNeighbourDiagonalRightUp() const {
    return neighbour_diagonal_right_up;
}

Particle *Particle::getNeighbourDiagonalLeftUp() const {
    return neighbour_diagonal_left_up;
}

void Particle::setNeighbourRight(Particle *neighbourRight) {
    neighbour_right = neighbourRight;
}

void Particle::setNeighbourLeft(Particle *neighbourLeft) {
    neighbour_left = neighbourLeft;
}

void Particle::setNeighbourUp(Particle *neighbourUp) {
    neighbour_up = neighbourUp;
}

void Particle::setNeighbourDown(Particle *neighbourDown) {
    neighbour_down = neighbourDown;
}

void Particle::setNeighbourDiagonalRightDown(Particle *neighbourDiagonalRightDown) {
    neighbour_diagonal_right_down = neighbourDiagonalRightDown;
}

void Particle::setNeighbourDiagonalLeftDown(Particle *neighbourDiagonalLeftDown) {
    neighbour_diagonal_left_down = neighbourDiagonalLeftDown;
}

void Particle::setNeighbourDiagonalRightUp(Particle *neighbourDiagonalRightUp) {
    neighbour_diagonal_right_up = neighbourDiagonalRightUp;
}

void Particle::setNeighbourDiagonalLeftUp(Particle *neighbourDiagonalLeftUp) {
    neighbour_diagonal_left_up = neighbourDiagonalLeftUp;
}

const std::vector<Particle> &Particle::getDiagonalNeighbours() const {
    return DiagonalNeighbours;
}


const std::vector<Particle> &Particle::getLateralNeighbours() const {
    return LateralNeighbours;
}
/**
 * calculates if this particle and the other particle are neighbours
 * @param p2 = the other particle
 * @result 0,1, or 2.
 * 0 means they're not neighbours,
 * 1 means they're lateral neighbours,
 * and 2 means they're diagonal neighbours*/
int Particle::isNeighbours(Particle &p2) {
    double x1 = this->getX()[0];
    double y1 = this->getX()[1];


    double x2 = p2.getX()[0];
    double y2 = p2.getX()[1];

    //if they're lateral neighbours
    if(x1 == x2 && (y1 == y2 + 1 || y2 == y1 +1 )){
        return 1;
    }
    if(y1 == y2 && (x1 == x2 + 1 || x2 == x1 +1 )){
        return 1;
    }
    if(x1 == x2 +1 || x2 == x1 +1){
        if(y1 == y2 + 1 || y2 == y1 +1){
            return 2;
        }

    }

    return 0;
}

std::vector<Particle*> setDiagonalNeighbours(Particle &p){
    std::vector<Particle*> diagNeighbours;

    diagNeighbours.push_back(p.getNeighbourDiagonalRightUp());
    diagNeighbours.push_back(p.getNeighbourDiagonalRightDown());
    diagNeighbours.push_back(p.getNeighbourDiagonalLeftUp());
    diagNeighbours.push_back(p.getNeighbourDiagonalLeftDown());

    return diagNeighbours;

}
std::vector<Particle*> setLateralNeighbours(Particle &p){
    std::vector<Particle*> latNeighbours;

    latNeighbours.push_back(p.getNeighbourRight());
    latNeighbours.push_back(p.getNeighbourLeft());
    latNeighbours.push_back(p.getNeighbourUp());
    latNeighbours.push_back(p.getNeighbourDown());

    return latNeighbours;
}












