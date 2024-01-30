/*
 * Particle.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include <array>
#include <string>
#include <vector>

class Particle {

private:
    /**
     * Position of the particle
     */
    std::array<double, 3> x{};

    /**
     * Velocity of the particle
     */
    std::array<double, 3> v{};

    /**
     * Force effective on this particle
     */
    std::array<double, 3> f{};

    /**
     * Force which was effective on this particle
     */
    std::array<double, 3> old_f{};

    /**
     * Mass of this particle
     */
    double m{};

    /**
     * Type of the particle. Use it for whatever you want (e.g. to separate
     * molecules belonging to different bodies, matters, and so on)
     */
    int type;

    /**
     * sigma value for Lennard-Jones force calculation
     */
    double sig{};

    /**
     * epsilon value for Lennard-Jones force calculation
     */
    double eps{};

    /**
     * neighbours of a particle
     */
     Particle *neighbour_right;
public:
    Particle *getNeighbourRight() const;

    Particle *getNeighbourLeft() const;

    Particle *getNeighbourUp() const;

    Particle *getNeighbourDown() const;

    Particle *getNeighbourDiagonalRightDown() const;

    Particle *getNeighbourDiagonalLeftDown() const;

    Particle *getNeighbourDiagonalRightUp() const;

    Particle *getNeighbourDiagonalLeftUp() const;

private:
    Particle *neighbour_left;
public:
    void setNeighbourRight(Particle *neighbourRight);

    void setNeighbourLeft(Particle *neighbourLeft);

    void setNeighbourUp(Particle *neighbourUp);

    void setNeighbourDown(Particle *neighbourDown);

    void setNeighbourDiagonalRightDown(Particle *neighbourDiagonalRightDown);

    void setNeighbourDiagonalLeftDown(Particle *neighbourDiagonalLeftDown);

    void setNeighbourDiagonalRightUp(Particle *neighbourDiagonalRightUp);

    void setNeighbourDiagonalLeftUp(Particle *neighbourDiagonalLeftUp);

private:
    Particle *neighbour_up ;
     Particle *neighbour_down;
    Particle *neighbour_diagonal_right_down;
    Particle *neighbour_diagonal_left_down;
    Particle *neighbour_diagonal_right_up;
    Particle *neighbour_diagonal_left_up;
    std::vector<Particle> LateralNeighbours;
    std::vector<Particle> DiagonalNeighbours;
public:



    const std::vector<Particle> &getDiagonalNeighbours() const;
    const std::vector<Particle> &getLateralNeighbours() const;

    std::vector<Particle*> setDiagonalNeighbours(Particle &p);
    std::vector<Particle*> setLateralNeighbours(Particle &p);



public:
    explicit Particle(int type = 0);

    Particle(
            const Particle &other);

    Particle(
            // for visualization, we need always 3 coordinates
            // -> in case of 2d, we use only the first and the second
            std::array<double, 3> x_arg, std::array<double, 3> v_arg,
            double m_arg, double sig, double eps, int type = 0);


    virtual ~Particle();

    [[nodiscard]] const std::array<double, 3> &getX() const;

    [[nodiscard]] const std::array<double, 3> &getV() const;

    [[nodiscard]] const std::array<double, 3> &getF() const;

    [[nodiscard]] const std::array<double, 3> &getOldF() const;

    [[nodiscard]] double getM() const;

    [[nodiscard]] int getType() const;

    [[nodiscard]] double getSig() const;

    [[nodiscard]] double getEps() const;

     void setNeighbours(Particle right, Particle left, Particle up, Particle down,
                                     Particle diagonal_r_down, Particle diagonal_r_up,
                                     Particle diagonal_l_down,
                                     Particle diagonal_l_up);


    bool operator==(Particle &other);

    [[nodiscard]]int isNeighbours( Particle &p2, double h);

    [[nodiscard]] std::string toString() const;

    void setX(std::array<double, 3> newX);
    void setV(std::array<double, 3> newV);
    void setF(std::array<double, 3> newF);
    void setOldF(std::array<double, 3> oldF);



};

