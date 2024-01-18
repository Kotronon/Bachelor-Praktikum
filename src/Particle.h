/*
 * Particle.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include <array>
#include <string>

class Particle {

private:
    /**
     * Position of the particle
     */
    std::array<double, 3> x{};

    /**
     * Previous position of the particle
     */
    std::array<double , 3> oldX{};

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
     * true if the particle is fixed on a specific position (e.g. wall)
     */
    bool fixed;

public:
    explicit Particle(int type = 0);

    Particle(
            const Particle &other);

    Particle(
            // for visualization, we need always 3 coordinates
            // -> in case of 2d, we use only the first and the second
            std::array<double, 3> x_arg, std::array<double, 3> v_arg,
            double m_arg, double sig, double eps, int type = 0, bool fixed = false);

    virtual ~Particle();

    [[nodiscard]] const std::array<double, 3> &getX() const;

    [[nodiscard]] const std::array<double, 3> &getOldX() const;

    [[nodiscard]] const std::array<double, 3> &getV() const;

    [[nodiscard]] const std::array<double, 3> &getF() const;

    [[nodiscard]] const std::array<double, 3> &getOldF() const;

    [[nodiscard]] double getM() const;

    [[nodiscard]] int getType() const;

    [[nodiscard]] double getSig() const;

    [[nodiscard]] double getEps() const;

    [[nodiscard]] bool getFixed() const;

    bool operator==(Particle &other);

    [[nodiscard]] std::string toString() const;

    void setX(std::array<double, 3> newX);

    void setOldX(std::array<double, 3> newOldX);

    void setV(std::array<double, 3> newV);

    void setF(std::array<double, 3> newF);

    void setOldF(std::array<double, 3> oldF);
};

