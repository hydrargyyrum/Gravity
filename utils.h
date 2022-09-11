#pragma once

#include <SFML/Graphics.hpp>
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <vector>

using namespace  std;



class vertex : public sf::Vertex
{
public:
    sf::Vector3f pos;
    float scale = 6.;

    vertex() { }

    vertex(float x, float y, float z, sf::Color color) {
        this->pos = sf::Vector3f(x, y, z);
        this->color = color;
    }
    
    void rotate(double & angle)
    {
        double i, j, k, x;
        i = pos.x;
        j = pos.y;
        k = pos.z;
        x = i*std::cos(angle)-j*std::sin(angle);
        auto vect = sf::Vector2f(x, k);
        position = scale*vect+sf::Vector2f(400,300);
    };

    vertex operator+(const vertex r) const
    {
        vertex l = *this;
        vertex res;
        res.color = l.color;
        res.pos = sf::Vector3f(
            l.pos.x + r.pos.x,
            l.pos.y + r.pos.y,
            l.pos.z + r.pos.z
        );
        return res;
    }


    void print() {
        std::cout << "x = " << pos.x << " y = " << pos.y << " z = " << pos.z << std::endl;
    }
};