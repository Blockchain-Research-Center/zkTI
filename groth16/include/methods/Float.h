#ifndef __float_h
#define __float_h

#include <cassert>
#include <iostream>
#include <cmath>
#include <vector>
#include <map>

#define __float_constant_
#define __BASE 2
#define __MAXIMUM_BASE_POW 32 // MUST BE EVEN
#define __MINIMUM_BASE_POW 256

// #define __BASE 10
// #define __MAXIMUM_BASE_POW 10 // MUST BE EVEN
// #define __MINIMUM_BASE_POW 20

using namespace std;

typedef unsigned long long u64;
typedef short int i16;

class Float {

public:
    i16 scale; // E
    u64 value; // S
    const u64 maximum_value = get_maximum_value();
    const float minimum_value = get_minimum_value();

    u64 zero_point = 0; // Z

    float real_value; // a = BASE^{-E} * (S - Z)

    static Float one() {
        return Float(__MAXIMUM_BASE_POW, get_maximum_value(), 1.0);
    }

    static Float zero() {
        return Float(__MINIMUM_BASE_POW, 0, 0.0);
    }

    static u64 get_maximum_value() {
        return (u64)pow(__BASE, __MAXIMUM_BASE_POW) - 1;
    }

    static float get_minimum_value() {
        return 0;
    }

public:
    Float() {
        this->value = 0;
        this->scale = __MINIMUM_BASE_POW;
        this->real_value = 0;
    }

    Float(i16 scale, u64 value, float real_value, u64 zero_point=0) {
        this->scale = scale;
        this->value = value;
        this->real_value = real_value;
        this->zero_point = zero_point;
    }

    Float(float x) {
        this->real_value = x;
        this->scale = 0;
        this->value = 1;

        while(x < 1) {
            x *= __BASE;
            this->scale++;
        }

        if(scale >= __MAXIMUM_BASE_POW) {
            // 0
            this->value = 0;
            this->scale = __MINIMUM_BASE_POW;
            return;
        }

        // Start with the largest possible scale
        this->scale += __MAXIMUM_BASE_POW;
        this->value = (u64)(x * maximum_value);

        // If value is too large, decrease scale
        while (this->value > maximum_value) {
            this->scale--;
            this->value /= __BASE;
        }
    }

    Float(int x) : Float((float)x) {
        
    }

    Float one_minus() {
        u64 value = ((get_maximum_value()) * (u64)pow(__BASE, this->scale - __MAXIMUM_BASE_POW)) - this->value;
        i16 scale = this->scale;

        while (value < maximum_value) {
            scale++;
            value *= __BASE;
        }

        while (value > maximum_value) {
            scale--;
            value /= __BASE;
        }

        return Float(scale, value, 1 - this->real_value);
    }

    void operator=(const Float& other)
    {
        real_value = other.real_value;
        value = other.value;
        scale = other.scale;
    }

    Float operator+(const Float& other)
    {
        float real_value = this->real_value + other.real_value;
        // if(abs(this->scale - other.scale) >32) {
        //     std::cout << abs(this->scale - other.scale) << *this << other << std::endl;
        // }
        // if scaling distance > (__MAXIMUM_BASE_POW / 2), omit smaller one
        if(this->scale >= other.scale + (__MAXIMUM_BASE_POW + 2)) {
            return Float(other.scale, other.value, real_value);
        } else if(other.scale >= this->scale + (__MAXIMUM_BASE_POW + 2)){
            return Float(this->scale, this->value, real_value);
        }

        // std::cout << this->scale << " " << other.scale << " ";

        i16 scale = max(this->scale, other.scale);
        u64 this_value = this->value;
        u64 other_value = other.value;

        if (this->scale < scale) {
            // std::cout << scale - this->scale;
            this_value *= pow(__BASE, scale - this->scale);
        } else {
            // std::cout << scale - other.scale;
            other_value *= pow(__BASE, scale - other.scale);
        }

        u64 value = this_value + other_value;
        
        int i = 0;
        while (value > maximum_value) {
            scale--;
            value /= __BASE;
            i++;
        }
        // std::cout << " " << i << " " << scale << std::endl;

        // std::cout << "scale: " << abs(this->scale - other.scale) << " shift: " << i << std::endl;

        return Float(scale, value, real_value);
    }

    Float operator*(const Float& other)
    {
        float real_value = this->real_value * other.real_value;
        i16 scale = this->scale + other.scale;
        u64 value = this->value * other.value;

        int i = 0;
        while(value > maximum_value) {
            scale--;
            value /= __BASE;
            i++;
        }

        // std::cout << pow(2, 32) * (this->value * other.value - (u64)(value * pow(2, i)));
        // std::cout << " ";
        // std::cout << this->value * other.value;
        // std::cout << " compare: " << (pow(2, 32) * (this->value * other.value - (u64)(value * pow(2, i))) < (this->value * other.value ));
        // std::cout << std::endl;
        // std::cout << i << " " << *this << " * " << other << std::endl;

        return Float(scale, value, real_value);
    }

    Float operator/(const Float& other)
    {
        float real_value = this->real_value / other.real_value;
        i16 scale = this->scale + __MAXIMUM_BASE_POW - other.scale;
        // if(scale < 0) {
        //     std::cout << "111: " << *this << other << std::endl;
        // }
        u64 value = this->value * maximum_value / other.value;

        // int i = 0;
        while(value > maximum_value) {
            scale--;
            value /= __BASE;
            // i++;
        }

        // std::cout << i << std::endl;

        return Float(scale, value, real_value);
    }

    friend ostream &operator<<(ostream &out, const Float &f )
    { 
        out << "M: " << f.value << "*" << __BASE << "^\{-" << f.scale << "\}  " << "R: " << f.real_value;
        return out;            
    }

    
};

#endif