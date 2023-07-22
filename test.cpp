#include <iostream>
#include <cmath>

struct Wector {
    private:
        long double x, y, z;
    public:
        Wector(long double x0, long double y0, long double z0) : x(x0), y(y0), z(z0) 
        {}
        Wector() : x(0.0), y(0.0), z(0.0)  
        {}
        friend const Wector operator-(const Wector &vector);
        friend const Wector operator+(const Wector& left, const Wector& right);
        friend const Wector operator-(const Wector& left, const Wector& right);
        friend const Wector operator*(const Wector& vector, const long double scalar);
        friend const Wector operator*(const long double scalar, const Wector& vector);
        
        void display() {std::cout << x << ", " << y << ", " << z << std::endl;}

        long double quad(const Wector &vector) {
            long double quadratus = pow(vector.x, 2) + pow(vector.y, 2) + pow(vector.z, 2);
            return quadratus;
        }

        long double abs(const Wector &vector) {
            return sqrt(quad(vector));
        }
};

const Wector operator-(const Wector &vector) {
    long double X = -vector.x;
    long double Y = -vector.y;
    long double Z = -vector.z;
    return Wector(X, Y, Z);
}

const Wector operator+(const Wector &left, const Wector &right) {
    long double X = left.x + right.x;
    long double Y = left.y + right.y;
    long double Z = left.z + right.z;
    return Wector(X, Y, Z);
}

const Wector operator-(const Wector &left, const Wector &right) {
    long double X = left.x - right.x;
    long double Y = left.y - right.y;
    long double Z = left.z - right.z;
    return Wector(X, Y, Z);
}

const Wector operator*(const Wector &vector, const long double scalar) {
    long double X = vector.x * scalar;
    long double Y = vector.y * scalar;
    long double Z = vector.z * scalar;
    return Wector(X, Y, Z);
}

const Wector operator*(const long double scalar, const Wector &vector) {
    long double X = vector.x * scalar;
    long double Y = vector.y * scalar;
    long double Z = vector.z * scalar;
    return Wector(X, Y, Z);
}

int main() {
    Wector w1 = Wector(1, 2, 3);
    Wector w2 = Wector(1, 2, 3);
    Wector w;
    long double n = abs(w1);
    std::cout << n;
    // w.display();

    return 0;
}