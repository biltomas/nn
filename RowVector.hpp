#pragma once
#include "Matrix.hpp"
#include <iostream> 
#include <vector> 
#include <algorithm>
#include <string>
using namespace std; 

template <typename T = float>
class RowVector { 
    Matrix<T> vector_;
public: 

	RowVector(std::vector<T> vector) : vector_(vector, 1) {}

    RowVector(uint size) {
        std::vector<T> vector (size, T());
        vector_ = Matrix<T>(vector, 1);
    }

	void setValue(const uint pos, const float value) {
        vector_[{0, pos}] = value;
    }

    float coeffRef(const uint pos) const {
        if (pos >= length()) {
            throw std::out_of_range("Index out of range");
        }
        return vector_[{0, pos}];
    }

	T dot(const RowVector<T>& vector2) const {
        T result = 0;
        for (unsigned i = 0; i < length(); i++) {
            result += coeffRef(i) * vector2.coeffRef(i);
        }
        return result;
    }

    unsigned length() const { return vector_.size().second; }

    const Matrix<T>& data() const {
        return vector_;
    }
};

template <typename T>
RowVector<T> operator*(RowVector<T>& m1, const Matrix<T>& m2) {
    return RowVector<T>((m1.data() * m2).to_vector());
}

template <typename T>
RowVector<T> operator*(RowVector<T>& m1, const RowVector<T>& m2) {
    return RowVector<T>((m1.data() * m2).to_vector());
}


template <typename T>
RowVector<T> operator-(RowVector<T>& m1, const RowVector<T>& m2) {
	//return RowVector((m1.vector - m2.vector).to_vector());
    return RowVector<T>((m1.data() - m2.data()).to_vector());
}

template <typename T>
std::ostream &operator<<(std::ostream& out, const RowVector<T> &mt) {
    out << mt.data();
    return out;
}
