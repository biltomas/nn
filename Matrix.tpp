#include "Matrix.hpp"
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <vector>
#include <numeric>


template <typename T>
std::ostream &operator<<(std::ostream &out, const Matrix<T> &mt) {
	for (size_t row = 0; row < mt.rows(); row++) {
		const char *sep = "";
		for (size_t col = 0; col < mt.cols(); col++) {
			out << sep << mt[ {row, col} ];
			sep = ", ";
		}
		out << '\n';
	}
	return out;
}

template <typename T>
bool Matrix<T>::operator==(const Matrix<T>& m1) const {
	if (rows() != m1.rows() || cols() != m1.cols()) {
		return false;
	}
	return std::equal(matrix_.begin(), matrix_.end(),
			m1.Matrix_.begin(), m1.Matrix_.end());
}

template <typename T>
bool Matrix<T>::operator!=(const Matrix<T>& m1) const {
	return !(*this == m1);
}

template <typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& m2) {
	if (rows() != m2.rows()) {
		throw std::invalid_argument("Rows not matching");
	}
	if (cols() != m2.cols()) {
		throw std::invalid_argument("Columns not matching");
	}
	std::transform(matrix_.begin(), matrix_.end(), 
			m2.Matrix_.begin(), matrix_.begin(), std::plus<T>());
	return *this;
}

template <typename T>
Matrix<T>& operator+(Matrix<T> m1, const Matrix<T>& m2) {
	m1 += m2;
	return m1;
}

template <typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& m2) {
        if (rows() != m2.rows()) {
                throw std::invalid_argument("Rows not matching");
        }
        if (cols() != m2.cols()) {
                throw std::invalid_argument("Columns not matching");
        }
        std::transform(matrix_.begin(), matrix_.end(),
                        m2.matrix_.begin(), matrix_.begin(), std::minus<T>());
        return *this;
}


template <typename T>
Matrix<T> operator-(Matrix<T> m1, const Matrix<T>& m2) {
	m1 -= m2;
	return m1;
}


template <typename T>
Matrix<T>& Matrix<T>::operator*=(const T k) {
	std::transform(matrix_.begin(), matrix_.end(), matrix_.begin(),
			[k](auto x){ return k * x; });
	return *this;
}

template <typename T>
Matrix<T> operator*(Matrix<T> m1, const T k) {
	m1 *= k;
	return m1;
}

template <typename T>
Matrix<T> operator*(const T k, Matrix<T> m1) {
	m1 *= k;
	return m1;
}

template <typename T>
Matrix<T> operator*(const Matrix<T> m1, const Matrix<T>& m2) {
    if (m1.cols() != m2.rows()) {
        throw std::invalid_argument("Invalid matrix size");
    }
	Matrix<T> result(m1.rows(), m2.cols());
	std::pair <uint,uint> coordinates_target;
	std::pair <uint,uint> coordinates_m1;
	std::pair <uint,uint> coordinates_m2;
	// cout << "multiply" << endl;
	// cout << "m1:" << endl;
	// for (uint row = 0; row < m1.rows(); row++) {
	// 	for (uint col = 0; col < m1.cols(); col++) {
	// 		cout << m1[make_pair(row, col)] << " ";
	// 	}
	// 	cout << endl;
	// }
	// cout << "m2:" << endl;
	// for (uint row = 0; row < m2.rows(); row++) {
	// 	for (uint col = 0; col < m2.cols(); col++) {
	// 		cout << m2[make_pair(row, col)] << " ";
	// 	}
	// 	cout << endl;
	// }
	// cout << "m1 size: " << m1.size().first << ", " << m1.size().second << endl;
	// cout << "m2 size: " << m2.size().first << ", " << m2.size().second << endl;
	for (uint x = 0; x < m2.cols(); x++) {
		for (uint y = 0; y < m1.rows(); y++) {
			coordinates_target = std::make_pair(y,x); 
			result[coordinates_target] = 0;
			// cout << y << ", " << x << endl;
			for (uint xx = 0; xx < m1.cols(); xx++){
				coordinates_m1 = std::make_pair(y, xx); 
				coordinates_m2 = std::make_pair(xx, x); 
				// cout << "m1: " << coordinates_m1.first << ", " << coordinates_m1.second << endl;
				// cout << "m2: " << coordinates_m2.first << ", " << coordinates_m2.second << endl;
				result[coordinates_target] += m1[coordinates_m1] * m2[coordinates_m2];
			}
		}
	}
	// cout << "result:" << endl;
	// for (uint row = 0; row < result.rows(); row++) {
	// 	for (uint col = 0; col < result.cols(); col++) {
	// 		cout << result[make_pair(row, col)] << " ";
	// 	}
	// 	cout << endl;
	// }
	return result;
}

template <typename T>
T& Matrix<T>::operator[](std::pair<size_t, size_t> index){
	if (index.first >= rows() || index.second >= cols()) {
		throw std::out_of_range("Index out of range");
	}
	return matrix_[index.first * cols() + index.second];
}

template <typename T>
const T& Matrix<T>::operator[](std::pair<size_t, size_t> index) const {
	return matrix_[index.first * cols() + index.second];
}

template <typename T>
Matrix<T> Matrix<T>::transpose() {
	Matrix<T> newMatrix(cols_, rows_);
	for (unsigned col = 0; col < cols_; col++) {
		for (unsigned row = 0; row < rows_; row++) {
			newMatrix[{col, row}] = operator[]({row, col});
		}
	}
	return newMatrix;
	
};

template <typename T>
void Matrix<T>::setRandom() {
    std::pair<size_t, size_t> Matrix_size = size();

    std::vector<T> v(Matrix_size.first * Matrix_size.second);
    // std::generate(v.begin(), v.end(), (float) std::rand/RAND_MAX);
    for (size_t i = 0; i < v.size(); i++) {
        v[i] = static_cast<T>(std::rand()) / static_cast<T>(RAND_MAX);
        // cout << v[i] << " ";
    }
    *this = Matrix(v, Matrix_size.first);
}
