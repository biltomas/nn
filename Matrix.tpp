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
Matrix<T>& matmul(Matrix<T>& m1, Matrix<T>& m2, Matrix<T>& target) {
    if (m1.cols() != m2.rows()) {
        throw std::invalid_argument("Invalid matrix size");
    }
	std::pair <uint,uint> coordinates_m1;
	std::pair <uint,uint> coordinates_m2;
    /*
    const unsigned BLOCK_SIZE = 8;
	for (uint target_y = 0; target_y < m1.rows(); target_y++) {
	    for (uint target_x = 0; target_x < m2.cols(); target_x += BLOCK_SIZE) {
			for (uint shared_axis = 0; shared_axis < m1.cols(); shared_axis++){
				coordinates_m1 = std::make_pair(target_y, shared_axis); 
                for (unsigned block = target_x; block < target_x + BLOCK_SIZE && block < m2.cols(); block++) {
				    coordinates_m2 = std::make_pair(shared_axis, block); 
				    target[{target_y, block}] += m1[coordinates_m1] * m2[coordinates_m2];
                }
			}
		}
	}*/
    #pragma omp parallel for shared(m1, m2, target) num_threads(8) 
	for (uint target_x = 0; target_x < m2.cols(); target_x++) {
	    for (uint target_y = 0; target_y < m1.rows(); target_y++) {
			for (uint shared_axis = 0; shared_axis < m1.cols(); shared_axis++) {
                target[{target_y, target_x}] += m1[{target_y, shared_axis}] * m2[{shared_axis, target_x}];
			}
		}
	}
	return target;
}


template <typename T>
T& Matrix<T>::operator[](std::pair<size_t, size_t> index){
	if (index.first >= rows() || index.second >= cols()) {
		throw std::out_of_range("Index out of range");
	}
    if (transposed) {
        return matrix_[index.second * rows() + index.first];
    }
	return matrix_[index.first * cols() + index.second];
}

template <typename T>
const T& Matrix<T>::operator[](std::pair<size_t, size_t> index) const {
    if (transposed) {
        return matrix_[index.second * rows() + index.first];
    }
	return matrix_[index.first * cols() + index.second];
}

template <typename T>
Matrix<T>& Matrix<T>::transpose() {
    /*
	Matrix<T> newMatrix(cols_, rows_);
	for (unsigned col = 0; col < cols_; col++) {
		for (unsigned row = 0; row < rows_; row++) {
			newMatrix[{col, row}] = operator[]({row, col});
		}
	}
    */
    transposed = !transposed;
    std::swap(rows_, cols_);
    return *this;
};

template <typename T>
void Matrix<T>::setRandom() {
    // std::generate(v.begin(), v.end(), (float) std::rand/RAND_MAX);
    for (size_t i = 0; i < matrix_.size(); i++) {
        matrix_[i] = (static_cast<T>(std::rand()) / (static_cast<T>(RAND_MAX)))/100;
        // std::cout << matrix_[i] << " ";
    }
}

template <typename T>
void Matrix<T>::setZero() {
    // std::generate(v.begin(), v.end(), (float) std::rand/RAND_MAX);
    for (size_t i = 0; i < matrix_.size(); i++) {
        matrix_[i] = 0;
        // cout << v[i] << " ";
    }
}

template <typename T>
void Matrix<T>::setNumber(float number) {
    // std::generate(v.begin(), v.end(), (float) std::rand/RAND_MAX);
    for (size_t i = 0; i < matrix_.size(); i++) {
        matrix_[i] = number + ((static_cast<T>(std::rand()) / (static_cast<T>(RAND_MAX)))/100);
        // std::cout << matrix_[i] << " ";
    }
}
