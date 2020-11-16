#include "matrix.hpp"
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <vector>
#include <numeric>

using namespace std;

ostream &operator<<(ostream &out, const matrix &mt) {
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

bool matrix::operator==(const matrix& m1) const {
	if (rows() != m1.rows() || cols() != m1.cols()) {
		return false;
	}
	return std::equal(matrix_.begin(), matrix_.end(),
			m1.matrix_.begin(), m1.matrix_.end());
}

bool matrix::operator!=(const matrix& m1) const {
	return !(*this == m1);
}

matrix& matrix::operator+=(const matrix& m2) {
	if (rows() != m2.rows()) {
		throw std::invalid_argument("Rows not matching");
	}
	if (cols() != m2.cols()) {
		throw std::invalid_argument("Columns not matching");
	}
	std::transform(matrix_.begin(), matrix_.end(), 
			m2.matrix_.begin(), matrix_.begin(), std::plus<value_type>());
	return *this;
}

matrix operator+(matrix m1, const matrix& m2) {
	m1 += m2;
	return m1;
}

matrix& matrix::operator-=(const matrix& m2) {
        if (rows() != m2.rows()) {
                throw std::invalid_argument("Rows not matching");
        }
        if (cols() != m2.cols()) {
                throw std::invalid_argument("Columns not matching");
        }
        std::transform(matrix_.begin(), matrix_.end(),
                        m2.matrix_.begin(), matrix_.begin(), std::minus<value_type>());
        return *this;
}


matrix operator-(matrix m1, const matrix& m2) {
	m1 -= m2;
	return m1;
}

matrix& matrix::operator*=(int k) {
	std::transform(matrix_.begin(), matrix_.end(), matrix_.begin(),
			[k](auto x){ return k * x; });
	return *this;
}

matrix operator*(matrix m1, const int k) {
	m1 *= k;
	return m1;
}

matrix operator*(const int k, matrix m1) {
	m1 *= k;
	return m1;
}

matrix operator*(matrix m1, const matrix& m2) {
	matrix result = matrix(m1.size().first, m2.size().second);
	std::pair <uint,uint> coordinates_target;
	std::pair <uint,uint> coordinates_m1;
	std::pair <uint,uint> coordinates_m2;
	for (uint x = 0; x < m1.size().first; x++) {
		for (uint y = 0; y < m2.size().second; y++) {
			coordinates_target = std::make_pair(x,y); 
			result[coordinates_target] = 0;
			for (uint xx = 0; xx < m1.size().first; xx++){
				coordinates_m1 = std::make_pair((x + xx)%m1.size().first,y); 
				coordinates_m2 = std::make_pair(x,(y + xx)%m1.size().first); 
				result[coordinates_target] += m1[coordinates_m1] * m2[coordinates_m2];
			}
		}
	}
	return result;
}

matrix::value_type& matrix::operator[](std::pair<size_t, size_t> index){
	if (index.first >= rows() || index.second >= cols()) {
		throw std::out_of_range("Index out of range");
	}
	return matrix_[index.first * cols() + index.second];
}

const matrix::value_type& matrix::operator[](std::pair<size_t, size_t> index) const {
	return matrix_[index.first * cols() + index.second];
}
