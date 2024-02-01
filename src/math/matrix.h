#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <math.h>
#include <optional>

namespace ptracer {
namespace math {

// clang-format off
template<typename T>
class Matrix4x4 {
public:
    Matrix4x4() : m{{{1, 0, 0, 0}, 
                     {0, 1, 0, 0}, 
                     {0, 0, 1, 0}, 
                     {0, 0, 0, 1}}} {}

    Matrix4x4(T m00, T m01, T m02, T m03, 
              T m10, T m11, T m12, T m13, 
              T m20, T m21, T m22, T m23, 
              T m30, T m31, T m32, T m33)
        : m{{{m00, m01, m02, m03},
             {m10, m11, m12, m13},
             {m20, m21, m22, m23},
             {m30, m31, m32, m33}}} {}

    T& operator()(unsigned int i, unsigned int j) {
        return m[i][j];
    }
    const T& operator()(unsigned int i, unsigned int j) const {
        return m[i][j];
    }

    Matrix4x4<T> operator*(const Matrix4x4<T>& mat) const {
        Matrix4x4<T> ret;
        for (unsigned int i = 0; i < 4; ++i)
            for (unsigned int j = 0; j < 4; ++j)
                ret.m[i][j] = m[i][0] * mat.m[0][j] +
                              m[i][1] * mat.m[1][j] +
                              m[i][2] * mat.m[2][j] +
                              m[i][3] * mat.m[3][j];
        return ret;
    }

    Matrix4x4<T> operator*(T scalar) const {
        Matrix4x4<T> ret;
        for (unsigned int i = 0; i < 4; ++i)
            for (unsigned int j = 0; j < 4; ++j)
                ret[i][j] = m[i][j] * scalar;
        return ret;
    }

    std::optional<Matrix4x4<T>> inverse() const;

private:
    std::array<std::array<T, 4>, 4> m;
};

template<typename T>
inline Matrix4x4<T> operator*(T scalar, const Matrix4x4<T>& mat) {
    return mat * scalar;
}

template<typename T>
inline Matrix4x4<T> operator*(const Matrix4x4<T>& lhs, const Matrix4x4<T>& rhs) {
    return lhs * rhs;
}

template<typename T>
inline std::optional<Matrix4x4<T>> Inverse(const Matrix4x4<T>& mat) {
    return mat.inverse();
}

template<typename T>
inline Matrix4x4<T> Transpose(const Matrix4x4<T>& mat) {
    return {mat(0, 0), mat(1, 0), mat(2, 0), mat(3, 0),
            mat(0, 1), mat(1, 1), mat(2, 1), mat(3, 1),
            mat(0, 2), mat(1, 2), mat(2, 2), mat(3, 2),
            mat(0, 3), mat(1, 3), mat(2, 3), mat(3, 3)};
}

template<typename T>
inline T Det3x3(const Matrix4x4<T>& mat) {
    return mat(0, 0) * (mat(1, 1) * mat(2, 2) - mat(1, 2) * mat(2, 1)) +
           mat(0, 1) * (mat(1, 2) * mat(2, 0) - mat(1, 0) * mat(2, 2)) +
           mat(0, 2) * (mat(1, 0) * mat(2, 1) - mat(1, 1) * mat(2, 0));
}
// clang-format on

template<typename T>
std::optional<Matrix4x4<T>> Matrix4x4<T>::inverse() const {
    Matrix4x4<T> ret = *this;

    std::array<int, 4> indxc, indxr;
    std::array ipiv{0, 0, 0, 0};

    for (int i = 0; i < 4; i++) {
        int irow = -1, icol = -1;
        T big = 0;
        for (int j = 0; j < 4; j++) {
            if (ipiv[j] != 1) {
                for (int k = 0; k < 4; k++) {
                    if (ipiv[k] == 0) {
                        if (std::abs(ret.m[j][k]) >= big) {
                            big  = std::abs(ret.m[j][k]);
                            irow = j;
                            icol = k;
                        }
                    } else if (ipiv[k] > 1) {
                        return {};
                    }
                }
            }
        }
        ++ipiv[icol];
        if (irow != icol) {
            for (int k = 0; k < 4; ++k)
                std::swap(ret.m[irow][k], ret.m[icol][k]);
        }
        indxr[i] = irow;
        indxc[i] = icol;
        if (ret.m[icol][icol] == 0)
            return {};
        T pivinv          = 1.0 / ret.m[icol][icol];
        ret.m[icol][icol] = 1.f;
        for (int j = 0; j < 4; j++)
            ret.m[icol][j] *= pivinv;
        for (int j = 0; j < 4; j++) {
            if (j != icol) {
                T save         = ret.m[j][icol];
                ret.m[j][icol] = 0;
                for (int k = 0; k < 4; k++)
                    ret.m[j][k] -= ret.m[icol][k] * save;
            }
        }
    }
    for (int j = 4 - 1; j >= 0; j--) {
        if (indxr[j] != indxc[j]) {
            for (int k = 0; k < 4; k++)
                std::swap(ret.m[k][indxr[j]], ret.m[k][indxc[j]]);
        }
    }

    return ret;
}

} // namespace math

using Mat4 = math::Matrix4x4<Float>;

} // namespace ptracer

#endif // __MATRIX_H__