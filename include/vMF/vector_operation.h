//
// Created by tin on 1/08/22.
//

#ifndef VMF_VECTOR_OPERATION_H
#define VMF_VECTOR_OPERATION_H

#include <array>

namespace vector_operation {

    template<size_t Dim, typename RealType = double>
    RealType dot(const std::array<RealType, Dim> &a, const std::array<RealType, Dim> &b) {
        RealType result = 0;
        for (size_t i = 0; i < Dim; ++i)
            result += a[i] * b[i];
        return result;
    }

    template<size_t Dim, typename RealType = double>
    std::array<RealType, Dim> sub(const std::array<RealType, Dim> &a, const std::array<RealType, Dim> &b) {
        std::array<RealType, Dim> result;
        for (size_t i = 0; i < Dim; ++i)
            result[i] = a[i] - b[i];
        return result;
    }

    template<size_t Dim, typename RealType = double>
    RealType norm(const std::array<RealType, Dim> &a) {
        RealType result = 0;
        for (size_t i = 0; i < Dim; ++i) {
            result += a[i] * a[i];
        }
        return std::sqrt(result);
    }
}

#endif //VMF_VECTOR_OPERATION_H
