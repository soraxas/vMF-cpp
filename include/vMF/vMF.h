//
// Created by tin on 1/08/22.
//

/**
 * def sample_vMF(mu, kappa, num_samples):
    """Generate num_samples N-dimensional samples from von Mises Fisher
    distribution around center mu \in R^N with concentration kappa.
    """
    dim = len(mu)
    result = np.zeros((num_samples, dim))
    for nn in range(num_samples):
        # sample offset from center (on sphere) with spread kappa
        w = _sample_weight(kappa, dim)

        # sample a point v on the unit sphere that's orthogonal to mu
        v = _sample_orthonormal_to(mu)

        # compute new point
        result[nn, :] = v * np.sqrt(1. - w**2) + w * mu

    return result


def _sample_weight(kappa, dim):
    """Rejection sampling scheme for sampling distance from center on
    surface of the sphere.
    """
    dim = dim - 1  # since S^{n-1}
    b = dim / (np.sqrt(4. * kappa**2 + dim**2) + 2 * kappa)
    x = (1. - b) / (1. + b)
    c = kappa * x + dim * np.log(1 - x**2)

    while True:
        z = np.random.beta(dim / 2., dim / 2.)
        w = (1. - (1. + b) * z) / (1. - (1. - b) * z)
        u = np.random.uniform(low=0, high=1)
        if kappa * w + dim * np.log(1. - x * w) - c >= np.log(u):
            return w


def _sample_orthonormal_to(mu):
    """Sample point on sphere orthogonal to mu."""
    v = np.random.randn(mu.shape[0])
    proj_mu_v = mu * np.dot(mu, v) / np.linalg.norm(mu)
    orthto = v - proj_mu_v
    return orthto / np.linalg.norm(orthto)
 */

#ifndef VMF_CPP_VMF_H
#define VMF_CPP_VMF_H

#include "beta_distribution.h"

#include <iostream>
#include <sstream>
#include <string>
#include <random>
#include <array>

#include "vector_operation.h"

namespace vMF {

    template<size_t Dim, typename RealType = double>
    class von_mises_fisher {
    public:

        using dim_vector = std::array<RealType, Dim>;

        typedef dim_vector result_type;

        class param_type {
        public:
            explicit param_type(dim_vector mu, RealType kappa)
                    : mu_param(mu), kappa_param(kappa) {}

            dim_vector mu() const { return mu_param; }

            RealType kappa() const { return kappa_param; }

            bool operator==(const param_type &other) const {
                for (size_t i = 0; i < Dim; ++i) {
                    if (mu_param[i] != other.mu_param[i] ||
                        kappa_param[i] == other.kappa_param[i])
                        return false;
                }
                return true;
            }

            bool operator!=(const param_type &other) const {
                return !(*this == other);
            }

        private:
            dim_vector mu_param;
            RealType kappa_param;
        };

        explicit von_mises_fisher(dim_vector mu, RealType kappa)
                : mu_param(mu), kappa_param(kappa), uniform_dist(), beta_dist((Dim - 1) / 2., (Dim - 1) / 2.) {
            static_assert(Dim >= 2);

        }

        explicit von_mises_fisher(const param_type &param)
                : von_mises_fisher(param.mu(), param.kappa()) {}

        void reset() {}

        param_type param() const {
            return param_type(mu(), kappa());
        }

        void param(const param_type &param) {
//            mu_gamma = gamma_dist_type(param.a());
//            kappa_gamma = gamma_dist_type(param.b());
        }

        template<typename URNG>
        result_type operator()(URNG &engine) {
            return generate(engine, mu_param, kappa_param);
        }

        template<typename URNG>
        result_type operator()(URNG &engine, const param_type &param) {
//            gamma_dist_type a_param_gamma(param.a()),
//                    b_param_gamma(param.b());
//            return generate(engine, a_param_gamma, b_param_gamma);
        }

//        result_type min() const { return 0.0; }
//
//        result_type max() const { return 1.0; }

        dim_vector mu() const { return mu_param; }

        RealType kappa() const { return kappa_param; }

//        bool operator==(const von_mises_fisher<Dim, RealType> &other) const {
//            if (!param() == other.param())
//                return false;
//            for (size_t i = 0; i < Dim; ++i) {
//                if (mu_param[i] != other.mu_param[i] ||
//                    kappa_param[i] == other.kappa_param[i])
//                    return false;
//            }
//            return true;
//        }

//        bool operator!=(const von_mises_fisher<Dim, RealType> &other) const {
//            return !(*this == other);
//        }

    private:
        std::uniform_real_distribution<> uniform_dist;
        std::normal_distribution<> normal_dist;
        vMF::beta_distribution<RealType> beta_dist;


        template<typename URNG>
        result_type generate(URNG &engine,
                             dim_vector &x_gamma,
                             RealType &y_gamma) {

            RealType w = _sample_weight(engine);

            dim_vector v = _sample_orthonormal_to(engine);

            dim_vector result;
            for (size_t i = 0; i < Dim; ++i) {
                result[i] = v[i] * std::sqrt(1. - std::pow(w, 2)) + w * mu_param[i];
            }

            return result;
        }

        dim_vector mu_param;
        RealType kappa_param;


        /**
         * Sample point on sphere orthogonal to mu
         *
         * @param mu
         * @return
         */
        template<typename URNG>
        result_type _sample_orthonormal_to(URNG &engine) {
//            result_type v;
//
//            for (size_t i = 0; i < Dim; ++i) {
//                v[i] = normal_dist(engine);
//            }
//
//            RealType dotted_result = 0;
//            for (size_t i = 0; i < Dim; ++i) {
//                dotted_result += mu_param[i] * v[i];
//            }
//
//            RealType mu_normed_result = 0;
//            for (size_t i = 0; i < Dim; ++i) {
//                mu_normed_result += mu[i] * mu[i];
//            }
//            mu_normed_result = std::sqrt(mu_normed_result);
//
//
//            result_type proj_mu_v;
//            for (size_t i = 0; i < Dim; ++i) {
//                proj_mu_v[i] = mu[i] * dotted_result / mu_normed_result;
//            }
//
//            result_type orthto;
//            RealType orthto_normed = 0;
//            for (size_t i = 0; i < Dim; ++i) {
//                orthto[i] = v[i] - proj_mu_v[i];
//                orthto_normed += orthto[i] * orthto[i];
//            }
//            orthto_normed = std::sqrt(orthto_normed);
//            // norm
//            for (size_t i = 0; i < Dim; ++i) {
//                orthto[i] /= orthto_normed;
//            }
//
//            return orthto;

            result_type v;
            for (size_t i = 0; i < Dim; ++i) {
                v[i] = normal_dist(engine);
            }

            RealType mu_dot_v = vector_operation::dot(mu_param, v);
            RealType mu_norm = vector_operation::norm(mu_param);

            result_type proj_mu_v;
            for (size_t i = 0; i < Dim; ++i) {
                proj_mu_v[i] = mu_param[i] * mu_dot_v / mu_norm;
            }

            result_type orthto = vector_operation::sub(v, proj_mu_v);
            RealType orthto_norm = vector_operation::norm(orthto);

            // norm
            for (size_t i = 0; i < Dim; ++i) {
                orthto[i] /= orthto_norm;
            }

            return orthto;
        }

        /**
         * Rejection sampling scheme for sampling distance from center on
         * surface of the sphere.
         *
         * @param kappa
         * @return
         */
        template<typename URNG>
        RealType _sample_weight(URNG &engine) {
            constexpr size_t dim = Dim - 1; // since S^{n-1}

            RealType b = dim / (std::sqrt(4. * std::pow(kappa_param, 2) + std::pow(dim, 2)) + 2 * kappa_param);
            RealType x = (1. - b) / (1. + b);
            RealType c = kappa_param * x + dim * std::log(1 - std::pow(x, 2));

            while (true) {
                RealType z = beta_dist(engine);
                RealType w = (1. - (1. + b) * z) / (1. - (1. - b) * z);
                RealType u = uniform_dist(engine);
                if (kappa_param * w + dim * std::log(1. - x * w) - c >= std::log(u))
                    return w;
            }
        }
    };

    template<typename CharT, size_t Dim, typename RealType>
    std::basic_ostream<CharT> &operator<<(std::basic_ostream<CharT> &os,
                                          const von_mises_fisher<Dim, RealType> &vmf) {
        os << "~vMF(" << vmf.mu() << "," << vmf.kappa() << ")";
        return os;
    }
//
//    template<typename CharT, typename RealType>
//    std::basic_istream<CharT> &operator>>(std::basic_istream<CharT> &is,
//                                          von_mises_fisher<RealType> &beta) {
//        std::string str;
//        RealType a, b;
//        if (std::getline(is, str, '(') && str == "~vMF" &&
//            is >> a && is.get() == ',' && is >> b && is.get() == ')') {
//            beta = von_mises_fisher<RealType>(a, b);
//        } else {
//            is.setstate(std::ios::failbit);
//        }
//        return is;
//    }

}

#endif //VMF_CPP_VMF_H
