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

#include <array>
#include <iostream>
#include <random>
#include <sstream>
#include <string>

#include "beta_distribution.h"
#include "vector_operation.h"

namespace vMF
{

    template <typename RealType = double>
    class von_mises_fisher
    {
    public:
        template <size_t Dim>
        using dim_vector = std::array<RealType, Dim>;

        template <size_t Dim>
        using result_t = dim_vector<Dim>;

        explicit von_mises_fisher()
        //          : uniform_dist(), beta_dist((dim_ - 1) / 2., (dim_ - 1) / 2.)
        {
        }

        template <typename URNG, size_t Dim>
        result_t<Dim> operator()(URNG &engine, const dim_vector<Dim> &mu, const RealType kappa)
        {
            return generate(engine, mu, kappa);
        }

    private:
        std::uniform_real_distribution<> uniform_dist;
        std::normal_distribution<> normal_dist;
        vMF::beta_distribution<RealType> beta_dist;

        template <typename URNG, size_t Dim>
        result_t<Dim> generate(URNG &engine, const dim_vector<Dim> &mu, const RealType kappa)
        //        result_t<Dim> generate(URNG &engine, const dim_vector<Dim> &mu, const RealType kappa)
        {
            RealType w = _sample_weight<URNG, Dim>(engine, kappa);
            //            RealType w = _sample_weight(engine, kappa, Dim);

            //            dim_vector<3> w;
            //            w[0] = _sample_weight<URNG, Dim>(engine, 3500);
            //            w[1] = _sample_weight<URNG, Dim>(engine, 3500);
            //            w[2] = _sample_weight<URNG, Dim>(engine, 3500);

            dim_vector<Dim> v = _sample_orthonormal_to(engine, mu);

            dim_vector<Dim> result;
            for (size_t i = 0; i < Dim; ++i)
            {
                result[i] = v[i] * std::sqrt(1. - std::pow(w, 2)) + w * mu[i];
                //                result[i] = v[i] * std::sqrt(1. - std::pow(w[i], 2)) + w[i] * mu[i];
            }
            //////////////////////
            //            RealType normed = vector_operation::norm(result);
            //            for (size_t i = 0; i < Dim; ++i)
            //            {
            //                result[i] /= normed;
            //            }
            //////////////////////

            return result;
        }

        /**
         * Sample point on sphere orthogonal to mu
         *
         * @param mu
         * @return
         */
        template <typename URNG, size_t Dim>
        result_t<Dim> _sample_orthonormal_to(URNG &engine, const dim_vector<Dim> &mu)
        {
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

            result_t<Dim> v;
            for (size_t i = 0; i < Dim; ++i)
            {
                v[i] = normal_dist(engine);
            }

            RealType mu_dot_v = vector_operation::dot(mu, v);
            RealType mu_norm = vector_operation::norm(mu);

            result_t<Dim> proj_mu_v;
            for (size_t i = 0; i < Dim; ++i)
            {
                proj_mu_v[i] = mu[i] * mu_dot_v / mu_norm;
            }

            result_t<Dim> orthto = vector_operation::sub(v, proj_mu_v);
            RealType orthto_norm = vector_operation::norm(orthto);

            // norm
            for (size_t i = 0; i < Dim; ++i)
            {
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
        template <typename URNG, size_t Dim>
        RealType _sample_weight(URNG &engine, const RealType kappa)
        {
            constexpr size_t dim = Dim - 1;  // since S^{n-1}

            RealType b = dim / (std::sqrt(4. * std::pow(kappa, 2) + std::pow(dim, 2)) + 2 * kappa);
            RealType x = (1. - b) / (1. + b);
            RealType c = kappa * x + dim * std::log(1 - std::pow(x, 2));

            while (true)
            {
                RealType z = beta_dist(engine);
                RealType w = (1. - (1. + b) * z) / (1. - (1. - b) * z);
                RealType u = uniform_dist(engine);
                if (kappa * w + dim * std::log(1. - x * w) - c >= std::log(u))
                    return w;
            }
        }
    };

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

}  // namespace vMF

#endif  // VMF_CPP_VMF_H
