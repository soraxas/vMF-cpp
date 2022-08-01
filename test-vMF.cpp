#include <stdio.h>

#include <cmath>

#include "include/vMF/beta_distribution.h"
#include "soraxas_toolbox/main.h"
#include "vMF/vMF.h"

using namespace sxs;

template <size_t Dim>
double sample_weight(double kappa)
{
    static_assert(Dim >= 2);

    constexpr size_t dim = Dim - 1;  // since S^{n-1}
    double b = dim / (std::sqrt(4. * std::pow(kappa, 2) + std::pow(dim, 2)) + 2 * kappa);
    double x = (1. - b) / (1. + b);
    double c = kappa * x + dim * std::log(1 - std::pow(x, 2));

    std::random_device r;
    std::seed_seq seed2{r(), r(), r(), r(), r(), r(), r(), r()};
    std::mt19937 e2(seed2);

    std::uniform_real_distribution<> uniform{};
    vMF::beta_distribution param{dim / 2., dim / 2.};

    vMF::beta_distribution beta{param};

    while (true)
    {
        double z = beta(e2);
        double w = (1. - (1. + b) * z) / (1. - (1. - b) * z);
        double u = uniform(e2);
        if (kappa * w + dim * std::log(1. - x * w) - c >= std::log(u))
            return w;
    }
}

int main()
{
    println("hi");
    print(sample_weight<2>(5));

    std::random_device r;
    std::seed_seq seed2{r(), r(), r(), r(), r(), r(), r(), r()};
    std::mt19937 e2(seed2);

    println();

    vMF::von_mises_fisher<3> vmf{{2, 2, 5}, 10000000000000};

    println("vmf(e2)");
    for (int i = 0; i < 100; ++i)
        println(vmf(e2), ",");
    //    println(vmf(e2));

    // def _sample_weight(kappa, dim):
    //     """Rejection sampling scheme for sampling distance from center on
    //     surface of the sphere.
    //     """
    //     dim = dim - 1  # since S^{n-1}
    //     b = dim / (np.sqrt(4. * kappa**2 + dim**2) + 2 * kappa)
    //     x = (1. - b) / (1. + b)
    //     c = kappa * x + dim * np.log(1 - x**2)
    //
    //     while True:
    //         z = np.random.beta(dim / 2., dim / 2.)
    //         w = (1. - (1. + b) * z) / (1. - (1. - b) * z)
    //         u = np.random.uniform(low=0, high=1)
    //         if kappa * w + dim * np.log(1. - x * w) - c >= np.log(u):
    //             return w

    println("done");
}