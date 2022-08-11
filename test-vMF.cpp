#include <stdio.h>

#include <cmath>

#include "include/vMF/beta_distribution.h"
#include "soraxas_toolbox/main.h"
#include "vMF/vMF.h"
#include "vMF/vMF_constexpr.h"

using namespace sxs;

template <typename Name>
void test_print()
{
    // String as template parameter
    std::cout << Name::c_str();
}
// using seq = sequence<MACRO_GET_STR("Hello world!")>;

#include <type_traits>
#include <typeindex>
#include <unordered_map>
#include <unordered_set>

#include "soraxas_toolbox/compile_time_string.h"

// class CompileTimeDict {
// public:
//
//
//
//     template<typename CompileTimeString>
//     double &of() {
//         auto search_result = already_processed_ts.find(typeid(CompileTimeString));
//         if (search_result != already_processed_ts.end()) {
////            return 0;
////            std::cout << "Already processed " << typeid(CompileTimeString).name() << std::endl;
//        } else {
//            already_processed_ts[typeid(CompileTimeString)] = 0.;
////            std::cout << "Processing " << typeid(CompileTimeString).name() << "... \n";
////            return 1;
//        }
//        return already_processed_ts[typeid(CompileTimeString)];
////        return std::string(CompileTimeString::c_str());
//    }
//
//    template<typename CompileTimeString>
//    decltype(auto) increment() {
//        auto &&val = of<CompileTimeString>();
//        val += 1;
//        return val;
//    }
//
//    std::unordered_map<std::type_index, double> already_processed_ts;
//};

///**
// * A public facing compile time dict
// *
// * @tparam Tag
// */
// template<typename Tag>
// class NumericCompileTimeDict : public CompileTimeAnyTypeDict<Tag> {
// public:
//    using tag = Tag;
//
//
//    template<typename CompileTimeString, typename T=double>
//    static T &of() {
//        _CompileTimeDict_RegisterClass<Tag, CompileTimeString>::doRegister();
//        static T thing;
//        return thing;
//    }
//
//
////    template<typename CompileTimeString>
////    static double &off() {
////        RegisterClass<Tag, CompileTimeString>::doRegister();
////        auto storage = get_double_storage();
////        storage
////        return thing;
////    }
//
//    template<typename CompileTimeString>
//    static decltype(auto) increment() {
//        auto &&val = of<CompileTimeString>();
//        val += 1;
//        return val;
//    }
//
//    static auto &get_double_storage() {
//        static std::unordered_map<std::string, double> double_storage_;
//        return double_storage_;
//    }
//
////    static void print_dict() {
////        std::cout << "===== Dict: " << tag::c_str() << " =====" << std::endl;
////        for (auto &&key: keys()) {
////            std::cout << "- " << key << std::endl;
////        }
////        std::cout << "===== ===== ===== =====" << std::endl;
////    }
//
// private:
//
//
//};

#include <typeindex>
#include <typeinfo>

#include "soraxas_toolbox/compile_time_dict.h"
#include "soraxas_toolbox/globals.h"

int main()
{
    //
    //    using numdict = NumericCompileTimeDict<CT_STR("storage") >;
    //
    //    numdict::

    // std::unordered_set<std::type_index>()

    using numdict = sxs::CompileTimeMappingTypeDict<CT_STR("yes storage")>;

    numdict::print_dict();

    numdict::get_map<std::type_index>();

    numdict::of<CT_STR("ok")>() += 1;
    numdict::increment<CT_STR("ok")>();
    numdict::increment<CT_STR("ok")>();

    numdict::increment<CT_STR("ok"), int>();

    numdict::of<CT_STR("ok"), int>() += 1;
    numdict::of<CT_STR("ok"), int>() += 1;

    numdict::of<CT_STR("ok"), std::string>() = "ha";

    numdict::print_dict();

    using dict = sxs::CompileTimeAnyTypeDict<CT_STR("storage")>;

    dict::print_dict();

    dict::of<CT_STR("one key"), std::string>();

    dict::print_dict();
    println("====");

    //    println(dict::of<CT_STR("ha") >());
    for (auto &&key : dict::keys())
    {
        println(key);
    }

    dict::of<CT_STR("name"), std::string>() = "HAHA this is a stored string";

    println(dict::of<CT_STR("name"), std::string>());

    dict::of<CT_STR("my stat"), double>() = 1;
    dict::of<CT_STR("my stat"), double>() += 25;
    println(dict::of<CT_STR("my stat"), double>());

    double &lets_keep_track_of_this = dict::of<CT_STR("my stat"), double>();
    lets_keep_track_of_this *= -1;
    println(dict::of<CT_STR("my stat"), double>());

    //
    //    dict yes{};
    //    CompileTimeDict<CT_STR("dict 1") > yes2{};
    //
    //    println(yes.of<CT_STR("qq") >());
    //    println(yes.increment<CT_STR("qq") >());
    //    println(yes.increment<CT_STR("qq") >());
    //    println(yes.increment<CT_STR("qq") >());
    //    println(yes.increment<CT_STR("qq2") >());
    ////    println(yes.tkeys);
    //
    //    println(yes2.of<CT_STR("qq") >());
    //    prin

    //    for (auto &&key: CompileTimeDict<CT_STR("dict 1") >::keys) {
    //        println(key);
    //    }

    //
    //    for (auto &&key: yes.keys()) {
    //        println(key);
    //    }

    //    auto &&d = get<CT_STR("dict"), int>();
    //    println(d);
    //    d += 2;
    //    println(d);
    //    println(get<CT_STR("dict"), int>());
    //
    //
    //    constexpr auto myArray{[]() constexpr {
    //        std::array<double, 5> result{};
    //        for (int i = 0; i < 5; ++i) {
    //            result[i] = i * 4;
    //        }
    //        return result;
    //    }()};
    //
    ////    test_print<CT_STR("Hello World!") >();
    //
    //    auto dict = CompileTimeDict{};
    //
    //    println(dict.increment<CT_STR("Yes") >());
    //    println(dict.increment<CT_STR("Yes") >());
    //    println(dict.increment<CT_STR("Yes") >());

    return 1;

    println("hi");

    std::random_device r;
    std::seed_seq seed2{r(), r(), r(), r(), r(), r(), r(), r()};
    std::mt19937 e2(seed2);

    println();

    constexpr double kappa = 100;

    vMF::von_mises_fisher_constexpr<3> vmfc{{2, 2, 5}, kappa};
    vMF::von_mises_fisher vmf{};

    println("vmf(e2)");
    for (int i = 0; i < 100; ++i)
    {
        std::array<double, 3> mu = {2, 2, 5};
        println_spaced(vmfc(e2), vmf(e2, mu, kappa), ",");
    }
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