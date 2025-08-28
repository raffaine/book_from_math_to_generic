#include <random>

#include "../common_utility.h"
#include "../common_math.h"

template <typename T, Integer N, typename Op, typename R = T>
void run_operation(T min_bound, T max_bound, N min_bound2, N max_bound2, std::function<T(N,T)> base_op, const int NUM = 1000000) {
    // Seed with a real random value, if available
    std::random_device r;
 
    // Choose a random mean between min_bound and max_bound
    std::default_random_engine e(r());
    std::uniform_int_distribution<R> u1dist(min_bound, max_bound);
    std::uniform_int_distribution<N> u2dist(min_bound2, max_bound2);

    // Seed the input vectors with a fresh set of random numbers
    std::vector<T> input(NUM);
    std::vector<N> input2(NUM);
    for (int i=0; i < NUM; i++) {
        input[i] = u1dist(e);
        input2[i] = u2dist(e);
    }

    // Create 2 output vectors, one for each approach
    std::vector<T> output(NUM), output2(NUM);

    // Time N operations using our custom approach
    Timer t;
    for (int i=0; i < NUM; i++) {
        output2[i] = (input[i] == T(0))? T(0) : power_group(input[i], input2[i], Op());
    }
    auto custom_duration = t.clock();
    
    // Time N multiplications using the hardware implementation
    for (int i=0; i < NUM; i++) {
        output[i] = base_op(input[i], input2[i]);
    }
    auto base_duration = t.clock();
    
    std::cout << "Custom Operation took " << custom_duration << std::endl;    
    std::cout << "Base Operation took " << base_duration << std::endl;
    std::cout << "Custom runs " << ((double) base_duration.count()) / custom_duration.count() << " times faster." << std::endl;
    auto failed_ind = std::invoke([&] {
        for (int i=0; i < NUM; i++) {
            if (!compare(output[i], output2[i])) return i;
        }
        return -1;
    });

    if (failed_ind < 0 ) {
        std::cout << "Results match" << std::endl;
    } else {
        std::cout << "Results don't match on i=" << failed_ind << ", " << +output[failed_ind] << " != " << +output2[failed_ind] << std::endl;
        std::cout << "Values multiplied: " << +input[failed_ind] << Op() << +input2[failed_ind] << std::endl;
    }
}

int8_t my_multiply(int8_t a, int b) {
    return static_cast<int8_t>(b) * a;
}

template <typename T, Integer N>
T my_power(T a, N n) {
    if (a == T(0)) return a;
    if (n < 0) {
        a = inverse(a, std::multiplies<T>());
        n = -n;
    }
    T res = T(1);
    while(n-- > 0) {
        res = std::multiplies<T>()(res, a);
    }
    return res;
}

template <typename T, typename Op>
void runFibs(T min_bound, T max_bound, Op base_op, Op custom_op, const int NUM=1000000) {
    // Seed with a real random value, if available
    std::random_device r;
 
    // Choose a random mean between min_bound and max_bound
    std::default_random_engine e(r());
    std::uniform_int_distribution<T> udist(min_bound, max_bound);

    // Seed the input vectors with a fresh set of random numbers
    std::vector<T> input(NUM);
    for (int i=0; i < NUM; i++) {
        input[i] = udist(e);
    }

    // Create 2 output vectors, one for each approach
    std::vector<T> output(NUM), output2(NUM);

    // Time N multiplications using the hardware implementation
    Timer t;
    for (int i=0; i < NUM; i++) {
        output[i] = base_op(input[i]);
    }
    auto base_duration = t.clock();
    
    // Time N operations using our custom approach
    for (int i=0; i < NUM; i++) {
        output2[i] = custom_op(input[i]);
    }
    auto custom_duration = t.clock();
    
    
    std::cout << "Custom Operation took " << custom_duration << std::endl;    
    std::cout << "Base Operation took " << base_duration << std::endl;
    std::cout << "Custom runs " << ((double) base_duration.count()) / custom_duration.count() << " times faster." << std::endl;
    auto failed_ind = std::invoke([&] {
        for (int i=0; i < NUM; i++) {
            if (!compare(output[i], output2[i])) return i;
        }
        return -1;
    });

    if (failed_ind < 0 ) {
        std::cout << "Results match" << std::endl;
    } else {
        std::cout << "Results don't match on i=" << failed_ind << ", " << +output[failed_ind] << " != " << +output2[failed_ind] << std::endl;
        std::cout << "Value for Nth Fibonacci, with n=" << +input[failed_ind] << std::endl;
    }
}

int fibonacci_iterative(int n) {
    if (n == 0) return 0;
    std::pair<int, int> v = {0, 1};
    for (int i = 1; i < n; ++i) {
        v = {v.second, v.first + v.second};
    }
    return v.second;
}

int main(int argc, char** argv) {

    std::cout << "Results for 64bit signed integer" << std::endl;
    run_operation<int64_t, int64_t, std::plus<int64_t>>(-500000, 500000, -500000, 500000, std::multiplies<int64_t>());
    std::cout << std::endl;
    std::cout << "Results for 32bit signed integer" << std::endl;
    run_operation<int, int, std::plus<int>>(-5000, 5000, -5000, 5000, std::multiplies<int>());
    std::cout << std::endl;
    std::cout << "Results for 32bit unsigned integer" << std::endl;
    run_operation<uint32_t, uint32_t, std::plus<uint32_t>>(1000U, 10000U, 1000U, 10000U, std::multiplies<uint32_t>());
    std::cout << std::endl;
    std::cout << "Results for 8bit signed integer" << std::endl;
    run_operation<int8_t, int, std::plus<int8_t>>(-10, 10, -10, 10, my_multiply);

    std::cout << std::endl;
    std::cout << "Results for 32bit signed float with 32bit integer exponents" << std::endl;
    run_operation<float, int, std::multiplies<float>, int>(-15, 15, -10, 10, my_power<float, int>);

    std::cout << std::endl;
    std::cout << "Results for 32bit signed integer Fibonacci" << std::endl;
    runFibs(10, 100, fibonacci_iterative, mat_fibonacci);



    return 0;
}
