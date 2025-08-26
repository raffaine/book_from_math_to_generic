#include <iostream>
#include <chrono>
#include <random>
#include <functional>
#include <limits>
#include <concepts>
#include <string>

// Simple timer to clock the implementations

#define PRINT_TIME_AS_MS 1 // In case measures are more familiar in milliseconds

using namespace std::chrono;

class Timer {
    public:
        Timer() : _start(high_resolution_clock::now()) {}

        microseconds clock() {
            auto end = high_resolution_clock::now();
            auto diff = end - _start;
            _start = end;
            return duration_cast<microseconds>(diff);
        }
private:
    high_resolution_clock::time_point _start;    
};

std::ostream& operator<<(std::ostream& os, microseconds t) {
#if PRINT_TIME_AS_MS
    os << duration_cast<milliseconds>(t).count() << "ms";
#else
    os << t.count() << " microseconds";
#endif
    return os;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, std::plus<T>) {
    os << "x";
    return os;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, std::multiplies<T>) {
    os << "^";
    return os;
}

////// Math functions

template <typename T>
concept Integer = std::integral<T>;

template <Integer T>
bool odd(T n) {
    return n & 1;
}

template <Integer T>
T half(T n) {
    return n >> 1;
}

template <typename T>
concept Regular = std::regular<T>;

template <typename Op, typename T>
concept BinaryOperation = Regular<T> && std::invocable<Op, T, T>;

template <typename Op, typename T>
concept SemigroupOperation = BinaryOperation<Op, T> && requires (Op op, T a, T b, T c) {
    op(a, op(b, c)) == op(op(a, b), b);
};

template <Regular A, Integer N, typename Op>
requires SemigroupOperation<Op, A>
A power_accumulate_semigroup(A r, A a, N n, Op op) {
    // A precondition for correctness is that n >= 0
    if (n == 0) return r;
    while(true) {
        if (odd(n)) {
            r = op(r, a);
            if (n == 1) return r;
        }
        n = half(n);
        a = op(a, a);
    }
}

template <Regular A, Integer N, typename Op>
requires SemigroupOperation<Op, A>
A power_semigroup(A a, N n, Op op) {
    // Precondition n > 0
    while (!odd(n)) {
        a = op(a, a);
        n = half(n);
    }
    if (n == 1) return a;
    return power_accumulate_semigroup(a, op(a, a), half(n - 1), op);
}

template <typename T>
T identity_element(std::plus<T>) {
    return T(0);
}
template <typename T>
T identity_element(std::multiplies<T>) {
    return T(1);
}

template <typename Op, typename T>
concept MonoidOperation = SemigroupOperation<Op, T> && requires (Op op, T a) {
    op(a, identity_element(op)) == a;
    op(identity_element(op), a) == a;
};

template <Regular A, Integer N, typename Op>
requires MonoidOperation<Op, A>
A power_monoid(A a, N n, Op op) {
    // Precondition n > 0
    if (n == 0) return identity_element(op);
    return power_semigroup(a, n, op);
}

template <typename T>
T inverse(const T& a, std::plus<T>) {
    return -a;
}

template <typename T>
T inverse(const T& x, std::multiplies<T>) {
    return T(1) / x;
}

template <typename Op, typename T>
concept GroupOperation = MonoidOperation<Op, T> && requires (Op op, T a) {
    op(a, inverse(a, op)) == identity_element(op);
    op(inverse(a, op), a) == identity_element(op);
};

template <Regular A, Integer N, typename Op>
requires GroupOperation<Op, A>
A power_group(A a, N n, Op op) {
    if (n < 0) {
        n = -n;
        a = inverse(a, op);
    }
    return power_monoid(a, n, op);
}

template <typename T>
bool compare(const T& a, const T& b) {
    return a == b;
}

template <>
bool compare(const float& x, const float& y) {
    const float eps = std::numeric_limits<float>::epsilon()*100;
    if (y == 0.) return x == 0.;
    auto diff = x / y;
    diff = (diff < 1.) ? 1. - diff : diff - 1.;
    return diff <= eps;
}

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

// Uses a Union to better compact accessors
union Matrix {
    // Offers two ways of accessing the same data, either component-wise or array form
    struct {
        int a11, a12, a21, a22;
    } c;
    int arr[4];

    // Operations required so Matrix is Regular (Default Constructible, Copyable and Equality Check)
    Matrix() = default;

    Matrix(int a11, int a12, int a21, int a22) {
        c.a11 = a11;
        c.a12 = a12;
        c.a21 = a21;
        c.a22 = a22;
    }

    Matrix(const int elem[4]) {
        for (int i=0; i<4; i++)
            arr[i] = elem[i];
    }

    Matrix(const Matrix& M) : Matrix(M.arr) {}

    bool operator==(const Matrix& A) const {
        for (int i=0; i<4; i++) {
            if (arr[i] != A.arr[i]) return false;
        }
        return true;
    }

    // Operation required so that Matrix can be negated and therefore is a Monoid over Addition
    friend Matrix operator*(int s, const Matrix& A) {
        return Matrix(s * A.c.a11, s * A.c.a12, s * A.c.a21, s * A.c.a22);
    }
};

namespace std {
    // Makes the Matrix valid for Multiplicative Operations
    template <>
    struct multiplies<Matrix> {
        Matrix operator()(const Matrix& A, const Matrix& B) const {
            return Matrix(
                A.c.a11 * B.c.a11 + A.c.a12 * B.c.a21,
                A.c.a11 * B.c.a12 + A.c.a12 * B.c.a22,
                A.c.a21 * B.c.a11 + A.c.a22 * B.c.a21,
                A.c.a21 * B.c.a21 + A.c.a22 * B.c.a22
            );
        }
    };    
    // Makes the Matrix valid for Additive Operations
    template <>
    struct plus<Matrix> {
        Matrix operator()(const Matrix& A, const Matrix& B) const {
            return Matrix(
                A.c.a11 + B.c.a11,
                A.c.a12 + B.c.a12,
                A.c.a21 + B.c.a21,
                A.c.a22 + B.c.a22
            );
        }
    };
}

// Identity Element for Additive Monoids
template <>
Matrix identity_element(std::plus<Matrix>) {
    return Matrix({0,0,0,0});
}

// Identity Element for Multiplicative Monoids
template <>
Matrix identity_element(std::multiplies<Matrix>) {
    return Matrix({1,0,0,1});
}

// Inverse Element for Additive Groups
template <>
Matrix inverse(const Matrix& A, std::plus<Matrix>) {
    return -1 * A;
}

//TODO: Identity Element for Multiplicative Groups
template <>
Matrix inverse(const Matrix& A, std::multiplies<Matrix>) {
    return A; // need to calculate inverse
}

//
// Calculate Fibonacci using a convenient Matrix Representation for its calculation
//            v(n) = exp(M, (n-1))*v(1) 
//    where v(n) is a vector containing the nth value and its predecessor {fib(n), fib(n-1)}
//     and M is a 2x2 square Matrix
//
int mat_fibonacci(int n) {
    std::pair<int, int> v = {1, 0};
    if (n > 0) {
        Matrix M = power_group(Matrix({1,1,1,0}), n - 1, std::multiplies<Matrix>());
        v.first = M.c.a11;
    }

    return v.first;
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
