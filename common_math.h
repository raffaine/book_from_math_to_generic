#include <iostream>
#include <limits>
#include <concepts>
#include <functional>
#include <initializer_list>

template <typename T>
concept Integer = std::integral<T>;

template <typename T>
concept Number = std::integral<T> || std::floating_point<T>;

template <Integer T>
bool odd(T n) {
    return n & 1;
}

template <Integer T>
T half(T n) {
    return n >> 1;
}

template <typename T>
bool compare(const T& a, const T& b) {
    return a == b;
}

template <>
bool compare(const float& x, const float& y) {
    const float eps = std::numeric_limits<float>::epsilon()*100; //nudges epsilon a little since powers are really high
    if (y == 0.) return x == 0.;
    auto diff = x / y;
    diff = (diff < 1.) ? 1. - diff : diff - 1.;
    return diff <= eps;
}

template <typename T>
T identity_element(std::plus<T>) {
    return T(0);
}
template <typename T>
T identity_element(std::multiplies<T>) {
    return T(1);
}

template <typename T>
T inverse(const T& a, std::plus<T>) {
    return -a;
}

template <typename T>
T inverse(const T& x, std::multiplies<T>) {
    return 1 / x;
}

////////////////////
////// Concepts
////////////////////

template <typename T>
concept Regular = std::regular<T>;

template <typename Op, typename T>
concept BinaryOperation = Regular<T> && std::invocable<Op, T, T>;

template <typename Op, typename T>
concept SemigroupOperation = BinaryOperation<Op, T> && requires (Op op, T a, T b, T c) {
    op(a, op(b, c)) == op(op(a, b), b);
};

template <typename Op, typename T>
concept MonoidOperation = SemigroupOperation<Op, T> && requires (Op op, T a) {
    op(a, identity_element(op)) == a;
    op(identity_element(op), a) == a;
};

template <typename Op, typename T>
concept GroupOperation = MonoidOperation<Op, T> && requires (Op op, T a) {
    op(a, inverse(a, op)) == identity_element(op);
    op(inverse(a, op), a) == identity_element(op);
};

template <typename T>
concept InputIterator = std::input_iterator<T>;

template <typename T>
concept Semiring = std::regular<T>;

////////////////////
////// Algorithms
////////////////////

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

template <Regular A, Integer N, typename Op>
requires MonoidOperation<Op, A>
A power_monoid(A a, N n, Op op) {
    // Precondition n > 0
    if (n == 0) return identity_element(op);
    return power_semigroup(a, n, op);
}

template <Regular A, Integer N, typename Op>
requires GroupOperation<Op, A>
A power_group(A a, N n, Op op) {
    if (n < 0) {
        n = -n;
        a = inverse(a, op);
    }
    return power_monoid(a, n, op);
}

template <Regular A, Integer N>
A power(A a, N n) {
    return power_group(a, n, std::multiplies<A>());
}

/////////////////////////
///  New Types
/////////////////////////

// 2x2 Matrices (Uses a Union to better compact accessors)
union Matrix2x2 {
    // Offers two ways of accessing the same data, either component-wise or array form
    struct {
        int a11, a12, a21, a22;
    } c;
    int arr[4];

    // Operations required so Matrix2x2 is Regular (Default Constructible, Copyable and Equality Check)
    Matrix2x2(int d = 0) {
        c.a11 = c.a22 = d;
        c.a12 = c.a21 = 0;
    }

    Matrix2x2(int a11, int a12, int a21, int a22) {
        c.a11 = a11;
        c.a12 = a12;
        c.a21 = a21;
        c.a22 = a22;
    }

    Matrix2x2(const int elem[4]) {
        for (int i=0; i<4; i++)
            arr[i] = elem[i];
    }

    Matrix2x2(const Matrix2x2& M) : Matrix2x2(M.arr) {}

    bool operator==(const Matrix2x2& A) const {
        for (int i=0; i<4; i++) {
            if (arr[i] != A.arr[i]) return false;
        }
        return true;
    }

    //TODO: Not implemented yet
    Matrix2x2 inverse() const {
        return *this;
    }

    // Operation required so that Matrix2x2 can be negated and therefore is a Monoid over Addition
    friend Matrix2x2 operator*(int s, const Matrix2x2& A) {
        return Matrix2x2(s * A.c.a11, s * A.c.a12, s * A.c.a21, s * A.c.a22);
    }

    friend Matrix2x2 operator-(const Matrix2x2& A) {
        return -1 * A;
    }

    friend Matrix2x2 operator/(int s, const Matrix2x2& A) {
        return s * A.inverse();
    }

    // Makes the Matrix2x2 valid for Multiplicative Operations
    friend Matrix2x2 operator*(const Matrix2x2& A, const Matrix2x2& B) {
        return Matrix2x2(
            A.c.a11 * B.c.a11 + A.c.a12 * B.c.a21,
            A.c.a11 * B.c.a12 + A.c.a12 * B.c.a22,
            A.c.a21 * B.c.a11 + A.c.a22 * B.c.a21,
            A.c.a21 * B.c.a21 + A.c.a22 * B.c.a22
        );
    }

    // Makes the Matrix2x2 valid for Additive Operations
    friend Matrix2x2 operator+(const Matrix2x2& A, const Matrix2x2& B) {
        return Matrix2x2(
            A.c.a11 + B.c.a11,
            A.c.a12 + B.c.a12,
            A.c.a21 + B.c.a21,
            A.c.a22 + B.c.a22
        );
    }
};

// Polynomials
template <Number T>
class Polynomial {
    template <Number U>
    friend class Polynomial;

public:
    Polynomial() = default;
    Polynomial(std::initializer_list<T> args) : coeficients(std::move(args)) {}
    
    template <Number U>
    requires std::constructible_from<T, U>
    Polynomial(Polynomial<U> p) : coeficients(p.coeficients.size()) {
        std::copy(std::begin(p.coeficients), std::end(p.coeficients), std::begin(coeficients));
    }

    template <typename U>
    requires Semiring<U>
    U evaluate(const U& x) {
        auto first = std::begin(coeficients);
        const auto last = std::end(coeficients);
        if (first == last) return U(0);
        U sum(*first);
        while (++first != last) {
            sum *= x;
            sum += *first;
        }
        return sum;
    }

private:
    std::vector<T> coeficients;
};

// Complex Numbers
template <Number T>
class complex {
public:
    complex() = default;
    complex(T a,  T b = T(0)) : a(a), b(b) {}
    complex(const complex<T>& z) : a(z.a), b(z.b) {}

    bool operator==(const complex<T>& z) const {
        return a == z.a && b == z.b;
    }

    complex<T> conjugate() const { return complex(a, -b); }
    friend complex<T> operator~(const complex<T>& z) { return z.conjugate(); }

private:
    T a, b;
};

// Gauss Integers are Complex numbers with integer coeficients
template <Integer T>
using gauss_integer = complex<T>;


// N dimensional Square Matrices
template <typename T, int N>
requires Regular<T> && requires (T a, T b) {
    a * b;
    a + b;
}
class SqMatrix {
public:
    SqMatrix() {
        a.fill(T());
    }
    SqMatrix(T diag_value) : SqMatrix() {
        for (int i = 0; i < N; i++) {
            a[i * N + i] = diag_value;
        }
    }

    SqMatrix(std::initializer_list<T> vs) {        
        std::copy(std::begin(vs), std::end(vs), std::begin(a));
    }
    
    template <Number U>
    requires std::constructible_from<T, U>
    SqMatrix(const SqMatrix<U, N>& M) {
        std::copy(std::begin(M.a), std::end(M.a), std::begin(a));
    }

    friend bool operator==(const SqMatrix<T, N> A, const SqMatrix<T, N> B) {
        for (int i = 0; i < N * N; i++) {
            if (A.a[i] != B.a[i]) return false;
        }
        return true;
    }

    friend SqMatrix<T, N> operator*(T s, const SqMatrix<T, N>& M) {
        SqMatrix<T,N> R(M);
        std::transform(std::begin(M.a), std::end(M.a), std::begin(R.a), [s](const T& v) {
            return s * v;
        });
        return R;
    }

    // TODO: Not Implemented as it's not needed yet
    SqMatrix<T, N> inverse() const {
        return *this; 
    }

    friend SqMatrix<T, N> transpose(const SqMatrix<T, N>& M) {
        SqMatrix<T, N> R;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                R.a[i * N + j] = M.a[j * N + i];
            }
        }
        return R;
    }

    friend SqMatrix<T, N> operator/(T s, const SqMatrix<T, N>& M) {
        return s * M.inverse();
    }
    
    friend SqMatrix<T, N> operator*(const SqMatrix<T, N>& A, const SqMatrix<T, N>& B) {
        SqMatrix<T,N> C;
        
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                const int ij = i * N + j;
                for (int k = 0; k < N; k++) {
                    C.a[ij] += A.a[i * N + k] * B.a[k * N + j];
                }
            }
        }

        return C;
    }

    friend std::ostream& operator<<(std::ostream& os, SqMatrix<T, N> M) {
        for (int i = 0; i < N; i++) {
            os << "| ";
            for (int j = 0; j < N; j++) {
                os << +M.a[i * N + j] << " ";
            }
            os << "|" << std::endl;
        }
        return os;
    }

    T operator()(int i, int j) {
        const int ij = i * N + j;
        return (i >= 0 && j >= 0 && ij < a.size())? a[ij] : T(0);
    }
    
    void operator()(int i, int j, T v) {
        const int ij = i * N + j;
        if (i >= 0 && j >= 0 && ij < a.size()) {
            a[ij] = v;
        }
    }

private:
    std::array<T, N * N> a;
};



//
// Calculate Fibonacci using a convenient Matrix2x2 Representation for its calculation
//            v(n) = exp(M, (n-1))*v(1) 
//    where v(n) is a vector containing the nth value and its predecessor {fib(n), fib(n-1)}
//     and M is a 2x2 square Matrix2x2
//
int mat_fibonacci(int n) {
    std::pair<int, int> v = {1, 0};
    if (n > 0) {
        auto M = power(Matrix2x2({1,1,1,0}), n - 1);
        v.first = M.c.a11;
    }

    return v.first;
}
