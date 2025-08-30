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

// Small Utility for comparison of Floats
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

template <>
bool compare(const double& x, const double& y) {
    const double eps = std::numeric_limits<double>::epsilon()*100; //nudges epsilon a little since powers are really high
    if (y == 0.) return x == 0.;
    auto diff = x / y;
    diff = (diff < 1.) ? 1. - diff : diff - 1.;
    return diff <= eps;
}

// Partial Specialization for a Requirement for Monoids T (e, an identity element over operations +, *)
template <typename T>
T identity_element(std::plus<T>) {
    return T(0);
}

template <typename T>
T identity_element(std::multiplies<T>) {
    return T(1);
}

// Partial Specialization for Requirement for Groups T (a^-1, where a Op a^-1 == e specialized for Op in [+,*])
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

// A Type T that has an Equality Comparison, a Default Constructor and Copy Constructor (or assignment)
template <typename T>
concept Regular = std::regular<T>;

// An operation Op over two elements from type T, so that given a, b in T, Op(a,b) also produces a T (closed)
template <typename Op, typename T>
concept BinaryOperation = Regular<T> && std::invocable<Op, T, T> && requires (Op op, T a, T b) {
    { op(a, b) } -> std::same_as<T>;
};

// A BinaryOperation Op that is also associative over its elements
template <typename Op, typename T>
concept AssociativeOperation = BinaryOperation<Op, T> && requires (Op op, T a, T b, T c) {
    op(a, op(b, c)) == op(op(a, b), b);
};

// An Associative operation Op over elements in Semigroup T
template <typename Op, typename T>
concept SemigroupOperation = AssociativeOperation<Op, T>;

// A SemigroupOperation Op over elements in Monoid T have an identity element e
template <typename Op, typename T>
concept MonoidOperation = SemigroupOperation<Op, T> && requires (Op op, T a) {
    op(a, identity_element(op)) == a;
    op(identity_element(op), a) == a;
};

// A MonoidOperation Op over elements in Group T have an inverse element a^-1
template <typename Op, typename T>
concept GroupOperation = MonoidOperation<Op, T> && requires (Op op, T a) {
    op(a, inverse(a, op)) == identity_element(op);
    op(inverse(a, op), a) == identity_element(op);
};

// A commuatative operation is that which doesn't depend on the order of the operand
template <typename Op, typename T>
concept Commutative = BinaryOperation<Op, T> && requires (Op op, T a, T b) {
    op(a, b) == op(b, a);
};

// A GroupOperation Op over elements in Abelian Group T are commutative
template <typename Op, typename T>
concept AbelianGroupOperation = GroupOperation<Op, T> && Commutative<Op, T>;

// Types T and U that admits an Addition Operation between elements T and U
template <typename T, typename U = T>
concept Additive = requires (T a, U b) {
    a + b;
    b + a;
};

// Types T and U that admits a Multiplication Operation between elements T and U
template <typename T, typename U = T>
concept Multiplicative = requires (T a, U b) {
    a * b;
    b * a;
};

// A Monoid T where the operation is Addition
template <typename T>
concept AdditiveMonoid = MonoidOperation<std::plus<T>, T>;

// A Monoid T where the operation is Multiplication
template <typename T>
concept MultiplicativeMonoid = MonoidOperation<std::multiplies<T>, T>;

// A Group T where the operation is Addition
template <typename T>
concept AdditiveGroup = GroupOperation<std::plus<T>, T>;

// A Group T where the operation is Multiplication
template <typename T>
concept MultiplicativeGroup = GroupOperation<std::multiplies<T>, T>;

// Given an Additive type T that is Multiplicative with Type U, the Multiplication distributes over the Addition
template <typename T, typename U = T>
concept Distributive = Additive<T> && Multiplicative<T, U> && requires (U x, T y, T z) {
    // Distributive Property
    x * (y + z) == (x * y) + (x * z);
    (y + z) * x == (y * x) + (z * x);
};

// Given T is a Monoid over Addition (commuative) and Multiplication, it's a SemiRing if they distribute over another
template <typename T>
concept Semiring =  Commutative<std::plus<T>, T> && Distributive<T> &&
                    AdditiveMonoid<T> && MultiplicativeMonoid<T> && 
requires (T x) {
    // There are two different elements T(1) and T(0)
    T(1) != T(0);
    // T(0) works as the identity for the addition
    x + T(0) == x;
    T(0) + x == x;
    // T(1) works as the identity for the multiplication
    x * T(1) == x;
    T(1) * x == x;
    // T(0) when multiplied by any element results in T(0)
    x * T(0) == T(0);
    T(0) * x == T(0);
};

// A Semiring T where the elements can be inversed over addition (which means they form a group)
template <typename T>
concept Ring = Semiring<T> && AdditiveGroup<T>;

// A Ring T with a Commutative Multiplication over its elements
template <typename T>
concept CommutativeRing = Ring<T> && Commutative<std::multiplies<T>, T>;

// A Commutative Ring T with no zero divisors
template <typename T>
concept IntegralDomain = CommutativeRing<T> && requires (T a, T b) {
    // i.e. there is no a, b in T that when multiplied result in 0, unless either a or b are zero
    a != T(0) && b != T(0) && (a * b) != T(0);
};

template <Integer I>
I quotient(I a, I b) {
    return a / b;
}

template <Integer I>
I remainder(I a, I b) {
    return a % b;
}

using Natural = int;
template <Integer I>
Natural norm(I a) {
    return std::abs(a);
}

template <typename T>
concept EuclidianDomain = IntegralDomain<T> && requires (T a, T b) {    
    // Quotient q and Remainder r where for any integer a and non-zero b where a == q * b + r
    b != T(0) && a == quotient(a, b) * b + remainder(a, b);
    // Norm
    a == T(0) && norm(a) == T(0);
    a != T(0) && norm(a) != T(0);
    // Norm of the product of non-zero elements is bigger or equal to the norm of the elements
    b != T(0) && norm(a * b) >= norm(a);
    // Norm of the remainder of a and b is less than norm of b
    norm(remainder(a,b)) < norm(b);
};

// An Integral Domain T where every non-zero element is invertible
template <typename T>
concept Field = IntegralDomain<T> && requires (T a) {
    a != T(0) && inverse(std::multiplies<T>(), a) * a == T(1);
};

// An Additive Group G and a Ring R where R distributes over elements of G
template <typename G, typename R>
concept Module = AdditiveGroup<G> && Ring<R> && Distributive<G, R>;

// A Module on G and R where R is a Field
template <typename G, typename R>
concept VectorSpace = Module<G, R> && Field<R>;

////////////////////
////// Algorithms
////////////////////

// Given that for an Euclidian Domain E,
//    norm(remainder(a,b)) < norm(b)
//  the function below is guaranteed to terminate (on each iteration b will have a smaller norm until it inevitably reaches zero)
template <EuclidianDomain E>
E gcd(E a, E b) {
    while (b != E(0)) {
        a = remainder(a, b);
        std::swap(a, b);
    }
    return a;
}

// Technically the type A doesn't require equality or default constructible, maybe Semiregular?
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
    // Precondition n >= 0
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

// Partial specialization that applies the power_group over the Multiplication Operation (*) over A
template <Regular A, Integer N>
A power(A a, N n) {
    return power_group(a, n, std::multiplies<A>());
}

/////////////////////////
///  New Types
/////////////////////////

// 2x2 Matrices (Uses a Union to better compact accessors and faster operations with unrolled loops)
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


// A Helper to allow specialization of the Divide function (for the integer case)
template <Number T>
struct divide_impl;

// Complex Numbers
template <Number T>
class complex {
    friend struct divide_impl<T>;
public:
    complex() = default;
    complex(T a,  T b = T(0)) : a(a), b(b) {}
    complex(const complex<T>& z) : a(z.a), b(z.b) {}

    bool operator==(const complex<T>& z) const {
        return a == z.a && b == z.b;
    }

    complex<T> conjugate() const { return complex(a, -b); }
    friend complex<T> operator~(const complex<T>& z) { return z.conjugate(); }
    
    friend complex<T> operator+(const complex<T>& z) { return z; }
    friend complex<T> operator-(const complex<T>& z) { return complex<T>(-z.a, -z.b); }

    friend complex<T> operator-(const complex<T>& x, const complex<T>& y) {
        return x + (-y);
    }

    friend complex<T> operator+(const complex<T>& x, const complex<T>& y) {
        return complex<T>(x.a + y.a, x.b + y.b);
    }

    friend complex<T> operator*(const complex<T>& x, const complex<T>& y) {
        return complex<T>(x.a * y.a - x.b * y.b, x.a * y.b + x.b * y.a);
    }

    friend complex<T> operator/(const complex<T>& x, const complex<T>& y) {
        return divide_impl<T>::divide(x, y);
    }

    friend std::ostream& operator<<(std::ostream& os, const complex<T>& z) {
        os << +z.a << " + " << +z.b << "i";
        return os;
    }

private:
    T a, b;
};

// Avoid dividing Twice if Type is floating point (so that a/b == a * (1/b))
template <Number T>
struct divide_impl {
    static complex<T> divide(const complex<T>& x, const complex<T>& y) {
        const T d( T(1) / (y.a * y.a + y.b * y.b));
        return complex<T>(d * (x.a * y.a + x.b * y.b), d * (x.b * y.a + x.a * y.b));
    }
};

// Gauss Integers are Complex numbers with integer coeficients
template <Integer I>
using gauss_integer = complex<I>;

// TODO: Not sure this is how quotient is calculated for Gauss Integers
// At least performs an integer division (the above optimization doesn't apply since 1/b == 0 for integers)
template <Integer I>
struct divide_impl<I> {
    static complex<I> divide(const complex<I>& x, const complex<I>& y) {
        const I d(y.a * y.a + y.b * y.b);
        return complex<I>((x.a * y.a + x.b * y.b) / d, (x.b * y.a + x.a * y.b) / d);
    }
};

// TODO: Not really sure how to analitically calculate this
template <Integer I>
complex<I> remainder(const complex<I>& a, const complex<I>& b) {
    return complex<I>(1);
}
    

// N dimensional Square Matrices
template <typename T, int N>
requires Regular<T> && Additive<T> && Multiplicative<T>
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
