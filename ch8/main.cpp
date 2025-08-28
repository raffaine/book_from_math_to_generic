#include <iomanip>
#include "../common_utility.h"
#include "../common_math.h"

using namespace std;

constexpr array friends = {"Ari","Bev","Cal","Don","Eva","Fay","Gia"};
constexpr array places  = {"A","B","C","D","E","F","G"};

template <typename T, typename V>
void setv(T& M, const char* n1, const char* n2, V val) {
    const auto bf(begin(friends)), ef(end(friends));
    setv(M, find(bf, ef, n1) - bf, find(bf, ef, n2) - bf, val);
}

template <typename T>
void setv(T& M, int i, int j, bool val) {
    M(i,j,val);
    M(j,i,val);
}

template <typename T>
void setv(T& M, vector<pair<const char*, vector<const char*>>> vs) {
    for (auto [p, ps] : vs ) {
        for(auto f : ps) {
            setv(M, p, f, true);
        }
    }
}

template <typename T>
void print(T mb, const char* message) {
    for (int i = 0; i < friends.size(); i++) {
        cout << message << friends[i] << ": ";
        for (int j = 0; j < friends.size(); j++) {
            if(i != j && mb(i, j)) cout << friends[j] << " ";
        }
        cout << endl;
    }
}

template <typename T>
struct Dist {
    Dist() = default;
    Dist(T nv, bool i = false) : v(nv), inf(i) {}
    Dist(const Dist<T>& o) : inf(o.inf), v(o.v) {}

    friend bool operator==(const Dist& d1,const Dist& d2) {
        return (d1.inf && d2.inf) || (d1.v == d2.v);
    }
    
    friend Dist<T> inverse(const Dist<T>& x, multiplies<Dist<T>>) {
        return Dist((x.inf || x.v == 0)? 0 : 1 / x.v, x.v == 0);
    }
    
    friend Dist<T> operator*(const Dist<T>& a, const Dist<T>& b) {
        return Dist(a.v + b.v, a.inf || b.inf);
    }

    friend Dist<T> operator+(const Dist<T>& a, const Dist<T>& b) {
        Dist<T> c(a);
        c += b;
        return c;
    }

    const Dist<T>& operator+=(const Dist<T>& b) {
        if (inf) {
            inf = inf && b.inf;
            v = inf?  0 : b.v;
        } else if (!b.inf) {
            v = min(v, b.v);
        }

        return *this;
    }
    
    friend Dist<T> operator+(const Dist<T>& d) {
        return d;
    }
    friend ostream& operator<<(ostream& os, Dist<T> d) {
        if (d.inf) os << "inf";
        else os << setw(3) << +d.v;
        return os;
    }

    bool inf = true;
    T v = 0;
};

int main(int argc, char** argv) {
    Polynomial({1,2,3,4,5}).evaluate(1);

    SqMatrix<bool, friends.size()> mb(true);
    setv(mb, {{"Ari", {"Bev", "Don"}}, {"Bev", {"Fay"}}, {"Cal", {"Don"}}, {"Don", {"Fay"}}, {"Eva", {"Gia"}}});

    print(mb, "Initial friends to ");
    
    auto mf = power(mb, friends.size()-1);
    print(mf, "All socially connected people to ");

    SqMatrix<Dist<int>, places.size()> mp(0);
    for (int i = 0; i < places.size(); i++) {
        for (int j = 0; j < places.size(); j++) {
            if (i == j) continue;
            mp(i,j, Dist(0, true));
        }
    }
    mp(0, 1, 6);
    mp(0, 3, 3);
    mp(1, 4, 2);
    mp(1, 5, 10);
    mp(2, 0, 7);
    mp(3, 2, 5);
    mp(3, 5, 4);
    mp(4, 6, 3);
    mp(5, 2, 6);
    mp(5, 4, 7);
    mp(5, 6, 8);
    mp(6, 1, 9);

    cout << mp << endl;

    // cout << Dist(0, true) << " * " << Dist(0) << " = " << Dist(0, true) * Dist(0)  << endl;
    // cout << Dist(0, true) << " * " << Dist(1) << " = " << Dist(0, true) * Dist(1)  << endl;
    // cout << Dist(0, true) << " * " << Dist(0, true) << " = " << Dist(0, true) * Dist(0, true)  << endl;
    // cout << Dist(0, true) << " + " << Dist(0, true) << " = " << Dist(0, true) + Dist(0, true)  << endl;
    // cout << Dist(0, true) << " + " << Dist(0) << " = " << Dist(0, true) + Dist(0)  << endl;
    // cout << Dist(0, true) << " + " << Dist(1) << " = " << (Dist<int>(0, true) * Dist(0)) + Dist(1)  << endl;
    // cout << Dist(0) << " + " << Dist(0) << " = " << Dist(0) + Dist(0)  << endl;
    // cout << Dist(0) << " + " << Dist(2) << " = " << Dist(0) + Dist(2)  << endl;

    // SqMatrix<Dist<int>, places.size()> mpf;
    // const auto N = places.size();
    // for (int i = 0; i < N; i++) {
    //     for (int j = 0; j < N; j++) {
    //         for (int k = 0; k < N; k++) {
    //             cout << mp(i,k) << " * " << mp(k,j) << " (" << (mp(i,k) * mp(k,j)) << ") + " ;
    //             mpf(i, j, mpf(i, j) + (mp(i,k) * mp(k,j)));
    //         }
    //         cout << " = " << mpf(i, j) << endl;
    //     }
    // }
    // cout << mp * mp;
    const auto n = places.size();
    auto mpf = power(mp, n-1);
    cout << mpf << endl;
    cout << mpf * transpose(mp) << endl;

    return 0;
}
