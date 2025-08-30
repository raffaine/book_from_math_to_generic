#include <iomanip>
#include <stack>
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

// Used to Represent Distances that can either be a positive integer or inifinite representing no connection
template <typename T>
requires AdditiveGroup<T> && MultiplicativeGroup<T>
struct Dist {
    Dist(bool infinite, T value) : v((infinite) ? T(0) : value), inf(infinite) {}
    Dist() : Dist(true, T(0)) {}
    Dist(T value) : Dist(false, value) {}
    Dist(const Dist<T>& o) : inf(o.inf), v(o.v) {}

    friend bool operator==(const Dist& d1,const Dist& d2) {
        return (d1.inf && d2.inf) || (d1.v == d2.v);
    }
    
    bool operator<(const Dist& d2) {
        return (d2.inf && !inf) || (!d2.inf && inf) || (v < d2.v);
    }
    
    bool operator<=(const Dist& d2) {
        return *this < d2 || *this == d2;
    }
    
    friend Dist<T> inverse(const Dist<T>& x, multiplies<Dist<T>>) {
        return Dist(x.v == 0, (x.inf)? 0 : inverse(x.v));
    }
    
    friend Dist<T> operator*(const Dist<T>& a, const Dist<T>& b) {
        return Dist(a.inf || b.inf, a.v + b.v);
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
        if (d.inf) os << "  \u221E";
        else os << setw(3) << +d.v;
        return os;
    }

    bool inf = true;
    T v = 0;
};

int main(int argc, char** argv) {
    // General Tests
    cout << Polynomial({1,2,3,4,5}).evaluate(1) << endl;
    cout << ((((complex(2.) + complex(0., 2.)) / complex(2., 0.)) == ((complex<double>(1.) - complex(0., 1.)) * complex(0., 1.))) ? "Equal" : "Different") << endl;
    cout << gauss_integer(2)/gauss_integer(2) << endl;
    cout << gcd(148, 36) << endl;

    // Social Nets Exercise
    SqMatrix<bool, friends.size()> mb(true);
    setv(mb, {{"Ari", {"Bev", "Don"}}, {"Bev", {"Fay"}}, {"Cal", {"Don"}}, {"Don", {"Fay"}}, {"Eva", {"Gia"}}});

    print(mb, "Initial friends to ");
    
    auto mf = power(mb, friends.size()-1);
    print(mf, "All socially connected people to ");

    // Shortest Distance Exercise
    SqMatrix<Dist<int>, places.size()> mp(0);
    mp(0, 1, 6); mp(0, 3, 3);
    mp(1, 4, 2); mp(1, 5, 10);
    mp(2, 0, 7);
    mp(3, 2, 5); mp(3, 5, 4);
    mp(4, 6, 3);
    mp(5, 2, 6); mp(5, 4, 7); mp(5, 6, 8);
    mp(6, 1, 9);

    // Calculate the Shortest Distance by calculating above matrix power to n-1
    const auto n = places.size();
    auto mpf = power(mp, n-1);
    // For visual aid, prints the resulting matrix
    cout << mpf << endl;

    // Until user requests to stop, print the distance between two given cities in the set
    do {
        cout << endl << "Given the following Cities (in this order as rows and columns): ";
        for (auto p : places) {
            cout << p << " ";
        }
        cout << endl << " and City Distances given as " << endl << mp << endl;
        string v1, v2;
        cout << "city 1 (ctrl+c to terminate): ";
        cin >> v1;
        cout << "city 2: ";
        cin >> v2;
        
        const auto p1 = find(begin(places), end(places), v1);
        const auto p2 = find(begin(places), end(places), v2);
        if (p1 == p2 || p1 == end(places) || p2 == end(places)) continue;

        Dist shortest = mpf(p1 - begin(places), p2 - begin(places));
        if (shortest == Dist(true, 0)) {
            cout << "There is no connection between the 2 places." << endl;
            continue;
        }
        cout << "Shortest distance between: " << v1 << " and " << v2 << " is ";
        cout << shortest << endl << endl;
        
        // Search the resulting matrix for the element with minimum value for a given column (the best way to that node)
        // starting with destination and stopping once we found the initial node (use a stack to reverse it later)
        stack<int> path;
        path.push(p2-begin(places));
        while(path.top() != (p1 - begin(places))) {
            int pos = 0;
            for(int i=0; i < places.size(); i++) {
                auto cur = mpf(i, path.top());
                if (i != path.top() && cur < shortest) {
                    shortest = cur;
                    pos = i;
                }
            }
            
            path.push(pos);
        }

        cout << "It takes the following path" << endl;
        // Print path on the reverse order
        while(!path.empty()) {
            cout << places[path.top()] << " ";
            path.pop();
        }
        cout << endl;

    } while(true);

    return 0;
}
