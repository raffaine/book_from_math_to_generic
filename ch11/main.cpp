#include "../common_utility.h"
#include "../common_math.h"

#include <iterator>
#include <list>
#include <stack>

// Concepts
template <typename F>
concept ForwardIterator = std::forward_iterator<F>;


// Exercises
template <ForwardIterator F, Integer N>
void reverse_inplace_n(F f, N n) {
    std::stack<F> s;
    for(int i=0; i < half(n); i++) {
        s.push(f++);
    }

    if(odd(n)) {
        f++;
    }

    while(!s.empty()) {
        std::swap(*(s.top()), *f++);
        s.pop();
    }
}

template <ForwardIterator F>
void reverse_inplace(F f, F l) {
    reverse_inplace_n(f, std::distance(f, l));
}

// Utils
using namespace std;
template <typename C>
requires requires (C c) {
    { begin(c) };
    { end(c) };
}
void print(const C& c) {
    for(auto f = begin(c); f != end(c); f++) {
        cout << *f << " ";
    }
    cout << endl;
}

int main(int argc, char** argv) {
    list l1 = {1,2,3,4,5,6};
    print(l1);
    reverse_inplace(l1.begin(), l1.end());
    print(l1);
    
    list l2 = {1,2,3,4,5,6};
    print(l2);
    reverse_inplace_n(l2.begin(), 4);
    print(l2);

    vector v1 = {2,4,6,8,10};
    print(v1);
    reverse_inplace(v1.begin(), v1.end());
    print(v1);

    return 0;
}