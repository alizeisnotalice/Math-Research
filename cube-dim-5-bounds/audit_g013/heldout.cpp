#include "round012_graph_moments_gmp.hpp"
#include <iostream>
int main(){
 round012::Engine engine("","",10.0);
 round012::Counts k{{5,0,5,0,5,0,0,5,0,0}};
 std::cout<<engine.expectation(k).get_str()<<"\n";
}
