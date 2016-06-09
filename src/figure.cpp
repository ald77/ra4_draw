#include "figure.hpp"

#include "utilities.hpp"

using namespace std;

Figure::FigureComponent::FigureComponent(const Figure &figure,
                                         const shared_ptr<Process> &process):
  figure_(figure),
  process_(process){
}
