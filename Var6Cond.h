//
// Created by pkochetk on 11/20/16.
//

#ifndef SUPERCOMPUTING_VAR6COND_H
#define SUPERCOMPUTING_VAR6COND_H

#include "Condition.h"
#include <math.h>

class Var6Cond: public Condition {
    double F(double x, double y) {
        return (x*x+y*y)*sin(x*y);
    }

    double fi(double x, double y) {
        return 1+sin(x*y);
    }
};

#endif //SUPERCOMPUTING_VAR6COND_H
