//
// Created by pkochetk on 11/20/16.
//

#ifndef SUPERCOMPUTING_VAR6COND_H
#define SUPERCOMPUTING_VAR6COND_H

#include "Condition.h"
class Var6Cond: public Condition {
    double F(double x, double y) {
        return 4*(2-3*x*x-3*y*y);
    }

    double fi(double x, double y) {
        return (1-x*x)*(1-x*x) + (1-y*y)*(1-y*y);
    }
};

#endif //SUPERCOMPUTING_VAR6COND_H
