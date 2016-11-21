//
// Created by pkochetk on 11/20/16.
//

#ifndef SUPERCOMPUTING_CONDITION_H
#define SUPERCOMPUTING_CONDITION_H


class Condition {
public:
    virtual double F(double x, double y) = 0;
    virtual double fi(double x, double y) = 0;
};


#endif //SUPERCOMPUTING_CONDITION_H
