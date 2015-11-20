/*
 * Timer.hh
 *
 *  Created on: Nov 13, 2015
 *      Author: brandon
 */

#ifndef CODE_TIMER_HH_
#define CODE_TIMER_HH_

#include <iostream>
#include <chrono>

class Timer
{
public:
    Timer() : beg_(clock_::now()) {}
    void reset() { beg_ = clock_::now(); }
    double elapsed() const {
        return std::chrono::duration_cast<second_>
            (clock_::now() - beg_).count(); }

private:
    typedef std::chrono::high_resolution_clock clock_;
    typedef std::chrono::duration<double, std::ratio<1> > second_;
    std::chrono::time_point<clock_> beg_;
};






#endif /* CODE_TIMER_HH_ */
