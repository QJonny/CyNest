/*
    Copyright (C) 2012 The NEST Initiative
    This file is part of NEST.
*/

#include "config.h"

#ifndef GSL_BINOMIAL_RANDOMDEV_H
#define GSL_BINOMIAL_RANDOMDEV_H

#include "randomgen.h"
#include "randomdev.h"
#include "lockptr.h"
#include "dictdatum.h"
#include "gslrandomgen.h"

#ifdef HAVE_GSL

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


/*BeginDocumentation
Name: rdevdict::gsl_binomial - GSL binomial random deviate generator
Description:
This function returns a random integer from the binomial distribution,
the number of successes in n independent trials with probability
p. The probability distribution for binomial variates is,

   p(k) = (n! / k!(n-k)!) p^k (1-p)^(n-k)  , 0<=k<=n, n>0

Please note that the RNG used to initialize gsl_binomial has to be
from the GSL (prefixed gsl_ in rngdict)

Parameters:
   p - probability of success in a single trial (double)
   n - number of trials (positive integer)

SeeAlso: CreateRDV, RandomArray, rdevdict
Author: Jochen Martin Eppler
*/


namespace librandom {

/**
 Class GSL_BinomialRandomDev

 Generates an RNG which returns Binomial(k;p;n)
 distributed random numbers out of an RNG which returns
 binomially distributed random numbers:

    p(k) = (n! / k!(n-k)!) p^k (1-p)^(n-k)  , 0<=k<=n, n<0

 Arguments:
  - pointer to an RNG
  - parameter p (optional, default = 0.5)
  - parameter n (optional, default = 1)

 @see http://www.gnu.org/software/gsl/manual/html_node/The-Binomial-Distribution.html
 @ingroup RandomDeviateGenerators
*/

class GSL_BinomialRandomDev : public RandomDev
  {
  public:
    // accept only lockPTRs for initialization,
    // otherwise creation of a lock ptr would
    // occur as side effect---might be unhealthy
    GSL_BinomialRandomDev(RngPtr, double p_s = 0.5, unsigned int n_s=1);
    GSL_BinomialRandomDev(double p_s = 0.5, unsigned int n_s=1);

    /**
     * set parameters for p and n
     * @parameters
     * p - success probability for single trial
     * n - number of trials
     */
    void   set_p_n   (double, unsigned int);
    void   set_p     (double);                 //!<set p
    void   set_n     (unsigned int);           //!<set n

    /**
     * Import sets of overloaded virtual functions.
     * We need to explicitly include sets of overloaded
     * virtual functions into the current scope.
     * According to the SUN C++ FAQ, this is the correct
     * way of doing things, although all other compilers
     * happily live without.
     */
    using RandomDev::uldev;

    unsigned long uldev();         //!< draw integer
    unsigned long uldev(RngPtr) const;   //!< draw integer, threaded
    bool has_uldev() const { return true; }

    double operator()();           //!< return as double
    double operator()(RngPtr) const;     //!< return as double, threaded

    //! set distribution parameters from SLI dict
    void set_status(const DictionaryDatum&);

    //! get distribution parameters from SLI dict
    void get_status(DictionaryDatum&) const;

  private:
    double p_;                  //!<probability p of binomial distribution
    unsigned int n_;            //!<parameter n in binomial distribution

    gsl_rng *rng_;
  };

  inline
  double GSL_BinomialRandomDev::operator()()
  {
    return static_cast<double>(uldev());
  }

  inline
  double GSL_BinomialRandomDev::operator()(RngPtr rthrd) const
  {
    return static_cast<double>(uldev(rthrd));
  }

}

#endif

#endif

