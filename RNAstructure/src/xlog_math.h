#ifndef _XLOG_MATH_
#define _XLOG_MATH_

// Extended logarithmic/exponential functions from Mann 2006's work, numerically stable hidden markov model.
// Includes a log value for log of 0, implements summation, product for implementation of forward-backward algorithm.

#define LOG_OF_ZERO (-100000.0) // Log of zero.

// Convert probabilities into log space, with defaults base e.
double xlog(double prob);

double xexp(double log_value);

double xlog_sum(double log1, double log2);

double xlog_mul(double log1, double log2);

double xlog_div(double log1, double log2);

#endif // _XLOG_MATH_
