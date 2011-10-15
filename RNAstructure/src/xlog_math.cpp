#include <math.h>
#include "xlog_math.h"
#include <stdio.h>
#include <stdlib.h>

// Convert probabilities into log space, with defaults base e.
double xlog(double value)
{
	if(value == 0)
	{
		return(LOG_OF_ZERO);
	}
	else if (value > 0)
	{
		return(log(value));
	}
	else
	{
		printf("log of a negative number @ %s(%d): %.6f", __FILE__, __LINE__, value);
		int* p_null = NULL;
		*p_null = 0;
		exit(0);
	}
}

double xexp(double log_value)
{
	if(log_value == LOG_OF_ZERO)
	{
		return(0);
	}
	else
	{
		return(exp(log_value));
	}
}

double xlog_sum(double log1, double log2)
{
	if(log1 == LOG_OF_ZERO)
	{
		return(log2);
	}
	else if(log2 == LOG_OF_ZERO)
	{
		return(log1);
	}
	else
	{
		if(log1 > log2)
		{
			return( log1 + xlog(1 + xexp(log2-log1)) );
		}
		else
		{
			return( log2 + xlog(1 + xexp(log1-log2)) );
		}
	}
}

double xlog_mul(double log1, double log2)
{
	if(log1 == LOG_OF_ZERO || log2 == LOG_OF_ZERO)
	{
		return(LOG_OF_ZERO);
	}
	else
	{
		return(log1 + log2);
	}
}

double xlog_div(double log1, double log2)
{
	if(log2 == LOG_OF_ZERO)
	{
		printf("Division by zero error at %s(%d)\n", __FILE__, __LINE__);
		exit(0);
	}
	else
	{
		return(xlog_mul(log1, -1.0 * log2));
	}
}
