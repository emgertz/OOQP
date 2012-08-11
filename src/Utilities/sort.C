/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include <sortT.C>

template
void doubleLexSort<double>( int first[], int n, int second[], double data[] );

template
void doubleLexSort<float>( int first[], int n, int second[], float data[] );
