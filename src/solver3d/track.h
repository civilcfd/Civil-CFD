/*
 * track.h
 *
 * defines data structures to track the solution
 *
 */

#ifndef _TRACK_H
#define _TRACK_H

#include "solver_data.h"

struct track_data {
  double t; /* timestep */
  int n; /* integer representing simulation step.  can become non-consecutive if user deletes steps */
  
  struct track_data *next;
  struct track_data *first;

};

int track_add(double t);
int track_add_n(double t, int n);
int track_write();
int track_read();
int track_delete(int n);
int track_rewind();
int track_next();
double track_t();
int track_empty();

#endif
