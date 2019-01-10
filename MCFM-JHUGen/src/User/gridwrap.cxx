//
//   @file    new_gridwrap.cxx         
//   
//
//   @author M.Sutton
// 
//   Copyright (C) 2013 M.Sutton (sutt@cern.ch)    
//
//   $Id: new_gridwrap.cxx, v0.0   Thu 18 Jul 2013 20:24:29 CEST sutt $

//namespace mcfm_bridge { 

static const int mxpart = 12;    // mcfm parameter : max number of partons in event record. defined in Inc/constants.f

/// function pointer hooks - set to 0 when no functions defined and applgrid not linked
void (*book_gridptr)()                            = 0;
void (*fill_gridptr)(const double evt[][mxpart])  = 0;
void (*write_gridptr)(double& )                   = 0;

//}

/// mcfm user hooks - do nothing if external function pointers are not assigned

/// book grids
extern "C" void  book_grid_()  {  if ( book_gridptr ) book_gridptr(); }
extern "C" void  book_grid__() {  book_grid_(); }

/// fill grids
extern "C" void  fill_grid_(const double evt[][mxpart])     {   if ( fill_gridptr ) fill_gridptr(evt); }
extern "C" void  fill_grid__(const double evt[][mxpart])    {   fill_grid_(evt); }

/// write grids
extern "C" void  write_grid_(double& xstotal)      { if ( write_gridptr ) write_gridptr(xstotal); }
extern "C" void  write_grid__(double& xstotal)     { write_grid_(xstotal); }



