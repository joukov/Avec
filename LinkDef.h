/*********************************************************************
Copyright 2000 - 2006 Bruno Wittmer
This file is part of the Avec package.

Avec is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

Avec is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Avec; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
**********************************************************************/
#ifdef __CINT__

#define __AVECROOT__

#include <TArrayD.h>
#include <TH1.h>

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class Avec+;

#pragma link C++ class std::vector<Avec*>+;

#pragma link C++ function operator<<(std::ostream&, const Avec&);
#pragma link C++ function operator<<(std::ostream&, const std::vector<Avec*>&);
#pragma link C++ function operator&(const Avec&, const Avec&);

#pragma link C++ function operator|(Avec&, Avec&);
#pragma link C++ function operator|(Avec&, const vector<Avec*>&);
#pragma link C++ function operator|(const vector<Avec*>&, Avec&);

#pragma link C++ function operator+(double, const Avec&);
#pragma link C++ function operator-(double, const Avec&);
#pragma link C++ function operator*(double, const Avec&);
#pragma link C++ function operator/(double, const Avec&);

#pragma link C++ function avec_remove_nan(Avec&, dattype);
#pragma link C++ function avec_mask(const Avec&, const Avec&);
#pragma link C++ function vsum(const std::vector<dattype> &);
#pragma link C++ function vmax(const vector<dattype>&);
#pragma link C++ function vmin(const vector<dattype>&);
#pragma link C++ function lvmax(const vector<dattype>&);
#pragma link C++ function lvmin(const vector<dattype>&);

#pragma link C++ function vdiff(const Avec&);
#pragma link C++ function vdiff2(const Avec&);
#pragma link C++ function vint(const Avec&);
#pragma link C++ function vmean(const Avec&);
#pragma link C++ function vrms(const Avec&);
#pragma link C++ function avec_compress(Avec&, Avec&, Avec::size_type, bool);
#pragma link C++ function avec_compress(Avec&, Avec&, Avec&, Avec::size_type, bool);
#pragma link C++ function avec_compress(Avec&, std::vector<Avec*>, Avec::size_type, bool);
#pragma link C++ function avec_compress(Avec&, std::vector<Avec*> , std::vector<Avec*>, Avec::size_type, bool);
#pragma link C++ function avec_resize(const Avec&, Avec::size_type);
#pragma link C++ function rescale_size(const Avec&, double);
#pragma link C++ function avec_map(const Avec&, const Avec&, const Avec&);


  /********************/
  /* math.h functions */
  /********************/
#pragma link C++ function acos(const vector<dattype>&);
#pragma link C++ function asin(const vector<dattype>&);
#pragma link C++ function atan(const vector<dattype>&);
#pragma link C++ function atan2(const vector<dattype>&, const vector<dattype>&);
#pragma link C++ function ceil(const vector<dattype>&);
#pragma link C++ function cos(const vector<dattype>&);
#pragma link C++ function cosh(const vector<dattype>&);
#pragma link C++ function exp(const vector<dattype>&);
#pragma link C++ function fabs(const vector<dattype>&);
#pragma link C++ function floor(const vector<dattype>&);
#pragma link C++ function fmod(const vector<dattype>&, dattype);
#pragma link C++ function log(const vector<dattype>&);
#pragma link C++ function log10(const vector<dattype>&);
#pragma link C++ function pow(const vector<dattype>&, dattype);
#pragma link C++ function sin(const vector<dattype>&);
#pragma link C++ function sinh(const vector<dattype>&);
#pragma link C++ function sqrt(const vector<dattype>&);
#pragma link C++ function tan(const vector<dattype>&);
#pragma link C++ function tanh(const vector<dattype>&);

  /******************/
  /* Random Numbers */
  /******************/
#pragma link C++ function rand(vector<dattype>::size_type);
#pragma link C++ function rand(const vector<dattype>&);
#pragma link C++ function urand(const vector<dattype> &);
#pragma link C++ function grand(const vector<dattype> &, double);
#pragma link C++ function grand(Avec::size_type, double);

  /**********************/
  /* ASCII Input/Output */
  /**********************/
#pragma link C++ function vecwrite(const std::vector<dattype>&, const std::string&);
#pragma link C++ function vecwrite(const std::vector<Avec*>&, const std::string &);

#pragma link C++ function vecread(vector<dattype>&, const string &);
#pragma link C++ function vecread(const vector<Avec*>&, const string &);

#pragma link C++ function avec_get(const std::string&, TFile&);
#pragma link C++ function avec_get(const std::string&, const std::string&);

#pragma link C++ function TAD2AV(TArrayD &);
#pragma link C++ function TH12avec(TH1 &);
#pragma link C++ function avec2TH1D(Avec&, const std::string& , const std::string);

  /// DRAW_AVEC WILL BECOME OBSOLETE use avec_draw instead
//#pragma link C++ function draw_avec(Avec&, const std::string&, const std::string&, const std::string&, const std::string&);

//#pragma link C++ function draw_avec(const Avec&, const Avec&, const std::string&, const std::string&, const std::string&, const std::string&);

  //////////////////////
#pragma link C++ function avec_draw(const Avec&, const std::string&, const std::string&, const std::string&, const std::string&, int, const std::string&);

#pragma link C++ function avec_draw(const Avec&, const Avec&, const std::string&, const std::string&, const std::string&, const std::string&, int, const std::string&);

#pragma link C++ function avec_draw(const Avec&, const Avec&, const Avec&, const Avec&, const std::string&, const std::string&, const std::string&, const std::string&, int, const std::string&);

// Draw a set of Avecs against one Avec (const Avec)
#pragma link C++ function avec_draw(const Avec&, const std::vector<const Avec*>&, const std::string&, const std::string&, const std::string&, const std::string&, int, const std::string&);
#pragma link C++ function avec_draw(const Avec&, const std::vector<Avec*>&, const std::string&, const std::string&, const std::string&, const std::string&, int, const std::string&);
#pragma link C++ function avec_draw(const Avec&, const std::vector<Avec>&, const std::string&, const std::string&, const std::string&, const std::string&, int, const std::string&);

#pragma link C++ function avec_draw(const Avec&, const std::vector<const Avec*>&, const std::vector<const Avec*>&, const std::string&, const std::string&, const std::string&, const std::string&, int, const std::string&);

#pragma link C++ function avec_draw(const Avec&, const std::vector<Avec*>&, const std::vector<Avec*>&, const std::string&, const std::string&, const std::string&, const std::string&, int, const std::string&);

#pragma link C++ function avec_draw(const std::vector<const Avec*>&, const std::vector<const Avec*>&, const std::string&, const std::string&, const std::string&, const std::string&, int, const std::string&);

#pragma link C++ function avec_draw(const std::vector<Avec*>&, const std::vector<Avec*>&, const std::string&, const std::string&, const std::string&, const std::string&, int, const std::string&);

#pragma link C++ function avec_draw(const std::vector<Avec>&, const std::vector<Avec>&, const std::string&, const std::string&, const std::string&, const std::string&, int, const std::string&);

#pragma link C++ function avec_draw(const std::vector<Avec>&, const std::vector<Avec>&, const std::vector<Avec>&, const std::string&, const std::string&, const std::string&, const std::string&, int, const std::string&);

//#pragma link C++ function avec_draw(const std::vector<const Avec*>&, const std::vector<const Avec*>&, const std::vector<const Avec*>&, const std::string&, const std::string&, const std::string&, const std::string&, int, const std::string&);


#pragma link C++ function avec_plot(const Avec&, int, const std::string&, const std::string&);
#pragma link C++ function avec_plot_mask(const Avec&, const Avec&, int, const std::string&, const std::string&);
//#pragma link C++ function hfill(const Avec&, TH1& );
#pragma link C++ function hfill(const Avec&, TH1& , const Avec&);
#pragma link C++ function hfill(const Avec&, const Avec&, TH1& );
#pragma link C++ function hfill(const Avec&, const Avec&, TH2& , const Avec&);

#pragma link C++ function avec_plot(const Avec&, const Avec&, int, int, const std::string&, const std::string&);
#pragma link C++ function avec_plot_mask(const Avec&, const Avec&, const Avec&, int, int, const std::string&, const std::string&);

  /**********/
  /* Avec2D */
  /**********/
#pragma link C++ class std::vector<Avec>+;
#pragma link C++ class Avec2D+;
#pragma link C++ function Avec2D::size();
#pragma link C++ function Avec2D::operator[](vector<Avec>::size_type);
#pragma link C++ function Avec2D::push_back(const Avec &);

#pragma link C++ function operator<<(std::ostream&, const Avec2D&);

#pragma link C++ function operator&(const Avec2D&, const Avec2D&);

#pragma link C++ function vmax(const Avec2D&);
#pragma link C++ function vmaxc(const Avec2D&);
#pragma link C++ function vmaxr(const Avec2D&);
#pragma link C++ function vmin(const Avec2D&);
#pragma link C++ function vminc(const Avec2D&);
#pragma link C++ function vminr(const Avec2D&);
#pragma link C++ function vsum(const Avec2D&);
#pragma link C++ function vsumc(const Avec2D&);
#pragma link C++ function vsumr(const Avec2D&);
#pragma link C++ function vtrans(const Avec2D&);
#pragma link C++ function avec_resize(const Avec2D&, Avec2D::size_type, Avec::size_type);
//#pragma link C++ function rescale_size(const Avec2D&, double, double);


  /********************/
  /* math.h functions */
  /********************/
#pragma link C++ function acos(const Avec2D&);
#pragma link C++ function asin(const Avec2D&);
#pragma link C++ function atan(const Avec2D&);
#pragma link C++ function atan2(const Avec2D&, const Avec2D&);
#pragma link C++ function ceil(const Avec2D&);
#pragma link C++ function cos(const Avec2D&);
#pragma link C++ function cosh(const Avec2D&);
#pragma link C++ function exp(const Avec2D&);
#pragma link C++ function fabs(const Avec2D&);
#pragma link C++ function floor(const Avec2D&);
#pragma link C++ function fmod(const Avec2D&, dattype);
#pragma link C++ function log(const Avec2D&);
#pragma link C++ function log10(const Avec2D&);
#pragma link C++ function pow(const Avec2D&, dattype);
#pragma link C++ function sin(const Avec2D&);
#pragma link C++ function sinh(const Avec2D&);
#pragma link C++ function sqrt(const Avec2D&);
#pragma link C++ function tan(const Avec2D&);
#pragma link C++ function tanh(const Avec2D&);


#pragma link C++ function rand(const Avec2D&);
#pragma link C++ function urand(const Avec2D&);
#pragma link C++ function grand(const Avec2D&, double);

#pragma link C++ function vecwrite(const Avec2D&, const string&);
#pragma link C++ function vecread(Avec2D&, const string &);

#pragma link C++ function avec2d_get(const std::string&, TFile&);
#pragma link C++ function avec2d_get(const std::string&, const std::string&);

#pragma link C++ function avec_remove_nan(Avec2D&, dattype);

#pragma link C++ function avec_draw(const Avec&, const Avec2D&, const std::string&, const std::string&, const std::string&, const std::string&, int, const std::string&);
#pragma link C++ function avec_draw(const Avec2D&, const Avec2D&, const std::string&, const std::string&, const std::string&, const std::string&, int, const std::string&);
#pragma link C++ function avec_draw(const Avec2D&, const std::string&, const std::string&, const std::string&, const std::string&, int, std::string);

#pragma link C++ function hfill(const Avec2D&, TH1& , const Avec2D&);


#pragma link C++ global gROOT;
#pragma link C++ global gEnv;

#endif
