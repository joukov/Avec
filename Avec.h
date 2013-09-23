/*********************************************************************
Copyright 2000 - 2006 Bruno Wittmer
          2006 Michael Henke
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

#ifndef AVEC_H
#define AVEC_H

#include <vector>
#include <fstream>
#include <string>
#include <iostream>

// Root and CINT definitions
#ifdef __CINT__
#define __AVECROOT__
#endif

// Definitions for Windows
#ifdef WIN32
#ifdef WINAVEC_EXPORTS
#define WINAVEC_API __declspec(dllexport)
#else
#define WINAVEC_API __declspec(dllimport)
#endif

// this template function is not available in windows...
template <class T>
T min(T a,T b)
{return (a < b ? a : b);}
#else
#define WINAVEC_API 
#endif

// Forward declaration of 2D-Avec
class Avec2D;

//!  Data type used for the Analysis vector.
/*!
  This was deliberately not implemented as template to simplify interactive use with CINT.
*/
typedef double dattype;

// Helping Class used for Avec ranges
class AvecRange
{
 public:
  AvecRange(std::vector<dattype>::size_type f, std::vector<dattype>::size_type l):first(f),last(l){}
  std::vector<dattype>::size_type first;
  std::vector<dattype>::size_type last;
};

//!  Analysis Vector class. 
/*!
  This class inherits from std::vector<dattype>. This class is useable with CINT and ROOT.\n
Everything that std::vector is able to do, also works for this class.
*/
#ifdef __AVECROOT__
#include <TNamed.h>
#endif
class WINAVEC_API Avec : public std::vector<dattype>
#ifdef __AVECROOT__
,  public TNamed
#endif
{
 public:
  explicit Avec(size_t newsize=0, dattype init = 0) : std::vector<dattype>(newsize,init)  {}
  Avec(iterator first, iterator last) : std::vector<dattype>(first,last){}
  Avec(const_iterator first, const_iterator last) : std::vector<dattype>(first,last){}
  Avec(size_t newsize, dattype low, dattype high);
  Avec(const Avec2D & source);
  Avec(const std::vector<float> & source);
  Avec(size_t newsize, const dattype data[]);
  virtual ~Avec(){}

  static void help(const std::string s="");

  void print(std::ostream &out = std::cout, const char* sep="\n") const;
  friend std::ostream& operator<< (std::ostream& o, const Avec& v);

  Avec operator()(size_t first, size_t last) const;
  //Avec operator[](AvecRange r) const;

  //void operator=(const std::vector<float>& source);

  Avec operator+(const Avec& source) const;
  Avec operator-(const Avec& source) const;
  Avec operator*(const Avec& source) const;
  Avec operator/(const Avec& source) const;
  Avec operator%(const Avec& source) const;

  Avec operator+(double offset) const;
  Avec operator-(double offset) const;
  Avec operator*(double factor) const;
  Avec operator/(double factor) const;
  Avec operator%(long int factor) const;

  friend Avec operator+(double offset, const Avec& source);
  friend Avec operator-(double offset, const Avec& source);
  friend Avec operator*(double factor, const Avec& source);
  friend Avec operator/(double factor, const Avec& source);
  friend Avec operator%(double factor, const Avec& source);

  Avec operator-() const;

  void operator+=(double offset);
  void operator-=(double offset);
  void operator*=(double factor);
  void operator/=(double factor);
  void operator%=(long int factor);

  void operator+=(const Avec& source);
  void operator-=(const Avec& source);
  void operator*=(const Avec& source);
  void operator/=(const Avec& source);
  void operator%=(const Avec& source);

  Avec operator<(const Avec& source) const;
  Avec operator>(const Avec& source) const;
  Avec operator<=(const Avec& source) const;
  Avec operator>=(const Avec& source) const;
  Avec operator!=(const Avec& source) const;
  Avec operator==(const Avec& source) const;
  
  Avec operator<(const dattype value) const;
  Avec operator>(const dattype value) const;
  Avec operator<=(const dattype value) const;
  Avec operator>=(const dattype value) const;
  Avec operator!=(const dattype value) const;
  Avec operator==(const dattype value) const;

  Avec operator&&(const Avec& source) const;
  Avec operator||(const Avec& source) const;

  // Overloaded Methods of TObject
  virtual void	Draw(Option_t* option = "");
  virtual void Browse(TBrowser* b);

  // Static members to manipulate the general behavior of the class
  static int VERBOSE_FLAG; //! If zero, suppress output
  static int OUTPUT_PRECISION; //! Precision for text output
  static int OUTPUT_WIDTH; //! Not used yet
  static std::vector< std::pair < int, int> > MGRAPH_PALETTE;
  
  ClassDef(Avec,1)
};

// operators
WINAVEC_API Avec operator&(const Avec& vec1, const Avec& vec2);

//WINAVEC_API vector<const Avec*> operator|(const Avec& first,const Avec& second);
//std::vector<const Avec*> operator|(const std::vector<const Avec*>& first, const Avec& second);
//std::vector<const Avec*> operator|(const Avec& first, const std::vector<const Avec*>& second);

WINAVEC_API std::vector<Avec*> operator|(Avec& first, Avec& second);
WINAVEC_API std::vector<Avec*> operator|(Avec& first, const std::vector<Avec*>& second);
WINAVEC_API std::vector<Avec*> operator|(const std::vector<Avec*>& first, Avec& second);

///////////////////////////////////
// Some functions of general use //
///////////////////////////////////
WINAVEC_API void avec_remove_nan(Avec& data, dattype replacement = 0);
WINAVEC_API Avec avec_mask(const Avec& data, const Avec& mask);
WINAVEC_API dattype vsum(const std::vector<dattype> &vec);
WINAVEC_API dattype vmax(const std::vector<dattype> &vec);
WINAVEC_API dattype vmin(const std::vector<dattype> &vec);
WINAVEC_API size_t lvmax(const std::vector<dattype> &vec);
WINAVEC_API size_t lvmin(const std::vector<dattype> &vec);
WINAVEC_API Avec vdiff(const Avec& v1);
WINAVEC_API Avec vdiff2(const Avec& v1);
WINAVEC_API Avec vint(const Avec& v1);
WINAVEC_API dattype vmean(const Avec& v1);
WINAVEC_API dattype vrms(const Avec& v1);
// Compress Avecs to a smaller number of entries by averaging over several original entries
WINAVEC_API void avec_compress(Avec& x, Avec& y, Avec::size_type n, bool remove_empty = false);
WINAVEC_API void avec_compress(Avec& x, Avec& y, Avec& err_y, Avec::size_type n, bool remove_empty = false);
WINAVEC_API void avec_compress(Avec& x, std::vector<Avec*> yv, Avec::size_type n, bool remove_empty = false);
WINAVEC_API void avec_compress(Avec& x, std::vector<Avec*> yv, std::vector<Avec*> err_yv, Avec::size_type n, bool remove_empty = false);
// Resize the Avec
WINAVEC_API Avec avec_resize(const Avec& v, Avec::size_type nbins);
// Resize the Avec by factor
WINAVEC_API Avec rescale_size(const Avec& v, double factor);
Avec avec_map(const Avec& x1, const Avec& y1, const Avec& x2);


// math.h functions
WINAVEC_API Avec acos(const std::vector<dattype> &vec);
WINAVEC_API Avec asin(const std::vector<dattype> &vec);
WINAVEC_API Avec atan(const std::vector<dattype> &vec);
WINAVEC_API Avec atan2(const std::vector<dattype> &y, const std::vector<dattype> &x);
WINAVEC_API Avec ceil(const std::vector<dattype> &vec);
WINAVEC_API Avec cos(const std::vector<dattype> &vec);
WINAVEC_API Avec cosh(const std::vector<dattype> &vec);
WINAVEC_API Avec exp(const std::vector<dattype> &vec);
WINAVEC_API Avec fabs(const std::vector<dattype> &vec);
WINAVEC_API Avec floor(const std::vector<dattype> &vec);
WINAVEC_API Avec fmod(const std::vector<dattype> &x, dattype y);
WINAVEC_API Avec log(const std::vector<dattype> &vec);
WINAVEC_API Avec log10(const std::vector<dattype> &vec);
WINAVEC_API Avec pow(const std::vector<dattype> &vec, dattype y);
WINAVEC_API Avec sin(const std::vector<dattype> &vec);
WINAVEC_API Avec sinh(const std::vector<dattype> &vec);
WINAVEC_API Avec sqrt(const std::vector<dattype> &vec);
WINAVEC_API Avec tan(const std::vector<dattype> &vec);
WINAVEC_API Avec tanh(const std::vector<dattype> &vec);

// Random Numbers
WINAVEC_API Avec rand(std::vector<dattype>::size_type size);
WINAVEC_API Avec rand(const std::vector<dattype> &vec);
WINAVEC_API Avec urand(const std::vector<dattype> &vec);
WINAVEC_API Avec grand(Avec::size_type size, double sigma = 1);
WINAVEC_API Avec grand(const std::vector<dattype> &vec, double sigma=1);

// ASCII Input/Oputput
WINAVEC_API std::ostream& operator<< (std::ostream& o, const std::vector<Avec*>& list);
WINAVEC_API void vecwrite(const std::vector<dattype> &vec, const std::string &filename);
WINAVEC_API void vecwrite(const std::vector<Avec*> &list, const std::string &filename);
WINAVEC_API int vecread(std::vector<dattype> &vec, const std::string &filename);
WINAVEC_API int vecread(const std::vector<Avec*> &list, const std::string &filename);

// Root extensions
#ifdef __CINT__
#define __AVECROOT__
#endif

#ifdef __AVECROOT__
class TArrayD;
class TH1;
class TH2;
class TH1D;
class TGraph;
class TMultiGraph;
class TGraphErrors;
class TFile;

Avec& avec_get(const std::string& name, TFile& file);
Avec& avec_get(const std::string& vecname, const std::string& filename);

// Conversion of Avec to and from histograms and arrays
WINAVEC_API Avec TAD2AV(TArrayD & source);
WINAVEC_API Avec TH12avec(TH1& h1);
WINAVEC_API TH1D* avec2TH1D(Avec& v1, const std::string& histo_title, const std::string histo_name);

// Drawing Avecs

// Helping function to manage ROOT object names
#include "TDirectory.h"
#include "sstream"
template <class ObjectType>
std::string GetValidObjectName(const std::string& name="", int replace=0)
{
  std::string base = (name == "")? "avec" : name;
  
  if(replace != 0){
    // Check if object with this name already exists, if yes remove it
    ObjectType* g_test;
    gDirectory->GetObject(base.c_str(),g_test);
    //   g_test = (ObjectType*)(gDirectory->FindObject( base.c_str() ));
    if( g_test != 0 ){
      g_test->Delete();
      if(Avec::VERBOSE_FLAG) std::cout << "existing Object " << base << " deleted!" << std::endl;
    }
    return base;
  }
  // Check if object with this name already exists, if yes add version number until non-existing name is found
  //ObjectType* g_test;
  std::string try_name=base;
  int ctr=1;
  while ( (ObjectType*)(gDirectory->FindObject( try_name.c_str() )) != 0 && ctr < 100){
    std::ostringstream oss;
    oss << base << "_" << ctr;
    ctr++;
    try_name = oss.str();
  }
  return try_name;
}

// Draw a single Avec
namespace AVEC{
const std::string no_axis_title = "__NO_AXIS_TITLE__";
}
WINAVEC_API TGraph *avec_draw(const Avec& y, const std::string& title="", const std::string & xTitle="", const std::string & yTitle=AVEC::no_axis_title, const std::string& options="AL", int draw=1, const std::string & name="");

// Draw two Avecs against each other
WINAVEC_API TGraph *avec_draw(const Avec& x, const Avec& y, const std::string& title="", const std::string & xTitle=AVEC::no_axis_title, const std::string & yTitle=AVEC::no_axis_title, const std::string& options="AL", int draw=1, const std::string & name="");

// Draw two Avecs against each other, with error bars
WINAVEC_API TGraphErrors *avec_draw(const Avec& x, const Avec& y, const Avec& e_x, const Avec& e_y, const std::string& title="", const std::string & xTitle="", const std::string & yTitle="", const std::string& options="AL", int draw=1, const std::string & name="");


// Draw a set of Avecs against one Avec (const Avec)
WINAVEC_API TMultiGraph *avec_draw(const Avec& x, const std::vector<const Avec*>& y, const std::string& title="", const std::string & xTitle="", const std::string & yTitle="", const std::string& options="AL", int draw=1, const std::string & name="avec_draw_multi");
WINAVEC_API TMultiGraph *avec_draw(const Avec& x, const std::vector<Avec*>& y, const std::string& title="", const std::string & xTitle="", const std::string & yTitle="", const std::string& options="AL", int draw=1, const std::string & name="avec_draw_multi");
WINAVEC_API TMultiGraph *avec_draw(const Avec& x, const std::vector<Avec>& y, const std::string& title="", const std::string & xTitle="", const std::string & yTitle="", const std::string& options="AL", int draw=1, const std::string & name="avec_draw_multi");


// Draw a set of Avecs and a set of error Avecs against one Avec
WINAVEC_API TMultiGraph *avec_draw(const Avec& x, const std::vector<const Avec*>& y, const std::vector<const Avec*>& e_y, const std::string& title="", const std::string & xTitle="", const std::string & yTitle="", const std::string & options="AP", int draw=1, const std::string & name="avec_draw_multi");

// Draw a set of Avecs and a set of error Avecs against one Avec
WINAVEC_API TMultiGraph *avec_draw(const Avec& x, const std::vector<Avec*>& y, const std::vector<Avec*>& e_y, const std::string& title="", const std::string & xTitle="", const std::string & yTitle="", const std::string & options="AP", int draw=1, const std::string & name="avec_draw_multi");

// Draw a set of Avecs against another set of Avecs
WINAVEC_API TMultiGraph *avec_draw(const std::vector<const Avec*>& x, const std::vector<const Avec*>& y, const std::string& title="", const std::string & xTitle="", const std::string & yTitle="", const std::string& options="AL", int draw=1, const std::string & name="avec_draw_multi");

// Draw a set of Avecs against another set of Avecs
WINAVEC_API TMultiGraph *avec_draw(const std::vector<Avec*>& x, const std::vector<Avec*>& y, const std::string& title="", const std::string & xTitle="", const std::string & yTitle="", const std::string& options="AL", int draw=1, const std::string & name="avec_draw_multi");

// Draw a set of Avecs against another set of Avecs
WINAVEC_API TMultiGraph *avec_draw(const std::vector<Avec>& x, const std::vector<Avec>& y, const std::string& title="", const std::string & xTitle="", const std::string & yTitle="", const std::string& options="AL", int draw=1, const std::string & name="avec_draw_multi");

// Draw a set of Avecs and a set of errors against another set of Avecs
WINAVEC_API TMultiGraph *avec_draw(const std::vector<Avec>& x, const std::vector<Avec>& y, const std::vector<Avec>& ey, const std::string& title="", const std::string & xTitle="", const std::string & yTitle="", const std::string& options="AL", int draw=1, const std::string & name="avec_draw_multi");

// Draw a set of Avecs and a set of error Avecs against another set of Avecs
//WINAVEC_API TMultiGraph *avec_draw(const std::vector<const Avec*>& x, const std::vector<const Avec*>& y, const std::vector<const Avec*>& e_y, const std::string& title="", const std::string & xTitle="", const std::string & yTitle="", const std::string & options="AP", int draw=1, const std::string & name="avec_draw");

//WINAVEC_API void AddLegend(TMultiGraph* mgr); 

////////////////////////// Avecs and histograms ///////////////////////////////////////////////////////////////

// Fill the values of one Avec into a histogram weighting with a second Avec
WINAVEC_API void hfill(const Avec& v, TH1& h, const Avec& w=Avec());

// Fill the values of one Avec into a histogram weighting with a second Avec (OBSOLETE)
WINAVEC_API void hfill(const Avec& v, const Avec& w, TH1& h);

// Fill the values of two Avecs into a 2D histogram
WINAVEC_API void hfill(const Avec& x, const Avec& y, TH2& h, const Avec& mask = Avec());


// Fill the values of Avec v for mask!=0 into a histogram and draw it
WINAVEC_API TH1* avec_plot_mask(const Avec& v, const Avec& mask, int nbins=100, const std::string & title="avec_plot", const std::string & name="avec_plot");

// Fill the values of one Avec into a histogram and draw it
WINAVEC_API TH1* avec_plot(const Avec& v, int nbins=100, const std::string & title="avec_plot", const std::string & name="avec_plot");

// Fill the values of two Avecs into a 2D histogram and draw it
WINAVEC_API TH2* avec_plot(const Avec& x, const Avec& y, int nbinsx=100, int nbinsy=100, const std::string & title="avec_plot2D", const std::string & name="avec_plot2D");

// Fill the values of x and y for mask!=0 into a 2D histogram and draw it
WINAVEC_API TH2* avec_plot_mask(const Avec& x, const Avec& y, const Avec& mask, int nbinsx=100, int nbinsy=100, const std::string & title="avec_plot2D", const std::string & name="avec_plot2D");

#else
#error missed root
#endif

#endif
