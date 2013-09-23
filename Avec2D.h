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

#ifndef AVEC2D_H
#define AVEC2D_H

#include "Avec.h"
  //#include <TObject.h>

class Avec2D : public std::vector<Avec>
,  public TNamed
{
 public:
  explicit Avec2D(size_type cols=0, size_type rows=0, const dattype init=0) : std::vector<Avec>(cols,Avec(rows,init)){;}

  Avec2D(const Avec & source, size_type cols);
  Avec2D(size_type cols, const Avec &v) : std::vector<Avec>(cols,v){}
    
    //Avec2D& operator=(const Avec2D&);

  void print(std::ostream &out = std::cout, const char* sep="\n") const;
  friend std::ostream& operator<< (std::ostream& o, const Avec2D& v);

  Avec2D operator+(double offset) const;
  Avec2D operator-(double offset) const;
  Avec2D operator*(double factor) const;
  Avec2D operator/(double factor) const;

  Avec2D operator+(const Avec& source) const;
  Avec2D operator-(const Avec& source) const;
  Avec2D operator*(const Avec& source) const;
  Avec2D operator/(const Avec& source) const;

  Avec2D operator+(const Avec2D& source) const;
  Avec2D operator-(const Avec2D& source) const;
  Avec2D operator*(const Avec2D& source) const;
  Avec2D operator/(const Avec2D& source) const;

  Avec2D operator-() const;

  void operator+=(double offset);
  void operator-=(double offset);
  void operator*=(double factor);
  void operator/=(double factor);

  void operator+=(const Avec& source);
  void operator-=(const Avec& source);
  void operator*=(const Avec& source);
  void operator/=(const Avec& source);

  void operator+=(const Avec2D& source);
  void operator-=(const Avec2D& source);
  void operator*=(const Avec2D& source);
  void operator/=(const Avec2D& source);

  Avec2D operator<(const dattype value) const;
  Avec2D operator>(const dattype value) const;
  Avec2D operator!=(const dattype value) const;

  Avec2D operator&&(const Avec2D& source) const;
  Avec2D operator||(const Avec2D& source) const;

  // Overloaded Methods of TObject
  virtual void	Draw(Option_t* option = "");
  virtual void Browse(TBrowser* b);

  ClassDef(Avec2D,1)
};

// operators
WINAVEC_API Avec2D operator&(const Avec2D& vec1, const Avec2D& vec2);

// Maxima of all columns
Avec vmaxc(const Avec2D& v1);
// Maxima of all rows
Avec vmaxr(const Avec2D& v1);
// Maximum of all elements
dattype vmax(const Avec2D &vec);
// Minima of all columns
Avec vminc(const Avec2D& v1);
// Minima of all rows
Avec vminr(const Avec2D& v1);
// Minimum of all elements
dattype vmin(const Avec2D &vec);

//WINAVEC_API dattype vmin(const Avec2D &vec);
//WINAVEC_API size_t lvmax(const Avec2D &vec);
//WINAVEC_API size_t lvmin(const Avec2D &vec);

// Sum of all components
dattype vsum(const Avec2D& v1);
// Sum all columns
Avec vsumc(const Avec2D& v1);
// Sum all rows
Avec vsumr(const Avec2D& v1);
// Transpose Avec
Avec2D vtrans(const Avec2D& v1);

// Resize the Avec2D in both dimensions 
Avec2D avec_resize(const Avec2D& v, Avec2D::size_type nsizex, Avec::size_type nsizey=0);
// Resize the Avec2D in both dimensions by factorx and factory
Avec2D rescale_size(const Avec2D& v, double factorx, double factory=0);

void avec_remove_nan(Avec2D& data, dattype replacement = 0);


// math.h functions
WINAVEC_API Avec2D acos(const Avec2D &vec);
WINAVEC_API Avec2D asin(const Avec2D &vec);
WINAVEC_API Avec2D atan(const Avec2D &vec);
WINAVEC_API Avec2D atan2(const Avec2D &y, const Avec2D &x);
WINAVEC_API Avec2D ceil(const Avec2D &vec);
WINAVEC_API Avec2D cos(const Avec2D &vec);
WINAVEC_API Avec2D cosh(const Avec2D &vec);
WINAVEC_API Avec2D exp(const Avec2D &vec);
WINAVEC_API Avec2D fabs(const Avec2D &vec);
WINAVEC_API Avec2D floor(const Avec2D &vec);
WINAVEC_API Avec2D fmod(const Avec2D &x, dattype y);
WINAVEC_API Avec2D log(const Avec2D &vec);
WINAVEC_API Avec2D log10(const Avec2D &vec);
WINAVEC_API Avec2D pow(const Avec2D &vec, dattype y);
WINAVEC_API Avec2D sin(const Avec2D &vec);
WINAVEC_API Avec2D sinh(const Avec2D &vec);
WINAVEC_API Avec2D sqrt(const Avec2D &vec);
WINAVEC_API Avec2D tan(const Avec2D &vec);
WINAVEC_API Avec2D tanh(const Avec2D &vec);

// Random Numbers
WINAVEC_API Avec2D rand(const Avec2D &vec);
WINAVEC_API Avec2D urand(const Avec2D &vec);
WINAVEC_API Avec2D grand(const Avec2D &vec, double sigma=1);

// ASCII Input/Oputput
WINAVEC_API void vecwrite(const Avec2D &vec, const std::string &filename);
WINAVEC_API int vecread(Avec2D &vec, const std::string &filename);


// Root extensions
#ifdef __CINT__
#define __AVECROOT__
#endif

#ifdef __AVECROOT__

class TH2D;

Avec2D& avec2d_get(const std::string& name, TFile& file);
Avec2D& avec2d_get(const std::string& vecname, const std::string& filename);

// Draw a set of Avecs against another set of Avecs
//WINAVEC_API TMultiGraph *avec_draw(const std::vector<Avec*>& x, const std::vector<Avec*>& y, const std::string& title="", const std::string & xTitle="", const std::string & yTitle="", const std::string& options="AL", int draw=1, const std::string & name="avec_draw");

// Draw columns of Avec2D against one Avec (must be of same size)
WINAVEC_API TMultiGraph *avec_draw(const Avec& x, const Avec2D& y, const std::string& title="", const std::string & xTitle="", const std::string & yTitle="", const std::string& options="AL", int draw=1, const std::string & name="avec_draw");

// Draw columns of two Avec2D against each other
WINAVEC_API TMultiGraph *avec_draw(const Avec2D& x, const Avec2D& y, const std::string& title="", const std::string & xTitle="", const std::string & yTitle="", const std::string& options="AL", int draw=1, const std::string & name="avec_draw");

// Draw columns of two Avec2D against each other with errors in y
WINAVEC_API TMultiGraph *avec_draw(const Avec2D& x, const Avec2D& y, const Avec2D& err_y, const std::string& title="", const std::string & xTitle="", const std::string & yTitle="", const std::string& options="AL", int draw=1, const std::string & name="avec_draw");

// Draw columns of Avec2D with errors against one Avec (must be of same size)
WINAVEC_API TMultiGraph *avec_draw(const Avec& x, const Avec2D& y, const Avec2D& ey, const std::string& title="", const std::string & xTitle="", const std::string & yTitle="", const std::string& options="AL", int draw=1, const std::string & name="avec_draw");

// Draw a single Avec2D into a 2D histogram, 1 value per bin
WINAVEC_API TH2D *avec_draw(const Avec2D& z, const std::string& title="", const std::string & xTitle="", const std::string & yTitle="", const std::string& options="LEGO2", int draw=1, std::string name="avec_draw");


// Fill all values of one Avec2D into a histogram weighting with a second Avec2D
WINAVEC_API void hfill(const Avec2D& v, TH1& h, const Avec2D& w=Avec2D());

// Fill all values of one Avec2D into a histogram and draw it
WINAVEC_API TH1* avec_plot(const Avec2D& v, int nbins=100, const std::string & title="avec_plot", const std::string & name="avec_plot");


#else
#error missed root
#endif


#endif
